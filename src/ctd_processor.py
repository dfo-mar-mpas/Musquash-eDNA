import os
import pandas as pd
import numpy as np
from datetime import datetime
import glob
import re
from geopy.distance import geodesic
import warnings
warnings.filterwarnings('ignore')

class CTDProcessor:
    def __init__(self, ctd_directory, samples_file):
        self.ctd_directory = ctd_directory
        self.samples_file = samples_file
        self.ctd_data = []
        self.samples_df = None
        
    def parse_ctd_file(self, file_path):
        """Parse a single CTD file and extract metadata and water parameters"""
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            # Initialize metadata dictionary
            metadata = {}
            data_start_line = None
            
            # Parse header metadata
            for i, line in enumerate(lines):
                if line.startswith('%'):
                    if ',' in line:
                        key, value = line[1:].strip().split(',', 1)
                        metadata[key.strip()] = value.strip().split(',')[0].strip()
                else:
                    # Found the data header
                    data_start_line = i
                    break
            
            # Check if coordinates exist and are valid
            start_lat = metadata.get('Start latitude', '')
            start_lon = metadata.get('Start longitude', '')
            
            if not start_lat or not start_lon or start_lat == '' or start_lon == '':
                print(f"Skipping corrupt file (missing coordinates): {file_path}")
                return None
                
            try:
                lat = float(start_lat)
                lon = float(start_lon)
            except (ValueError, TypeError):
                print(f"Skipping corrupt file (invalid coordinates): {file_path}")
                return None
            
            # Read the water parameter data
            if data_start_line is not None:
                # Read the CSV data starting from the header line
                df = pd.read_csv(file_path, skiprows=data_start_line)
                
                # Clean column names
                df.columns = df.columns.str.strip()
                
                if len(df) == 0:
                    print(f"Skipping file with no data: {file_path}")
                    return None
                
                # Extract file info
                file_info = {
                    'file_path': file_path,
                    'file_name': metadata.get('File name', os.path.basename(file_path)),
                    'device': metadata.get('Device', ''),
                    'cast_time_utc': metadata.get('Cast time (UTC)', ''),
                    'cast_time_local': metadata.get('Cast time (local)', ''),
                    'latitude': lat,
                    'longitude': lon,
                    'cast_duration': metadata.get('Cast duration (Seconds)', ''),
                    'data': df
                }
                
                return file_info
            else:
                print(f"No data section found in file: {file_path}")
                return None
                
        except Exception as e:
            print(f"Error processing file {file_path}: {str(e)}")
            return None
    
    def load_ctd_files(self):
        """Load all valid CTD files from directory"""
        print("Loading CTD files...")
        
        # Find all CTD files (assuming they have specific extensions or patterns)
        ctd_files = []
        for ext in ['*.txt', '*.ctd', '*.csv']:
            ctd_files.extend(glob.glob(os.path.join(self.ctd_directory, ext)))
        
        if not ctd_files:
            # Try to find files without specific extensions
            ctd_files = [f for f in glob.glob(os.path.join(self.ctd_directory, '*')) 
                        if os.path.isfile(f)]
        
        print(f"Found {len(ctd_files)} potential CTD files")
        
        valid_files = 0
        for file_path in ctd_files:
            ctd_info = self.parse_ctd_file(file_path)
            if ctd_info is not None:
                self.ctd_data.append(ctd_info)
                valid_files += 1
        
        print(f"Successfully loaded {valid_files} valid CTD files")
        return valid_files > 0
    
    def load_samples_file(self):
        """Load the samples metadata file"""
        try:
            # Try different separators
            for sep in [' ', '\t', ',']:
                try:
                    df = pd.read_csv(self.samples_file, sep=sep)
                    if len(df.columns) >= 6:  # Should have at least 6 columns
                        self.samples_df = df
                        break
                except:
                    continue
            
            if self.samples_df is None:
                raise Exception("Could not parse samples file with any separator")
            
            # Clean column names
            self.samples_df.columns = self.samples_df.columns.str.strip()
            
            # Ensure we have the required columns
            required_cols = ['Sample_ID', 'File_Name', 'Water_Sample', 'Station', 'Sample_Type', 'Latitude', 'Longitude']
            
            if not all(col in self.samples_df.columns for col in required_cols):
                print("Warning: Some expected columns missing. Available columns:")
                print(self.samples_df.columns.tolist())
            
            print(f"Loaded {len(self.samples_df)} samples")
            return True
            
        except Exception as e:
            print(f"Error loading samples file: {str(e)}")
            return False
    
    def calculate_distance(self, lat1, lon1, lat2, lon2):
        """Calculate distance between two coordinates using geodesic distance"""
        try:
            return geodesic((lat1, lon1), (lat2, lon2)).meters
        except:
            # Fallback to simple Euclidean distance if geodesic fails
            return np.sqrt((lat1 - lat2)**2 + (lon1 - lon2)**2) * 111000  # rough conversion to meters
    
    def find_nearest_ctd(self, sample_lat, sample_lon, max_distance_m=5000):
        """Find the nearest CTD measurement to a sample location"""
        min_distance = float('inf')
        nearest_ctd = None
        
        for ctd in self.ctd_data:
            distance = self.calculate_distance(sample_lat, sample_lon, 
                                             ctd['latitude'], ctd['longitude'])
            
            if distance < min_distance and distance <= max_distance_m:
                min_distance = distance
                nearest_ctd = ctd
        
        return nearest_ctd, min_distance
    
    def get_water_parameters_summary(self, ctd_data, depth_range=(0, 5)):
        """Extract water parameters from CTD data, averaging over specified depth range"""
        df = ctd_data['data']
        
        # Filter data by depth range
        if 'Depth (Meter)' in df.columns:
            mask = (df['Depth (Meter)'] >= depth_range[0]) & (df['Depth (Meter)'] <= depth_range[1])
            filtered_df = df[mask]
        else:
            # Use pressure as proxy for depth if depth column not available
            if 'Pressure (Decibar)' in df.columns:
                mask = (df['Pressure (Decibar)'] >= depth_range[0]) & (df['Pressure (Decibar)'] <= depth_range[1])
                filtered_df = df[mask]
            else:
                # Use first few rows if no depth/pressure info
                filtered_df = df.head(5)
        
        if len(filtered_df) == 0:
            filtered_df = df.head(1)  # Use first measurement if no data in range
        
        # Calculate mean values for water parameters
        parameters = {}
        
        param_columns = {
            'Temperature (Celsius)': 'Temperature_C',
            'Conductivity (MicroSiemens per Centimeter)': 'Conductivity_uS_cm',
            'Specific conductance (MicroSiemens per Centimeter)': 'Specific_Conductance_uS_cm',
            'Salinity (Practical Salinity Scale)': 'Salinity_PSS',
            'Sound velocity (Meters per Second)': 'Sound_Velocity_m_s',
            'Density (Kilograms per Cubic Meter)': 'Density_kg_m3',
            'Pressure (Decibar)': 'Pressure_dbar',
            'Depth (Meter)': 'Depth_m'
        }
        
        for ctd_col, param_name in param_columns.items():
            if ctd_col in filtered_df.columns:
                parameters[param_name] = filtered_df[ctd_col].mean()
            else:
                parameters[param_name] = np.nan
        
        # Add metadata
        parameters['CTD_File'] = ctd_data['file_name']
        parameters['CTD_Latitude'] = ctd_data['latitude']
        parameters['CTD_Longitude'] = ctd_data['longitude']
        parameters['CTD_Cast_Time'] = ctd_data['cast_time_utc']
        parameters['Measurements_Count'] = len(filtered_df)
        parameters['Depth_Range_Used'] = f"{depth_range[0]}-{depth_range[1]}m"
        
        return parameters
    
    def process_samples(self, max_distance_m=5000, depth_range=(0, 5)):
        """Process all samples and match with nearest CTD data"""
        if self.samples_df is None or len(self.ctd_data) == 0:
            print("Error: Samples file or CTD data not loaded")
            return None
        
        print(f"Processing {len(self.samples_df)} samples...")
        print(f"Maximum search distance: {max_distance_m/1000:.1f} km")
        print(f"Depth range for water parameters: {depth_range[0]}-{depth_range[1]} m")
        
        results = []
        matches_found = 0
        
        for idx, sample in self.samples_df.iterrows():
            try:
                sample_lat = float(sample['Latitude'])
                sample_lon = float(sample['Longitude'])
                
                # Find nearest CTD
                nearest_ctd, distance = self.find_nearest_ctd(sample_lat, sample_lon, max_distance_m)
                
                if nearest_ctd is not None:
                    # Get water parameters
                    water_params = self.get_water_parameters_summary(nearest_ctd, depth_range)
                    water_params['Distance_to_CTD_m'] = distance
                    matches_found += 1
                    
                    print(f"Sample {sample['Sample_ID']}: Found CTD match at {distance:.0f}m distance")
                else:
                    # No match found
                    water_params = {param: np.nan for param in [
                        'Temperature_C', 'Conductivity_uS_cm', 'Specific_Conductance_uS_cm',
                        'Salinity_PSS', 'Sound_Velocity_m_s', 'Density_kg_m3', 
                        'Pressure_dbar', 'Depth_m'
                    ]}
                    water_params.update({
                        'CTD_File': 'No match found',
                        'CTD_Latitude': np.nan,
                        'CTD_Longitude': np.nan,
                        'CTD_Cast_Time': 'No match found',
                        'Distance_to_CTD_m': np.nan,
                        'Measurements_Count': 0,
                        'Depth_Range_Used': f"{depth_range[0]}-{depth_range[1]}m"
                    })
                    print(f"Sample {sample['Sample_ID']}: No CTD match found within {max_distance_m/1000:.1f}km")
                
                # Combine original sample data with water parameters
                result = sample.to_dict()
                result.update(water_params)
                results.append(result)
                
            except Exception as e:
                print(f"Error processing sample {sample.get('Sample_ID', idx)}: {str(e)}")
                continue
        
        print(f"\nProcessing complete: {matches_found}/{len(self.samples_df)} samples matched with CTD data")
        
        # Create results DataFrame
        results_df = pd.DataFrame(results)
        return results_df
    
    def save_results(self, results_df, output_file):
        """Save results to CSV file"""
        try:
            results_df.to_csv(output_file, index=False)
            print(f"Results saved to: {output_file}")
            
            # Print summary statistics
            print("\n=== SUMMARY ===")
            total_samples = len(results_df)
            matched_samples = results_df['Distance_to_CTD_m'].notna().sum()
            
            print(f"Total samples: {total_samples}")
            print(f"Matched with CTD: {matched_samples}")
            print(f"Match rate: {matched_samples/total_samples*100:.1f}%")
            
            if matched_samples > 0:
                avg_distance = results_df['Distance_to_CTD_m'].mean()
                max_distance = results_df['Distance_to_CTD_m'].max()
                print(f"Average distance to CTD: {avg_distance:.0f}m")
                print(f"Maximum distance to CTD: {max_distance:.0f}m")
                
                # Show parameter availability
                param_cols = ['Temperature_C', 'Salinity_PSS', 'Conductivity_uS_cm']
                for col in param_cols:
                    if col in results_df.columns:
                        available = results_df[col].notna().sum()
                        print(f"{col}: {available}/{total_samples} samples ({available/total_samples*100:.1f}%)")
            
            return True
        except Exception as e:
            print(f"Error saving results: {str(e)}")
            return False

def main():
    """Main execution function"""
    # Configuration
    CTD_DIRECTORY = "/mnt/e/projects/Musquash_data/Metadata/Castaway-CTD" 
    SAMPLES_FILE = "/mnt/e/projects/Musquash_data/Metadata/all_years_metada.tsv" 
    OUTPUT_FILE = "/mnt/e/projects/Musquash_data/Metadata/Musq_CTD_water_parameters.csv"
    MAX_DISTANCE_M = 200  # Maximum distance to search for CTD data (meters)
    DEPTH_RANGE = (0.8, 1.2)   # Depth range to average water parameters (meters)
    
    print("=== CTD Water Parameters Processor ===")
    print(f"CTD Directory: {CTD_DIRECTORY}")
    print(f"Samples File: {SAMPLES_FILE}")
    print(f"Output File: {OUTPUT_FILE}")
    print()
    
    # Initialize processor
    processor = CTDProcessor(CTD_DIRECTORY, SAMPLES_FILE)
    
    # Load CTD files
    if not processor.load_ctd_files():
        print("Error: No valid CTD files found!")
        return
    
    # Load samples file
    if not processor.load_samples_file():
        print("Error: Could not load samples file!")
        return
    
    # Process samples and match with CTD data
    results = processor.process_samples(max_distance_m=MAX_DISTANCE_M, 
                                       depth_range=DEPTH_RANGE)
    
    if results is not None:
        # Save results
        processor.save_results(results, OUTPUT_FILE)
    else:
        print("Error: Processing failed!")

if __name__ == "__main__":
    main()