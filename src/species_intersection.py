import pandas as pd
import numpy as np
import re
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from upsetplot import UpSet, from_memberships
from typing import Optional, Tuple
from taxonomy_utils import NCBITaxonomyLookup
# import importlib
# importlib.reload(taxonomy_utils)


# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION  ── edit these paths / labels to match your environment
# ─────────────────────────────────────────────────────────────────────────────

ROOT       = "/mnt/e/projects/Musquash_data/Analysis"
INPUT_DIR  = f"{ROOT}/data2"
OUTPUT_DIR = f"{INPUT_DIR}/species_viz"
MARKER     = "mifish"

# eDNA QIIME2 collapsed-species tables  →  integer year label
EDNA_FILES = {
    f"{INPUT_DIR}/Musq-25-{MARKER}-filt-grpd-collapsed-species-table.tsv": 2025,
    f"{INPUT_DIR}/Musq-24-{MARKER}-filt-grpd-collapsed-species-table.tsv": 2024,
    f"{INPUT_DIR}/Musq-23-{MARKER}-filt-grpd-collapsed-species-table.tsv": 2023,
    f"{INPUT_DIR}/Musq-22-{MARKER}-filt-grpd-collapsed-species-table.tsv": 2022,
}

# Historical-record species-name files  →  study label string
HIST_FILES = {
    f"{ROOT}/Expected_species/Andrew_2021.txt"    : "Cooper et al. 2021",
    f"{ROOT}/Expected_species/Ipsen_2013.txt"     : "Ipsen et al. 2013",
    f"{ROOT}/Expected_species/Arens_2007.txt"     : "Arens et al. 2007",
    f"{ROOT}/Expected_species/singh_fish_2000.txt": "Singh et al. 2000",
}

STANDARD_ZONES = ["Zone_1", "Zone_2A", "Zone_2B", "Zone_3"]


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _extract_species(taxonomy: str) -> str:
    """Return binomial name from a full QIIME2 taxonomy string, or ''."""
    if pd.isna(taxonomy):
        return ""
    m = re.search(r"s__([A-Za-z]+\s+[A-Za-z]+)", str(taxonomy))
    return m.group(1).strip() if m else ""


def _year_sort_key(val):
    """
    Splits the string into an integer year and a string suffix.
    Returns a tuple: (year, suffix) for hierarchical sorting.
    """
    parts = str(val).split(' ')
    
    # 1. Extract the year and convert it to an integer
    year = int(parts[0])
    
    # 2. Extract the suffix if it exists, otherwise use an empty string.
    # Lowercasing ensures 'M' and 'm' sort predictably.
    suffix = "".join(parts[1:]).lower() if len(parts) > 1 else ""
    
    return (year, suffix)

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — load eDNA tables and outer-merge across years, keyed on Species
# ─────────────────────────────────────────────────────────────────────────────

def load_edna_tables(edna_files: dict) -> pd.DataFrame:
    """
    Read each QIIME2 TSV (header on row index 1), extract species names from
    the taxonomy string, prefix every zone column with the file's year label,
    and outer-merge all years on Species.

    Returns
    -------
    pd.DataFrame
        Columns: Taxonomy, Species, Zone_1_<yr>, Zone_2A_<yr>, … for every yr
    """
    frames = []

    for filepath, year in edna_files.items():
        df = pd.read_csv(filepath, sep="\t", header=1)
        df = df.rename(columns={"#OTU ID": "Taxonomy"})
        df["Species"] = df["Taxonomy"].apply(_extract_species)
        df = df[df["Species"] != ""].copy()

        zone_cols = [c for c in df.columns if c not in ("Taxonomy", "Species")]
        df = df.rename(columns={c: f"{c}_{year}" for c in zone_cols})
        frames.append(df)

    combined = frames[0]
    for df in frames[1:]:
        # Build a lookup so we can back-fill Taxonomy for new species
        tax_lookup = (
            df[["Species", "Taxonomy"]]
            .query("Taxonomy != ''")
            .drop_duplicates("Species")
            .set_index("Species")["Taxonomy"]
        )
        combined = pd.merge(combined, df.drop(columns="Taxonomy"),
                            on="Species", how="outer")

        missing = combined["Taxonomy"].isna() | (combined["Taxonomy"] == "")
        new_tax_mask = missing & combined["Species"].isin(tax_lookup.index)
        combined.loc[new_tax_mask, "Taxonomy"] = (
            combined.loc[new_tax_mask, "Species"].map(tax_lookup)
        )

    return combined


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — load historical presence files
# ─────────────────────────────────────────────────────────────────────────────

def load_hist_files(hist_files: dict) -> pd.DataFrame:
    """
    Load plain-text (or delimited) historical species files.
    Each file contributes one binary column (0/1) labelled with its study string.

    Returns
    -------
    pd.DataFrame
        Columns: Species, <study_label1>, <study_label2>, …
    """
    presence: dict[str, set] = {}

    for filepath, label in hist_files.items():
        try:
            raw = pd.read_csv(filepath, sep=None, engine="python", header=0)
        except Exception:
            raw = pd.read_csv(filepath, header=0)

        species = (
            raw.iloc[:, 0]
            .astype(str).str.strip()
            .pipe(lambda s: s[s.notna() & (s != "") & (s != "nan")])
            .tolist()
        )
        presence[label] = set(species)

    all_species = sorted(set().union(*presence.values()))
    rows = {"Species": all_species}
    for label, sp_set in presence.items():
        rows[label] = [1 if sp in sp_set else 0 for sp in all_species]

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — build master wide table (zones × years + historical columns)
# ─────────────────────────────────────────────────────────────────────────────

def build_master_table(
    edna_files: dict,
    hist_files: Optional[dict] = None,
    standard_zones: list = STANDARD_ZONES,
) -> pd.DataFrame:
    """
    Outer-join the eDNA wide table with optional historical presence columns.
    Column order: Taxonomy, Species, zone-year blocks (ascending year),
    then historical study columns.
    Rows sorted by descending total count; NaN → 0.
    """
    edna_df   = load_edna_tables(edna_files)
    # all_years = sorted(set(edna_files.values()), key=_year_sort_key)
    all_years = edna_files.values()

    if hist_files:
        hist_df = load_hist_files(hist_files)
        master  = pd.merge(edna_df, hist_df, on="Species", how="outer")
    else:
        master = edna_df.copy()

    master["Taxonomy"] = master["Taxonomy"].fillna("")

    # Coerce all non-identifier columns to integer
    for col in [c for c in master.columns if c not in ("Taxonomy", "Species")]:
        master[col] = pd.to_numeric(master[col], errors="coerce").fillna(0).astype(int)

    # ── column ordering ──────────────────────────────────────────────────────
    ordered = ["Taxonomy", "Species"]
    for year in all_years:
        for zone in standard_zones:
            col = f"{zone}_{year}"
            if col in master.columns:
                ordered.append(col)

    for lbl in (list(hist_files.values()) if hist_files else []):
        if lbl in master.columns:
            ordered.append(lbl)

    # Append any unexpected extra columns at the end
    ordered += [c for c in master.columns if c not in ordered]
    master = master[[c for c in ordered if c in master.columns]]

    # ── row ordering: descending total ───────────────────────────────────────
    count_cols   = [c for c in master.columns if c not in ("Taxonomy", "Species")]
    master["_t"] = master[count_cols].sum(axis=1)
    master = (master.sort_values("_t", ascending=False)
                    .drop(columns="_t")
                    .reset_index(drop=True))

    return master


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — collapse zone columns into per-year totals
# ─────────────────────────────────────────────────────────────────────────────

def aggregate_to_year_abundance(
    master: pd.DataFrame,
    edna_files: dict,
    hist_files: Optional[dict] = None,
) -> pd.DataFrame:
    # all_years = sorted(set(edna_files.values()), key=_year_sort_key)
    all_years = list(edna_files.values())if edna_files else []
    hist_labels = list(hist_files.values()) if hist_files else []
    
    edna_strs = [str(y) for y in all_years]
    all_source_cols = edna_strs + hist_labels

    rows: dict = {"Species": master["Species"].values}

    for year in all_years:
        zone_cols = [c for c in master.columns if c.endswith(f"_{year}")]
        rows[str(year)] = master[zone_cols].sum(axis=1).values if zone_cols else 0

    for lbl in hist_labels:
        if lbl in master.columns:
            rows[lbl] = master[lbl].values

    agg = pd.DataFrame(rows)

    # Calculate binary presence (1 or 0) for sorting
    presence_binary = agg[all_source_cols].map(lambda x: 1 if x > 0 else 0)
    
    # 1. Total number of sources (Years + Studies)
    agg["_nz"] = presence_binary.sum(axis=1)
    
    # 2. Unique pattern string (e.g., "1-1-0-0-1...")
    agg["_pattern"] = presence_binary.apply(lambda row: "".join(row.values.astype(str)), axis=1)
    
    # 3. Total abundance for tie-breaking
    agg["_tot"] = agg[edna_strs].sum(axis=1)

    # SORTING: 
    # _nz (Desc): Most hits at the top
    # _pattern (Desc): Groups identical patterns; 1s before 0s (staircase effect)
    # _tot (Desc): Highest reads within the same group
    agg = (agg.sort_values(["_nz", "_pattern", "_tot"], ascending=[False, False, False])
              .drop(columns=["_nz", "_pattern", "_tot"])
              .reset_index(drop=True))

    return agg


def sort_species_dataframe(
    df: pd.DataFrame,
    source_cols: list[str],
    abundance_cols: list[str],
) -> pd.DataFrame:
    """
    Sort a species DataFrame by detection pattern for heatmap staircase effect.

    Args:
        df:             DataFrame with a 'Species' column + numeric detection columns.
        source_cols:    All columns to consider for presence/absence (eDNA + historical).
        abundance_cols: Columns used for total abundance tie-breaking (eDNA only).

    Returns:
        Sorted DataFrame with temporary sort columns removed.
    """
    presence_binary = df[source_cols].map(lambda x: 1 if x > 0 else 0)

    df = df.copy()
    df["_nz"]      = presence_binary.sum(axis=1)
    df["_pattern"] = presence_binary.apply(lambda row: "".join(row.values.astype(str)), axis=1)
    df["_tot"]     = df[abundance_cols].sum(axis=1)

    return (
        df.sort_values(["_nz", "_pattern", "_tot"], ascending=[False, False, False])
          .drop(columns=["_nz", "_pattern", "_tot"])
          .reset_index(drop=True)
    )
# ─────────────────────────────────────────────────────────────────────────────
# STEP 5 — UpSet plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_upset(
    agg_df: pd.DataFrame,
    edna_cols: list[str],
    hist_cols: list[str],
    output_dir: str,
    name_prefix: str = "Musq",
    sort_by: str = "-degree",
    sort_categories_by: str = "input",
) -> None:
    """
    UpSet plot of species set intersections.
    A species is 'present' in an eDNA column if its count > 0;
    historical columns already hold 0/1 presence values.
    """
    os.makedirs(output_dir, exist_ok=True)

    all_cols = edna_cols + hist_cols
    presence = agg_df[all_cols].copy()
    for col in edna_cols:
        presence[col] = (presence[col] > 0).astype(int)

    memberships = [
        tuple(col for col in all_cols if row[col] == 1)
        for _, row in presence.iterrows()
    ]
    memberships = [m for m in memberships if m]  # drop empty rows

    if not memberships:
        print("UpSet: no memberships — skipping.")
        return

    fig = plt.figure(figsize=(14, 7))
    UpSet(
        from_memberships(memberships),
        subset_size="count",
        show_counts=True,
        sort_by=sort_by,
        sort_categories_by=sort_categories_by,
        min_subset_size=1,
        element_size=40,
    ).plot(fig=fig)

    plt.suptitle("Species Set Intersections", fontsize=14, fontweight="bold", y=0.98)
    out = f"{output_dir}/{name_prefix}_upset.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 6 — vertical heatmap (species on Y, surveys on X)
# ─────────────────────────────────────────────────────────────────────────────

def plot_COI_heatmap_vertical(
    agg_df: pd.DataFrame,
    edna_cols: list[str],
    hist_cols: list[str],
    output_dir: str,
    taxonomy_lookup: Optional[NCBITaxonomyLookup] = None,
    name_prefix: str = "Musq",
    edna_abundant_color: str = "#2E86AB",
    edna_rare_color: str = "#C8C8C8",
    hist_color: str = "#FF6B35",
    cell_size: float = 0.30,
    abundance_threshold: float = 0.01,
    max_rows_per_panel: Optional[int] = None,
    class_order: Optional[list[str]] = None,
    class_bar_width: float = 0.30,             # ◄ CHANGED: doubled from 0.18
) -> None:
    import math
    import itertools
    from collections import Counter  # ◄ CHANGED

    _DEFAULT_CLASS_ORDER = [
        "Actinopterygii", "Echinoidea", "Asteroidea", "Enteropneusta",
        "Malacostraca", "Thecostraca", "Maxillopoda", "Gastropoda",
        "Bivalvia", "Polychaeta", "Clitellata", "Gymnolaemata",
        "Pilidiophora", "Enoplea", "Tentaculata", "Hydrozoa", "Demospongiae",
    ]
    active_class_order = class_order if class_order is not None else _DEFAULT_CLASS_ORDER

    os.makedirs(output_dir, exist_ok=True)

    agg_df   = sort_species_dataframe(agg_df, edna_cols + hist_cols, edna_cols)
    all_cols = edna_cols + hist_cols
    n_cols   = len(all_cols)

    # ── CLASS LOOKUP + RE-SORT ────────────────────────────────────────────────
    species_to_class: dict[str, str] = {}
    show_class_bars = False

    if taxonomy_lookup is not None:
        try:
            sp_list = list(agg_df["Species"].values)
            raw = taxonomy_lookup.get_ranks_batch(sp_list, rank="class")
            species_to_class = {sp: (cls or "Unknown") for sp, cls in raw.items()}
            show_class_bars = True
        except Exception as exc:
            print(f"Warning: class lookup failed – class bars disabled. ({exc})")

    if show_class_bars:
        def _class_rank(i: int) -> int:
            cls = species_to_class.get(agg_df["Species"].iloc[i], "Unknown")
            try:
                return active_class_order.index(cls)
            except ValueError:
                return len(active_class_order)

        sorted_idx = sorted(range(len(agg_df)), key=_class_rank)
        agg_df = agg_df.iloc[sorted_idx].reset_index(drop=True)

    species   = agg_df["Species"].values
    n_species = len(species)

    # ── CLASS COLOUR MAP ──────────────────────────────────────────────────────
    _CLASS_PALETTE = [
        "#D6EAF8", "#D5F5E3", "#FDEBD0", "#E8DAEF",
        "#FADBD8", "#D1F2EB", "#FEF9E7", "#EAF2FF",
        "#F9EBEA", "#E9F7EF", "#FDF2F8", "#EBF5FB",
        "#F4ECF7", "#E8F8F5", "#FEF5E7", "#EAFAF1",
        "#F5EEF8",
    ]
    _seen_classes: list[str] = list(dict.fromkeys(
        species_to_class.get(sp, "Unknown") for sp in species
    ))
    class_color_map = {
        cls: _CLASS_PALETTE[i % len(_CLASS_PALETTE)]
        for i, cls in enumerate(_seen_classes)
    }

    # ── HELPER: compute (class_name, row_start, row_end) spans ───────────────
    def _class_spans(panel_species) -> list[tuple[str, int, int]]:
        spans: list[tuple[str, int, int]] = []
        for cls, grp in itertools.groupby(
            panel_species, key=lambda sp: species_to_class.get(sp, "Unknown")
        ):
            count = sum(1 for _ in grp)
            start = spans[-1][2] if spans else 0
            spans.append((cls, start, start + count))
        return spans

    # ── PANEL SPLIT ───────────────────────────────────────────────────────────
    use_two_panels = max_rows_per_panel is not None and n_species > max_rows_per_panel
    if use_two_panels:
        split          = math.ceil(n_species / 2)
        species_chunks = [species[:split], species[split:]]
    else:
        species_chunks = [species]

    # ── SCALING ───────────────────────────────────────────────────────────────
    scale        = cell_size / 0.30
    fs_ylabel    = max(5,  13 * scale)
    fs_xlabel    = max(5, 14 * scale)
    fs_axlabel   = max(6, 13 * scale)
    fs_legend    = max(6, 12 * scale)
    fs_classbar  = max(4,  11 * scale)
    marker_size  = max(6, 12 * scale)

    # ── RELATIVE ABUNDANCE ────────────────────────────────────────────────────
    rel_abundance: dict[str, dict] = {}
    for col in edna_cols:
        col_total = agg_df[col].sum()
        rel_abundance[col] = (
            (agg_df.set_index("Species")[col] / col_total).to_dict()
            if col_total > 0
            else {sp: 0.0 for sp in species}
        )

    def _cell_color(col: str, sp: str, val: float) -> str:
        if val <= 0:
            return "white"
        if col in edna_cols:
            return (
                edna_abundant_color
                if rel_abundance[col].get(sp, 0.0) >= abundance_threshold
                else edna_rare_color
            )
        return hist_color

    # ── COMMON NAMES (with uniqueness filter) ─────────────────────────────────
    common_names: dict[str, str] = {}
    if taxonomy_lookup is not None:
        try:
            raw_common = taxonomy_lookup.get_common_names_batch(
                list(species), hierarchical=True
            )
            # ◄ CHANGED: suppress common names shared by more than one species
            name_counts = Counter(v for v in raw_common.values() if v)
            common_names = {
                sp: name for sp, name in raw_common.items()
                if name and name_counts[name] == 1
            }
        except Exception as exc:
            print(f"Warning: batch common-name lookup failed: {exc}")

    def _make_ylabel(sp: str, common) -> str:
        italic = r"$\it{" + sp.replace(" ", r"\ ").replace("_", r"\ ") + r"}$"
        return italic + f"  ({common})" if common else italic

    # ── GEOMETRY (all in inches) ──────────────────────────────────────────────
    CB_W          = class_bar_width if show_class_bars else 0.0
    CB_GAP        = 0.06 if show_class_bars else 0.0

    # ◄ CHANGED: increased labelpad to prevent y-axis labels overlapping the class bar
    YLABEL_LABELPAD  = CB_W * 72 + 4 if show_class_bars else 4  # points (72 pt/inch)

    LEFT_PAD         = 0.15 + fs_ylabel * 0.10
    TOP_PAD          = 0.15 + fs_xlabel * 0.14
    BOT_AXIS_PAD     = 0.20
    LEGEND_PAD       = 0.55
    PANEL_GAP        = 0.50
    RIGHT_YLABEL_PAD = 0.10 + fs_ylabel * 0.10

    axes_w  = n_cols * cell_size
    axes_h1 = len(species_chunks[0]) * cell_size
    axes_h2 = len(species_chunks[1]) * cell_size if use_two_panels else 0.0
    axes_h  = axes_h1

    fig_h = LEGEND_PAD + BOT_AXIS_PAD + axes_h + TOP_PAD

    hmap_left1 = LEFT_PAD + CB_W + CB_GAP
    if use_two_panels:
        hmap_left2 = hmap_left1 + axes_w + PANEL_GAP + CB_W + CB_GAP
        fig_w = hmap_left2 + axes_w + RIGHT_YLABEL_PAD
    else:
        fig_w = hmap_left1 + axes_w + RIGHT_YLABEL_PAD

    plt.rcParams['svg.fonttype'] = 'path'
    fig = plt.figure(figsize=(fig_w, fig_h))

    top_inches = LEGEND_PAD + BOT_AXIS_PAD + axes_h

    def _add_axes(left_in: float, h_in: float, w_in: Optional[float] = None):
        if w_in is None:
            w_in = axes_w
        bot_in = top_inches - h_in
        return fig.add_axes([
            left_in / fig_w,
            bot_in  / fig_h,
            w_in    / fig_w,
            h_in    / fig_h,
        ])

    # ── BUILD PANEL LIST ──────────────────────────────────────────────────────
    ax0 = _add_axes(hmap_left1, axes_h1)
    panel_info = [(ax0, species_chunks[0], "left", hmap_left1, axes_h1)]

    if use_two_panels:
        ax1 = _add_axes(hmap_left2, axes_h2)
        panel_info.append((ax1, species_chunks[1], "right", hmap_left2, axes_h2))

    # ── DRAW CLASS BARS ───────────────────────────────────────────────────────
    def _draw_class_bar(cb_left: float, h_in: float, panel_species) -> None:
        spans = _class_spans(panel_species)
        if not spans:
            return
        n_rows = len(panel_species)
        cax = _add_axes(cb_left, h_in, w_in=CB_W)
        cax.set_xlim(0, 1)
        cax.set_ylim(n_rows, 0)
        cax.set_aspect("auto")
        cax.axis("off")

        for cls_name, row_start, row_end in spans:
            color = class_color_map.get(cls_name, "#F0F0F0")
            span_h = row_end - row_start

            # ◄ CHANGED: inset scaled to wider bar (keep ~5% margins on each side)
            cax.add_patch(Rectangle(
                (0.05, row_start + 0.04), 0.90, span_h - 0.08,
                facecolor=color, edgecolor="#888888",
                linewidth=0.6, clip_on=False,
            ))

            if span_h * cell_size > 0.18:
                # ◄ CHANGED: abbreviate label to first 3 letters + "." when span < 4 rows
                label = (
                    cls_name[:3] + "."
                    if span_h < 4
                    else cls_name
                )
                cax.text(
                    0.50, (row_start + row_end) / 2,
                    label,
                    ha="center", va="center",
                    rotation=90,
                    fontsize=fs_classbar, fontweight="bold",
                    color="#333333", clip_on=False,
                )

    if show_class_bars:
        for _ax, chunk, _side, hmap_left, h_in in panel_info:
            cb_left = hmap_left - CB_GAP - CB_W
            _draw_class_bar(cb_left, h_in, chunk)

    # ── DRAW HEATMAP PANELS ───────────────────────────────────────────────────
    for ax, panel_species, ylabel_side, _, _ in panel_info:
        n_rows = len(panel_species)

        for i, sp in enumerate(panel_species):
            row = agg_df.loc[agg_df["Species"] == sp, all_cols]
            for j, col in enumerate(all_cols):
                val = float(row[col].values[0]) if len(row) else 0.0
                ax.add_patch(Rectangle(
                    (j, i), 1, 1,
                    facecolor=_cell_color(col, sp, val),
                    edgecolor="lightgray", linewidth=0.5,
                ))

        ax.set_xlim(0, n_cols)
        ax.set_ylim(n_rows, 0)
        ax.set_aspect("equal", adjustable="box")

        y_labels = [_make_ylabel(sp, common_names.get(sp)) for sp in panel_species]
        ax.set_yticks(np.arange(n_rows) + 0.5)

        if ylabel_side == "right":
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")
            ax.set_yticklabels(y_labels, fontsize=fs_ylabel)
        else:
            # ◄ CHANGED: pad y-tick labels so they clear the class bar
            ax.set_yticklabels(y_labels, fontsize=fs_ylabel)
            ax.yaxis.set_tick_params(pad=YLABEL_LABELPAD)

        ax.set_xticks(np.arange(n_cols) + 0.5)
        ax.set_xticklabels(
            all_cols, rotation=45, ha="left",
            fontsize=fs_xlabel, fontweight="bold",
        )
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")

        for spine in ax.spines.values():
            spine.set_visible(False)

    # ── LEGEND ────────────────────────────────────────────────────────────────
    legend_handles = [
        Line2D([0], [0], marker="s", color="w",
               markerfacecolor=edna_abundant_color, markersize=marker_size,
               label=f"eDNA Detected (≥{abundance_threshold:.1%} rel. abund.)"),
        Line2D([0], [0], marker="s", color="w",
               markerfacecolor=edna_rare_color, markersize=marker_size,
               label=f"eDNA Detected (<{abundance_threshold:.1%} rel. abund.)"),
        Line2D([0], [0], marker="s", color="w",
               markerfacecolor=hist_color, markersize=marker_size,
               label="Historical Record"),
    ]
    # fig.legend(
    #     handles=legend_handles,
    #     loc="lower center",
    #     bbox_to_anchor=(0.5, 0.0),
    #     bbox_transform=fig.transFigure,
    #     ncol=3, frameon=True,
    #     fontsize=fs_legend, borderaxespad=0.5,
    # )

    png_path = f"{output_dir}/{name_prefix}_heatmap_vertical.png"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(png_path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png_path}")


def plot_heatmap_vertical(
    agg_df: pd.DataFrame,
    edna_cols: list[str],
    hist_cols: list[str],
    output_dir: str,
    taxonomy_lookup: Optional[NCBITaxonomyLookup] = None,
    name_prefix: str = "Musq",
    edna_abundant_color: str = "#2E86AB",
    edna_rare_color: str = "#C8C8C8",
    hist_color: str = "#FF6B35",
    cell_size: float = 0.30,
    abundance_threshold: float = 0.01,
) -> None:
    from collections import Counter  # ◄ CHANGED

    os.makedirs(output_dir, exist_ok=True)
    agg_df = sort_species_dataframe(agg_df, edna_cols + hist_cols, edna_cols)
    all_cols  = edna_cols + hist_cols
    species   = agg_df["Species"].values
    n_species = len(species)
    n_cols    = len(all_cols)

    scale        = cell_size / 0.30
    fs_ylabel    = max(5, 12  * scale)
    fs_xlabel    = max(5, 14 * scale)
    fs_axlabel   = max(6, 13 * scale)
    fs_legend    = max(6, 12 * scale)
    marker_size  = max(6, 12 * scale)

    rel_abundance: dict[str, dict] = {}
    for col in edna_cols:
        col_total = agg_df[col].sum()
        if col_total > 0:
            rel_abundance[col] = (agg_df.set_index("Species")[col] / col_total).to_dict()
        else:
            rel_abundance[col] = {sp: 0.0 for sp in species}

    def _cell_color(col: str, sp: str, val: float) -> str:
        if val <= 0:
            return "white"
        if col in edna_cols:
            ra = rel_abundance[col].get(sp, 0.0)
            return edna_abundant_color if ra >= abundance_threshold else edna_rare_color
        return hist_color

    # ── COMMON NAMES (with uniqueness filter) ─────────────────────────────────
    common_names = {}
    if taxonomy_lookup is not None:
        try:
            raw_common = taxonomy_lookup.get_common_names_batch(
                list(species), hierarchical=True
            )
            # ◄ CHANGED: suppress common names shared by more than one species
            name_counts = Counter(v for v in raw_common.values() if v)
            common_names = {
                sp: name for sp, name in raw_common.items()
                if name and name_counts[name] == 1
            }
        except Exception as exc:
            print(f"Warning: batch common-name lookup failed: {exc}")

    def _make_ylabel(sp, common):
        italic_part = r"$\it{" + sp.replace(" ", r"\ ").replace("_", r"\ ") + r"}$"
        return italic_part + f"  ({common})" if common else italic_part

    y_labels = [_make_ylabel(sp, common_names.get(sp)) for sp in species]

    LEFT_PAD   = 0.15 + fs_ylabel * 0.10
    RIGHT_PAD  = 2.80
    TOP_PAD    = 0.15 + fs_xlabel * 0.14
    BOT_PAD    = 0.40

    axes_w = n_cols    * cell_size
    axes_h = n_species * cell_size
    fig_w  = LEFT_PAD + axes_w + RIGHT_PAD
    fig_h  = TOP_PAD  + axes_h + BOT_PAD
    plt.rcParams['svg.fonttype'] = 'path'
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    ax.set_position([
        LEFT_PAD / fig_w,
        BOT_PAD  / fig_h,
        axes_w   / fig_w,
        axes_h   / fig_h,
    ])

    for i, sp in enumerate(species):
        row = agg_df.loc[agg_df["Species"] == sp, all_cols]
        for j, col in enumerate(all_cols):
            val = float(row[col].values[0]) if len(row) else 0.0
            fc  = _cell_color(col, sp, val)
            ax.add_patch(
                Rectangle((j, i), 1, 1,
                           facecolor=fc, edgecolor="lightgray", linewidth=0.5)
            )

    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_species, 0)
    ax.set_aspect("equal", adjustable="box")

    ax.set_yticks(np.arange(n_species) + 0.5)
    ax.set_yticklabels(y_labels, fontsize=fs_ylabel)

    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(all_cols, rotation=45, ha="left",
                       fontsize=fs_xlabel, fontweight="bold")
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")

    for spine in ax.spines.values():
        spine.set_visible(False)

    png_path = f"{output_dir}/{name_prefix}_heatmap_vertical.png"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(png_path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png_path}")


def plot_COI_heatmap_horizontal(
    agg_df: pd.DataFrame,
    edna_cols: list[str],
    hist_cols: list[str],
    output_dir: str,
    taxonomy_lookup: Optional[NCBITaxonomyLookup] = None,
    name_prefix: str = "Musq",
    edna_abundant_color: str = "#2E86AB",
    edna_rare_color: str = "#C8C8C8",
    hist_color: str = "#FF6B35",
    cell_size: float = 0.30,
    abundance_threshold: float = 0.01,
    max_cols_per_panel: Optional[int] = None,
    class_order: Optional[list[str]] = None,
) -> None:

    import os, math, itertools
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    from collections import Counter  # ◄ CHANGED

    os.makedirs(output_dir, exist_ok=True)

    agg_df = sort_species_dataframe(agg_df, edna_cols + hist_cols, edna_cols)

    all_cols = edna_cols + hist_cols
    species = agg_df["Species"].values
    n_species = len(species)
    n_samples = len(all_cols)

    species_to_class = {}
    show_class_bar = False

    if taxonomy_lookup is not None:
        try:
            raw = taxonomy_lookup.get_ranks_batch(list(species), rank="class")
            species_to_class = {sp: (cls or "Unknown") for sp, cls in raw.items()}
            show_class_bar = True
        except Exception:
            pass

    if show_class_bar:
        _DEFAULT = [
            "Actinopterygii","Echinoidea","Asteroidea","Enteropneusta",
            "Malacostraca","Thecostraca","Maxillopoda","Gastropoda",
            "Bivalvia","Polychaeta","Clitellata","Gymnolaemata",
            "Pilidiophora","Enoplea","Tentaculata","Hydrozoa","Demospongiae",
        ]
        order = class_order if class_order else _DEFAULT

        def rank(sp):
            cls = species_to_class.get(sp, "Unknown")
            return order.index(cls) if cls in order else len(order)

        sorted_idx = sorted(range(n_species), key=lambda i: rank(species[i]))
        agg_df = agg_df.iloc[sorted_idx].reset_index(drop=True)
        species = agg_df["Species"].values

    palette = [
        "#D6EAF8","#D5F5E3","#FDEBD0","#E8DAEF","#FADBD8","#D1F2EB",
        "#FEF9E7","#EAF2FF","#F9EBEA","#E9F7EF","#FDF2F8","#EBF5FB"
    ]
    seen_classes = list(dict.fromkeys(species_to_class.get(sp, "Unknown") for sp in species))
    class_color = {cls: palette[i % len(palette)] for i, cls in enumerate(seen_classes)}

    rel_abundance = {}
    for col in edna_cols:
        total = agg_df[col].sum()
        rel_abundance[col] = (
            (agg_df.set_index("Species")[col] / total).to_dict()
            if total > 0 else {sp: 0.0 for sp in species}
        )

    def cell_color(col, sp, val):
        if val <= 0:
            return "white"
        if col in edna_cols:
            return edna_abundant_color if rel_abundance[col][sp] >= abundance_threshold else edna_rare_color
        return hist_color

    # ── COMMON NAMES (with uniqueness filter) ─────────────────────────────────
    common_names = {}
    if taxonomy_lookup is not None:
        try:
            raw_common = taxonomy_lookup.get_common_names_batch(list(species), hierarchical=True)
            # ◄ CHANGED: suppress common names shared by more than one species
            name_counts = Counter(v for v in raw_common.values() if v)
            common_names = {
                sp: name for sp, name in raw_common.items()
                if name and name_counts[name] == 1
            }
        except Exception:
            pass

    def fmt_label(sp):
        italic = r"$\it{" + sp.replace(" ", r"\ ") + r"}$"
        return italic + f"\n({common_names.get(sp,'')})" if common_names.get(sp) else italic

    if max_cols_per_panel and n_species > max_cols_per_panel:
        split = math.ceil(n_species / 2)
        species_chunks = [species[:split], species[split:]]
    else:
        species_chunks = [species]

    fig = plt.figure(figsize=(max(10, n_species * 0.25), 6))
    gs = fig.add_gridspec(
        nrows=2,
        ncols=len(species_chunks),
        height_ratios=[0.25, 1],
        hspace=0.05,
        wspace=0.25
    )

    for p, chunk in enumerate(species_chunks):

        ax_bar = fig.add_subplot(gs[0, p])
        ax = fig.add_subplot(gs[1, p])

        for j, sp in enumerate(chunk):
            row = agg_df.loc[agg_df["Species"] == sp, all_cols]
            for i, col in enumerate(all_cols):
                val = float(row[col].values[0]) if len(row) else 0.0
                ax.add_patch(Rectangle(
                    (j, i), 1, 1,
                    facecolor=cell_color(col, sp, val),
                    edgecolor="lightgray", linewidth=0.4
                ))

        ax.set_xlim(0, len(chunk))
        ax.set_ylim(n_samples, 0)
        ax.set_aspect("auto")

        ax.set_yticks(np.arange(n_samples) + 0.5)
        ax.set_yticklabels(all_cols, fontsize=10, fontweight="bold")

        ax.set_xticks(np.arange(len(chunk)) + 0.5)
        ax.set_xticklabels([fmt_label(sp) for sp in chunk],
                           rotation=90, ha="center", fontsize=8)

        for spine in ax.spines.values():
            spine.set_visible(False)

        if show_class_bar:
            ax_bar.set_xlim(0, len(chunk))
            ax_bar.set_ylim(0, 1)
            ax_bar.axis("off")

            for cls, group in itertools.groupby(chunk, key=lambda sp: species_to_class.get(sp, "Unknown")):
                group = list(group)
                start = chunk.tolist().index(group[0])
                width = len(group)

                ax_bar.add_patch(Rectangle(
                    (start, 0), width, 1,
                    facecolor=class_color.get(cls, "#EEEEEE"),
                    edgecolor="gray", linewidth=0.5
                ))

                if width > 1:
                    # ◄ CHANGED: abbreviate when span < 4 columns
                    label = cls[:3] + "." if width < 4 else cls
                    ax_bar.text(
                        start + width / 2,
                        0.5,
                        label,
                        ha="center", va="center",
                        fontsize=8, fontweight="bold"
                    )

    legend_handles = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor=edna_abundant_color,
               label=f"eDNA ≥{abundance_threshold:.1%}"),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=edna_rare_color,
               label=f"eDNA <{abundance_threshold:.1%}"),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=hist_color,
               label="Historical"),
    ]

    fig.legend(handles=legend_handles, loc="lower center", ncol=3, frameon=True)

    png_path = f"{output_dir}/{name_prefix}_heatmap_horizontal.png"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(png_path.replace(".png", ".pdf"), bbox_inches="tight")

    plt.close(fig)
    print(f"Saved: {png_path}")

def plot_heatmap_horizontal(
    agg_df: pd.DataFrame,
    edna_cols: list[str],
    hist_cols: list[str],
    output_dir: str,
    taxonomy_lookup: Optional[NCBITaxonomyLookup] = None,
    name_prefix: str = "Musq",
    edna_abundant_color: str = "#2E86AB",
    edna_rare_color: str = "#C8C8C8",
    hist_color: str = "#FF6B35",
    cell_size: float = 0.30,
    abundance_threshold: float = 0.01,
) -> None:
    """
    Square-cell heatmap — horizontal layout:
      - Species on X-axis (bottom, italic), rotated 45 °
      - Survey/study columns on Y-axis (left, bold)

    eDNA columns use two colors:
      - edna_abundant_color (default blue)  : relative abundance >= abundance_threshold
      - edna_rare_color     (default lt-gray): relative abundance >0 but < abundance_threshold
    Relative abundance is computed per eDNA column (species reads / column total).

    Historical columns use hist_color (default orange).
    Common names (if lookup provided) are appended in roman (non-italic).

    cell_size controls the physical size of every cell in inches.  Font sizes
    scale with cell_size so the ratio is preserved regardless of dataset size.

    Saved as PNG (600 dpi) and PDF.
    """
    os.makedirs(output_dir, exist_ok=True)
    agg_df = sort_species_dataframe(agg_df, edna_cols + hist_cols, edna_cols)
    all_cols  = edna_cols + hist_cols
    species   = agg_df["Species"].values
    n_species = len(species)
    n_cols    = len(all_cols)

    # ------------------------------------------------------------------
    # Font sizes that scale linearly with cell_size.
    # ------------------------------------------------------------------
    scale        = cell_size / 0.30
    fs_xlabel    = max(5, 9  * scale)        # species tick labels (x-axis)
    fs_ylabel    = max(5, 11 * scale)        # survey/column tick labels (y-axis)
    fs_axlabel   = max(6, 13 * scale)        # axis title ("Species")
    fs_legend    = max(6, 12 * scale)        # legend text
    marker_size  = max(6, 12 * scale)        # legend square marker

    # ------------------------------------------------------------------
    # Pre-compute relative abundance for every eDNA column
    # ------------------------------------------------------------------
    rel_abundance: dict[str, dict] = {}
    for col in edna_cols:
        col_total = agg_df[col].sum()
        if col_total > 0:
            rel_abundance[col] = (agg_df.set_index("Species")[col] / col_total).to_dict()
        else:
            rel_abundance[col] = {sp: 0.0 for sp in species}

    def _cell_color(col: str, sp: str, val: float) -> str:
        if val <= 0:
            return "white"
        if col in edna_cols:
            ra = rel_abundance[col].get(sp, 0.0)
            return edna_abundant_color if ra >= abundance_threshold else edna_rare_color
        return hist_color

    # ------------------------------------------------------------------
    # Optional common-name batch lookup
    # ------------------------------------------------------------------
    common_names = {}
    if taxonomy_lookup is not None:
        try:
            common_names = taxonomy_lookup.get_common_names_batch(
                list(species), hierarchical=True
            )
        except Exception as exc:
            print(f"Warning: batch common-name lookup failed: {exc}")

    def _make_xlabel(sp, common):
        italic_part = r"$\it{" + sp.replace(" ", r"\ ").replace("_", r"\ ") + r"}$"
        return italic_part + f"  ({common})" if common else italic_part

    x_labels = [_make_xlabel(sp, common_names.get(sp)) for sp in species]

    # ------------------------------------------------------------------
    # Figure sizing: axes dimensions are EXACT multiples of cell_size.
    # X-tick labels are rotated 45° — estimate their vertical footprint.
    # ------------------------------------------------------------------
    LEFT_PAD   = 0.15 + fs_ylabel * 0.10    # y-tick labels (survey names)
    RIGHT_PAD  = 0.30                        # small right buffer
    TOP_PAD    = 0.50                        # small top buffer
    # Rotated 45° labels: approximate height ≈ longest label chars * font_pt / 144
    max_label_chars = max((len(lb) for lb in x_labels), default=10)
    BOT_PAD = 0.15 + max_label_chars * fs_xlabel / 144 #+ 1.80  # +legend rows

    axes_w = n_species * cell_size
    axes_h = n_cols    * cell_size
    fig_w  = LEFT_PAD + axes_w + RIGHT_PAD
    fig_h  = TOP_PAD  + axes_h + BOT_PAD

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    ax.set_position([
        LEFT_PAD / fig_w,
        BOT_PAD  / fig_h,
        axes_w   / fig_w,
        axes_h   / fig_h,
    ])

    # ------------------------------------------------------------------
    # Draw cells
    # X axis → species (index j), Y axis → surveys (index i, top-down)
    # ------------------------------------------------------------------
    for i, col in enumerate(all_cols):
        for j, sp in enumerate(species):
            row = agg_df.loc[agg_df["Species"] == sp, [col]]
            val = float(row[col].values[0]) if len(row) else 0.0
            fc  = _cell_color(col, sp, val)
            ax.add_patch(
                Rectangle((j, i), 1, 1,
                           facecolor=fc, edgecolor="lightgray", linewidth=0.5)
            )

    ax.set_xlim(0, n_species)
    ax.set_ylim(n_cols, 0)      # inverted: first survey row at top
    ax.set_aspect("equal", adjustable="box")

    # X-axis — species labels at bottom, rotated 45 °
    ax.set_xticks(np.arange(n_species) + 0.5)
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=fs_xlabel)
    ax.xaxis.tick_bottom()
    # ax.set_xlabel("Species", fontsize=fs_axlabel, fontweight="bold")

    # Y-axis — survey / study column labels at left
    ax.set_yticks(np.arange(n_cols) + 0.5)
    ax.set_yticklabels(all_cols, fontsize=fs_ylabel, fontweight="bold")
    ax.yaxis.tick_left()

    for spine in ax.spines.values():
        spine.set_visible(False)

    # Legend placed below the rotated species labels
    ax.legend(
        handles=[
            Line2D([0], [0], marker="s", color="w",
                markerfacecolor=edna_abundant_color, markersize=marker_size,
                label=f"eDNA Detected (≥{abundance_threshold:.2%} rel. abund.)"),
            Line2D([0], [0], marker="s", color="w",
                markerfacecolor=edna_rare_color, markersize=marker_size,
                label=f"eDNA Detected (<{abundance_threshold:.2%} rel. abund.)"),
            Line2D([0], [0], marker="s", color="w",
                markerfacecolor=hist_color, markersize=marker_size,
                label="Historical Record"),
        ],
        # loc="upper center",          # anchor point inside legend
        # bbox_to_anchor=(0.5, -1.3), # centered, below x-axis labels
        # bbox_transform=ax.transAxes,
        loc='lower center',
        bbox_to_anchor=(0.5, 0.02), # 2% from the very bottom of the image
        bbox_transform=fig.transFigure,
        frameon=True,
        fontsize=fs_legend,
        ncol=3,                      # horizontal layout (looks cleaner below plot)
    )
    plt.subplots_adjust(bottom=0.30)

    png_path = f"{output_dir}/{name_prefix}_heatmap_horizontal.png"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(png_path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png_path}")
    
# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def merge_and_visualize_species(
    edna_files: dict,
    hist_files: Optional[dict] = None,
    output_dir: str = OUTPUT_DIR,
    name_prefix: str = "Musq",
    standard_zones: list = STANDARD_ZONES,
    upset_sort_by: str = "-degree",
    upset_sort_categories_by: str = "input",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Full pipeline: calls merge_files() then visualize_species().
    Returns (master_df, abundance_df).
    """
    master, agg = merge_files(
        edna_files     = edna_files,
        hist_files     = hist_files,
        output_dir     = output_dir,
        name_prefix    = name_prefix,
        standard_zones = standard_zones,
    )

    # Pass the saved CSV path — consistent with standalone visualize_species() usage
    abundance_csv = f"{output_dir}/{name_prefix}_abundance_table.csv"

    visualize_species(
        abundance_csv            = abundance_csv,
        edna_files               = edna_files,
        hist_files               = hist_files,
        output_dir               = output_dir,
        name_prefix              = name_prefix,
        upset_sort_by            = upset_sort_by,
        upset_sort_categories_by = upset_sort_categories_by,
    )

    return master, agg

def merge_files(
    edna_files: dict,
    hist_files: Optional[dict] = None,
    output_dir: str = OUTPUT_DIR,
    name_prefix: str = "Musq",
    standard_zones: list = STANDARD_ZONES,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Steps 1–4: Load, merge, and aggregate all input files.
    Saves two CSVs and returns (master_df, abundance_df).

    Returns
    -------
    master_df     : wide table with every zone/year column + historical columns
    abundance_df  : per-year summed counts + historical presence, sorted for plotting
    """
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("STEP 1–3: Building master wide table …")
    master = build_master_table(edna_files, hist_files, standard_zones)
    master.to_csv(f"{output_dir}/{name_prefix}_master_wide.csv", index=False)
    print(f"  {len(master)} unique species | {len(master.columns)} columns")

    print("\nSTEP 4: Aggregating to per-year abundance …")
    agg = aggregate_to_year_abundance(master, edna_files, hist_files)
    agg.to_csv(f"{output_dir}/{name_prefix}_abundance_table.csv", index=False)

    # edna_year_strs = [str(y) for y in sorted(set(edna_files.values()), key=_year_sort_key)]
    edna_year_strs = list(edna_files.values())if edna_files else []
    hist_labels    = list(hist_files.values()) if hist_files else []

    print("\n--- Detection summary ---")
    for col in edna_year_strs:
        print(f"  {col}: {(agg[col] > 0).sum()} species detected")
    for lbl in hist_labels:
        if lbl in agg.columns:
            print(f"  {lbl}: {int(agg[lbl].sum())} species recorded")

    return master, agg

def visualize_species(
    abundance_csv: str,
    edna_files: dict,
    hist_files: Optional[dict] = None,
    output_dir: str = OUTPUT_DIR,
    name_prefix: str = "Musq",
    upset_sort_by: str = "-degree",
    upset_sort_categories_by: str = "input",
    abundance_threshold = 0.001
) -> None:
    """
    Steps 5–6: Generate UpSet plot and vertical heatmap from a saved abundance CSV.

    Parameters
    ----------
    abundance_csv : path to the abundance table CSV saved by merge_files()
    edna_files    : used only to derive the year-label column names
    hist_files    : used only to derive the study-label column names
    """
    # Load the abundance table from disk
    if not os.path.exists(abundance_csv):
        raise FileNotFoundError(f"Abundance CSV not found: {abundance_csv}")

    abundance_df = pd.read_csv(abundance_csv)
    print(f"Loaded {len(abundance_df)} species from: {abundance_csv}")

    # edna_year_strs = [str(y) for y in sorted(set(edna_files.values()), key=_year_sort_key, reverse = True)]
    edna_year_strs = list(edna_files.values())if edna_files else []
    hist_labels    = list(hist_files.values()) if hist_files else []

    print("=" * 60)
    print("STEP 5: Generating UpSet plot …")
    plot_upset(
        abundance_df,
        edna_cols          = edna_year_strs,
        hist_cols          = hist_labels,
        output_dir         = output_dir,
        name_prefix        = name_prefix,
        sort_by            = upset_sort_by,
        sort_categories_by = upset_sort_categories_by,
    )

    print("\nSTEP 6: Generating V & H heatmaps  …")

    lookup = NCBITaxonomyLookup(names_dmp_path="/mnt/e/projects/databases/tax_dump/ncbi_taxdump/names.dmp")
    
    plot_heatmap_vertical(
        abundance_df,
        edna_cols   = edna_year_strs,
        hist_cols   = hist_labels,
        output_dir  = output_dir,
        name_prefix = name_prefix,
        taxonomy_lookup = lookup,
        abundance_threshold = abundance_threshold
        # max_rows_per_panel = 50,
        # class_order=["Actinopteri", "Echinoidea", "Asteroidea", "Enteropneusta", "Malacostraca", "Hexanauplia", "Thecostraca", "Gastropoda", "Bivalvia", "Polychaeta", "Clitellata", "Gymnolaemata", "Pilidiophora", "Tentaculata", "Hydrozoa", "Demospongiae"]
    )

    plot_COI_heatmap_vertical(
        abundance_df,
        edna_cols   = edna_year_strs,
        hist_cols   = hist_labels,
        output_dir  = output_dir,
        name_prefix = name_prefix,
        taxonomy_lookup = lookup,
        abundance_threshold = abundance_threshold,
        max_rows_per_panel = 50,
        class_order=["Actinopteri", "Echinoidea", "Asteroidea", "Enteropneusta", "Malacostraca", "Hexanauplia", "Thecostraca", "Gastropoda", "Bivalvia", "Polychaeta", "Clitellata", "Gymnolaemata", "Pilidiophora", "Tentaculata", "Hydrozoa", "Demospongiae"]
    )


    # plot_COI_heatmap_horizontal(
    #     abundance_df,
    #     edna_cols   = edna_year_strs,
    #     hist_cols   = hist_labels,
    #     output_dir  = output_dir,
    #     name_prefix = name_prefix,
    #     taxonomy_lookup = lookup,
    #     abundance_threshold = abundance_threshold
    #     # max_rows_per_panel = 50,
    #     # class_order=["Actinopteri", "Echinoidea", "Asteroidea", "Enteropneusta", "Malacostraca", "Hexanauplia", "Thecostraca", "Gastropoda", "Bivalvia", "Polychaeta", "Clitellata", "Gymnolaemata", "Pilidiophora", "Tentaculata", "Hydrozoa", "Demospongiae"]
    # )

    print("\n✓ Visualization complete.")


# ─────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    master_df, abundance_df = merge_and_visualize_species(
        edna_files    = EDNA_FILES,
        hist_files    = HIST_FILES,
        output_dir    = OUTPUT_DIR,
        name_prefix   = f"Musq-{MARKER}",
        standard_zones= STANDARD_ZONES,
    )
    print(abundance_df.head(10).to_string(index=False))
