"""
End-to-end pipeline for sister chromatid composition inference.

This module provides the main entry point for running sister chromatid
state inference from BAM/BCF files with intermediate caching.

Pipeline Steps:
    1. BAM/BCF → Read Counts (extract ref/alt/other counts at SNP positions)
    2. Read Counts → Binning (optional: merge adjacent SNPs into units)
    3. Binned Data → Rough Classification (initial state labels)
    4. Rough Classification → Parameter Estimation (emission & transition)
    5. Parameters + Binned Data → HMM Inference (Viterbi decoding)
    6. HMM States → Blocks (contiguous segment grouping)
    7. Unit-level coordinates → Original SNP coordinates (conversion)
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from hic_basic.genome import Region, parse_ucsc_region
from hic_basic.plot.track import _sort_chromosome_position_table, plot_genome_phased_read_count_track

from .cache import CacheManager
from .detector_mii import (
    SisterChromatidBlock,
    SisterChromatidHMM,
    SisterChromatidResults,
    estimate_emissions_from_counts,
    estimate_transitions_rough,
    get_default_transition_params,
    rough_classification,
)

logger = logging.getLogger(__name__)


STATE_ORDER = ["MM", "PM", "PP", "U"]
STATE_COLORS = {
    "PP": "#d62728",
    "PM": "#2ca02c",
    "MM": "#1f77b4",
    "U": "#7f7f7f",
    "AMBIG": "#ffbf00",
}
STATE_TO_HMM_INDEX = {"PP": 0, "PM": 1, "MM": 2, "U": 3}
GAP_LABEL = "__GAP__"
AMBIG_LABEL = "AMBIG"
ROUGH_LABEL_ORDER = ["MM", "PM", "PP", "U", AMBIG_LABEL]


def bin_by_snp(read_counts_df: pd.DataFrame, n_snps: int = 10) -> pd.DataFrame:
    """
    Merge adjacent SNPs into units based on SNP count.
    
    Args:
        read_counts_df: DataFrame with columns [chromosome or chrom, pos or position, ref, alt, other]
        n_snps: Number of SNPs per unit (default: 10)
    
    Returns:
        DataFrame with aggregated units, columns: [chromosome, start_pos, end_pos, ref, alt, other, n_snps, position]
    """
    if len(read_counts_df) == 0:
        return read_counts_df
    
    # Normalize column names
    df = read_counts_df.copy()
    if 'chrom' in df.columns and 'chromosome' not in df.columns:
        df['chromosome'] = df['chrom']
    if 'pos' in df.columns and 'position' not in df.columns:
        df['position'] = df['pos']
    
    # Create unit IDs
    unit_ids = np.arange(len(df)) // n_snps
    df['unit_id'] = unit_ids
    
    # Aggregate
    aggregated = df.groupby('unit_id').agg({
        'chromosome': 'first',
        'position': ['first', 'last'],
        'ref': 'sum',
        'alt': 'sum',
        'other': 'sum',
        'unit_id': 'size',
    }).reset_index()
    
    # Flatten column names
    aggregated.columns = ['unit_id', 'chromosome', 'start_pos', 'end_pos', 'ref', 'alt', 'other', 'n_snps']
    
    # Use start_pos as position for HMM
    aggregated['position'] = aggregated['start_pos']
    
    return aggregated


def _write_plotly_png(fig, output_path: Union[str, Path], width: int = 1200, height: int = 700) -> None:
    """Write a Plotly figure to PNG, with a simple matplotlib fallback."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.write_image(str(output_path), width=width, height=height)
    except Exception as exc:  # pragma: no cover - only used when kaleido/chrome is unavailable
        import matplotlib.pyplot as plt

        fallback_title = getattr(fig.layout.title, "text", None) or output_path.name
        fallback_text = f"{fallback_title}\n\nPlotly image export failed:\n{exc}"
        fallback_fig, ax = plt.subplots(figsize=(width / 150, height / 150))
        ax.axis("off")
        ax.text(0.5, 0.5, fallback_text, ha="center", va="center", wrap=True)
        fallback_fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fallback_fig)


def _diagnostic_path(prefix: Union[str, Path], suffix: str) -> Path:
    return Path(f"{prefix}.{suffix}.png")


def _blank_diagnostic_png(title: str, output_path: Union[str, Path]) -> None:
    import plotly.graph_objects as go

    fig = go.Figure()
    fig.add_annotation(
        text=title,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(size=30),
    )
    fig.update_layout(
        template="plotly_white",
        title=title,
        font=dict(size=22),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )
    _write_plotly_png(fig, output_path)


def _prepare_counts_for_estimation(
    read_counts_df: pd.DataFrame,
    binning_method: Optional[str],
    binning_snp_number: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    counts_df = read_counts_df.copy()
    if "chrom" in counts_df.columns and "chromosome" not in counts_df.columns:
        counts_df["chromosome"] = counts_df["chrom"]
    if "pos" in counts_df.columns and "position" not in counts_df.columns:
        counts_df["position"] = counts_df["pos"]
    required = {"chromosome", "position", "ref", "alt", "other"}
    missing = required - set(counts_df.columns)
    if missing:
        raise ValueError(f"read_counts_df is missing required columns: {sorted(missing)}")
    counts_df = counts_df.sort_values(["chromosome", "position"]).reset_index(drop=True)
    original_counts_df = counts_df.copy()
    if binning_method == "by_snp":
        counts_df = bin_by_snp(counts_df, n_snps=binning_snp_number).reset_index(drop=True)
    elif binning_method is not None:
        raise ValueError("binning_method must be None or 'by_snp'")
    return counts_df, original_counts_df


def _effective_depth_threshold(
    counts_df: pd.DataFrame,
    depth_threshold: int,
    binning_method: Optional[str],
) -> pd.Series:
    """
    Convert per-SNP depth threshold to the current observation unit.

    In single-SNP mode this is unchanged. In SNP-binned mode, read counts are
    summed across SNPs, so the minimum total depth must scale with the number of
    SNPs that contributed to each bin.
    """
    if binning_method == "by_snp":
        if "n_snps" not in counts_df.columns:
            raise ValueError("Binned counts must include an n_snps column")
        return counts_df["n_snps"].astype(float) * float(depth_threshold)
    return pd.Series(float(depth_threshold), index=counts_df.index)


def _rough_classification_with_thresholds(
    read_counts: pd.DataFrame,
    effective_depth_threshold: pd.Series,
    other_threshold: float = 0.4,
    balance_threshold: float = 0.25,
) -> pd.Series:
    labels = []
    effective_depth_threshold = effective_depth_threshold.reset_index(drop=True)
    read_counts = read_counts.reset_index(drop=True)

    for i, row in read_counts.iterrows():
        ref = row["ref"]
        alt = row["alt"]
        other = row["other"]
        total = ref + alt + other

        if total < effective_depth_threshold.iloc[i]:
            labels.append("U")
            continue

        other_frac = other / total if total > 0 else 1.0
        if other_frac > other_threshold:
            labels.append("U")
            continue

        ref_frac = ref / total if total > 0 else 0.0
        alt_frac = alt / total if total > 0 else 0.0

        if ref_frac > 1 - balance_threshold:
            labels.append("PP")
        elif alt_frac > 1 - balance_threshold:
            labels.append("MM")
        elif abs(ref_frac - alt_frac) < balance_threshold and ref_frac > 0.2 and alt_frac > 0.2:
            labels.append("PM")
        else:
            labels.append(AMBIG_LABEL)

    return pd.Series(labels)


def _u_emission_prior_from_other(other_prob: float) -> Tuple[float, float, float]:
    """Build a symmetric ref/alt U-emission prior with modestly elevated other."""
    other_prob = float(other_prob)
    if not 0.0 < other_prob < 1.0:
        raise ValueError("u_emission_other_prior must be between 0 and 1")
    ref_alt_prob = (1.0 - other_prob) / 2.0
    return (ref_alt_prob, ref_alt_prob, other_prob)


def _rough_label_segments(counts_df: pd.DataFrame, rough_labels: pd.Series) -> pd.DataFrame:
    """Summarize chromosome-aware run-length segments in raw rough labels."""
    df = counts_df.reset_index(drop=True).copy()
    labels = rough_labels.reset_index(drop=True)
    df["rough_label"] = labels
    rows = []

    for chrom, chrom_df in df.groupby("chromosome", sort=False):
        run_start = None
        run_label = None
        prev_idx = None
        for idx, row in chrom_df.iterrows():
            label = row["rough_label"]
            if run_start is None:
                run_start = idx
                run_label = label
            elif label != run_label:
                run_df = df.iloc[run_start:prev_idx + 1]
                rows.append({
                    "chromosome": chrom,
                    "label": run_label,
                    "start_index": int(run_start),
                    "end_index": int(prev_idx),
                    "start_pos": int(run_df["position"].iloc[0]),
                    "end_pos": int(run_df["position"].iloc[-1]),
                    "length_bins": int(prev_idx - run_start + 1),
                })
                run_start = idx
                run_label = label
            prev_idx = idx

        if run_start is not None and prev_idx is not None:
            run_df = df.iloc[run_start:prev_idx + 1]
            rows.append({
                "chromosome": chrom,
                "label": run_label,
                "start_index": int(run_start),
                "end_index": int(prev_idx),
                "start_pos": int(run_df["position"].iloc[0]),
                "end_pos": int(run_df["position"].iloc[-1]),
                "length_bins": int(prev_idx - run_start + 1),
            })

    return pd.DataFrame(rows)


def _filter_short_rough_segments(
    rough_labels: pd.Series,
    segment_table: pd.DataFrame,
    min_segment_bins: int,
) -> pd.Series:
    """
    Mark short rough-label segments as gaps for transition estimation.

    Gaps are not joined across: A-gap-B contributes neither A->gap, gap->B,
    nor A->B transition counts.
    """
    if min_segment_bins < 1:
        raise ValueError("transition_min_segment_bins must be >= 1")
    filtered = rough_labels.reset_index(drop=True).astype(object).copy()
    if segment_table.empty:
        return filtered

    short_segments = segment_table[segment_table["length_bins"] < min_segment_bins]
    for _, segment in short_segments.iterrows():
        filtered.iloc[int(segment["start_index"]):int(segment["end_index"]) + 1] = GAP_LABEL
    return filtered


def _labels_with_chromosome_gap_separators(
    counts_df: pd.DataFrame,
    labels: pd.Series,
) -> pd.Series:
    """Insert gap separators between chromosomes before transition counting."""
    df = counts_df.reset_index(drop=True).copy()
    df["_label"] = labels.reset_index(drop=True).astype(object)
    values = []
    for _, chrom_df in df.groupby("chromosome", sort=False):
        values.extend(chrom_df["_label"].tolist())
        values.append(GAP_LABEL)
    return pd.Series(values, dtype=object)


def _apply_by_snp_u_transition_prior(
    transition_param_dict: Dict[str, float],
    transition_scale_info: Dict[str, Any],
    transition_u_gamma_prior: Optional[float],
    transition_u_gamma_factor: float,
    transition_u_rho_prior: float,
) -> None:
    """Use prior U transition parameters for SNP-binned transition estimation."""
    if transition_u_gamma_factor < 0:
        raise ValueError("transition_u_gamma_factor must be non-negative")
    if not 0 <= transition_u_rho_prior <= 0.32:
        raise ValueError("transition_u_rho_prior must be between 0 and 0.32")

    alpha = transition_param_dict["alpha"]
    beta = transition_param_dict["beta"]
    fallback_used = False
    if transition_u_gamma_prior is None:
        gamma_base = min(alpha, beta)
        if gamma_base <= 0:
            defaults = get_default_transition_params()
            gamma_base = min(defaults["alpha"], defaults["beta"])
            fallback_used = True
            print(
                "[mII estimate] by-SNP U gamma prior fallback: "
                "estimated alpha/beta has zero support "
                f"(alpha={alpha:.6f}, beta={beta:.6f}); "
                f"using default base={gamma_base:.6f}"
            )
        gamma_prior = transition_u_gamma_factor * gamma_base
    else:
        gamma_prior = float(transition_u_gamma_prior)
        if gamma_prior < 0:
            raise ValueError("transition_u_gamma_prior must be non-negative")

    max_gamma = max(0.0, 0.99 - max(alpha, beta))
    gamma_prior = min(gamma_prior, max_gamma)

    transition_param_dict["gamma"] = gamma_prior
    transition_param_dict["rho"] = transition_u_rho_prior
    transition_scale_info.update({
        "by_snp_u_prior_applied": True,
        "transition_u_gamma_prior": gamma_prior,
        "transition_u_gamma_factor": transition_u_gamma_factor,
        "transition_u_rho_prior": transition_u_rho_prior,
        "u_labels_treated_as_gaps_for_transition": True,
        "transition_u_gamma_fallback_used": fallback_used,
    })


def _median_observation_distance_mb(counts_df: pd.DataFrame) -> float:
    """Estimate one HMM step length in Mb from adjacent observations."""
    distances = []
    for _, chrom_df in counts_df.groupby("chromosome", sort=False):
        pos = chrom_df["position"].astype(float).sort_values()
        diffs = pos.diff().dropna()
        distances.extend(diffs[diffs > 0].tolist())
    if not distances:
        return 0.0
    return float(np.median(distances) / 1_000_000.0)


def _median_observation_distance_snps(
    counts_df: pd.DataFrame,
    binning_method: Optional[str],
) -> float:
    """Estimate one HMM step length in SNP units."""
    if binning_method == "by_snp":
        if "n_snps" not in counts_df.columns:
            raise ValueError("Binned counts must include an n_snps column")
        return float(counts_df["n_snps"].median())
    return 1.0


def _scale_transition_params(
    params: Dict[str, float],
    transition_param_unit: str,
    median_distance_mb: float,
    median_distance_snps: float,
) -> Tuple[Dict[str, float], Dict[str, Any]]:
    """
    Scale caller-specified transition parameters to the HMM's per-step unit.

    Estimated transition parameters are already per observation step and should
    not be passed through this helper.
    """
    valid_units = {"step", "mb", "snp"}
    if transition_param_unit not in valid_units:
        raise ValueError(f"transition_param_unit must be one of {sorted(valid_units)}")

    params = {k: float(v) for k, v in params.items()}
    if transition_param_unit == "step":
        scale_info = {
            "transition_param_unit": transition_param_unit,
            "median_distance_mb": median_distance_mb,
            "median_distance_snps": median_distance_snps,
            "scaling_applied": False,
        }
        return params, scale_info

    if transition_param_unit == "mb":
        scaled = {
            key: 1.0 - np.exp(-value * median_distance_mb)
            for key, value in params.items()
        }
    else:
        scaled = {
            key: 1.0 - (1.0 - value) ** median_distance_snps
            for key, value in params.items()
        }

    scaled["rho"] = min(scaled["rho"], 0.32)
    scale_info = {
        "transition_param_unit": transition_param_unit,
        "median_distance_mb": median_distance_mb,
        "median_distance_snps": median_distance_snps,
        "scaling_applied": True,
        "input_transition_params": params,
    }
    return scaled, scale_info


def _transition_params_from_source(
    transition_params: Union[str, Dict[str, float]],
    alpha: float,
    beta: float,
    gamma: float,
    rho: float,
    rough_labels_for_transition: Optional[pd.Series],
    transition_param_unit: str,
    median_distance_mb: float,
    median_distance_snps: float,
) -> Tuple[Dict[str, float], str, Dict[str, Any]]:
    if isinstance(transition_params, dict):
        required = {"alpha", "beta", "gamma", "rho"}
        missing = required - set(transition_params)
        if missing:
            raise ValueError(f"transition_params dict is missing keys: {sorted(missing)}")
        raw_params = {k: float(transition_params[k]) for k in ["alpha", "beta", "gamma", "rho"]}
        scaled_params, scale_info = _scale_transition_params(
            raw_params,
            transition_param_unit=transition_param_unit,
            median_distance_mb=median_distance_mb,
            median_distance_snps=median_distance_snps,
        )
        return scaled_params, "specified", scale_info

    if transition_params == "estimated":
        if rough_labels_for_transition is None:
            raise ValueError("rough_labels are required when transition_params='estimated'")
        alpha_est, beta_est, gamma_est, rho_est = estimate_transitions_rough(rough_labels_for_transition)
        scale_info = {
            "transition_param_unit": "step",
            "median_distance_mb": median_distance_mb,
            "median_distance_snps": median_distance_snps,
            "scaling_applied": False,
            "estimated_parameters_are_per_step": True,
        }
        return {
            "alpha": alpha_est,
            "beta": beta_est,
            "gamma": gamma_est,
            "rho": rho_est,
        }, "estimated", scale_info

    if transition_params == "default":
        raw_params = {"alpha": alpha, "beta": beta, "gamma": gamma, "rho": rho}
        scaled_params, scale_info = _scale_transition_params(
            raw_params,
            transition_param_unit=transition_param_unit,
            median_distance_mb=median_distance_mb,
            median_distance_snps=median_distance_snps,
        )
        return scaled_params, "specified", scale_info

    raise ValueError("transition_params must be 'default', 'estimated', or a parameter dict")


def _emission_stats_table(
    counts_df: pd.DataFrame,
    rough_labels: Optional[pd.Series],
    emission_params: Dict[str, Tuple[float, float, float]],
    source: str,
) -> pd.DataFrame:
    rows = []
    for label in STATE_ORDER:
        p_ref, p_alt, p_other = emission_params[label]
        row = {
            "label": label,
            "ref_ep": p_ref,
            "alt_ep": p_alt,
            "other_ep": p_other,
        }
        if source.startswith("estimated") and rough_labels is not None:
            mask = rough_labels.reset_index(drop=True) == label
            state_counts = counts_df.loc[mask, ["ref", "alt", "other"]]
            ref_count = int(state_counts["ref"].sum())
            alt_count = int(state_counts["alt"].sum())
            other_count = int(state_counts["other"].sum())
            row.update({
                "ref_count": ref_count,
                "alt_count": alt_count,
                "other_count": other_count,
                "total_count": ref_count + alt_count + other_count,
                "n_units": int(mask.sum()),
            })
        rows.append(row)
    df = pd.DataFrame(rows)
    if source.startswith("estimated"):
        cols = ["label", "n_units", "ref_count", "alt_count", "other_count", "total_count",
                "ref_ep", "alt_ep", "other_ep"]
        return df[[c for c in cols if c in df.columns]]
    return df[["label", "ref_ep", "alt_ep", "other_ep"]]


def _plot_table_png(table_df: pd.DataFrame, title: str, output_path: Union[str, Path]) -> None:
    import plotly.graph_objects as go

    display_df = table_df.copy()
    for col in display_df.columns:
        if pd.api.types.is_float_dtype(display_df[col]):
            display_df[col] = display_df[col].map(lambda x: f"{x:.4f}")
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=list(display_df.columns),
            fill_color="#e8eef7",
            align="center",
            font=dict(size=18),
            height=36,
        ),
        cells=dict(
            values=[display_df[col].tolist() for col in display_df.columns],
            fill_color="white",
            align="center",
            font=dict(size=17),
            height=32,
        ),
    )])
    fig.update_layout(template="plotly_white", title=title, font=dict(size=18), title_font=dict(size=26))
    _write_plotly_png(fig, output_path, width=1500, height=500)


def _plot_transition_matrix_png(
    transition_matrix: np.ndarray,
    transition_param_dict: Dict[str, float],
    source: str,
    output_path: Union[str, Path],
) -> None:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    param_df = pd.DataFrame({
        "parameter": ["alpha", "beta", "gamma", "rho"],
        "value": [transition_param_dict[k] for k in ["alpha", "beta", "gamma", "rho"]],
        "definition": [
            "Pure state (MM/PP) -> mixed PM probability",
            "Mixed PM -> pure state total probability",
            "Biological state (MM/PM/PP) -> U probability",
            "U -> each biological state probability",
        ],
    })
    param_df["value"] = param_df["value"].map(lambda x: f"{x:.6f}")
    ordered_idx = [STATE_TO_HMM_INDEX[state] for state in STATE_ORDER]
    ordered_matrix = transition_matrix[np.ix_(ordered_idx, ordered_idx)]

    fig = make_subplots(
        rows=1,
        cols=2,
        column_widths=[0.54, 0.46],
        specs=[[{"type": "heatmap"}, {"type": "table"}]],
        subplot_titles=("Transition matrix", "Parameters"),
    )
    fig.add_trace(
        go.Heatmap(
            z=ordered_matrix,
            x=STATE_ORDER,
            y=STATE_ORDER,
            zmin=0,
            zmax=1,
            colorscale="Blues",
            text=np.vectorize(lambda x: f"{x:.4f}")(ordered_matrix),
            texttemplate="%{text}",
            textfont=dict(size=18),
            colorbar=dict(title=dict(text="prob", font=dict(size=16)), tickfont=dict(size=14)),
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Table(
            columnwidth=[0.18, 0.18, 0.64],
            header=dict(
                values=list(param_df.columns),
                fill_color="#e8eef7",
                align="center",
                font=dict(size=18),
                height=38,
            ),
            cells=dict(
                values=[param_df[c].tolist() for c in param_df.columns],
                align=["center", "center", "left"],
                font=dict(size=17),
                height=34,
            ),
        ),
        row=1,
        col=2,
    )
    title_source = "estimated from rough labels" if source == "estimated" else "specified/fixed"
    fig.update_xaxes(title=dict(text="to state", font=dict(size=20)), tickfont=dict(size=18), row=1, col=1)
    fig.update_yaxes(
        title=dict(text="from state", font=dict(size=20)),
        tickfont=dict(size=18),
        autorange="reversed",
        row=1,
        col=1,
    )
    fig.update_layout(
        template="plotly_white",
        title=f"Transition matrix ({title_source})",
        title_font=dict(size=28),
        font=dict(size=18),
        legend=dict(font=dict(size=16)),
    )
    _write_plotly_png(fig, output_path, width=1600, height=650)


def _plot_rough_label_count_png(
    counts_df: pd.DataFrame,
    rough_labels: pd.Series,
    output_path: Union[str, Path],
    title: str = "Rough label classification",
) -> None:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    df = counts_df.copy().reset_index(drop=True)
    labels = rough_labels.reset_index(drop=True)
    df["rough_label"] = labels
    df["total_count"] = df["ref"] + df["alt"] + df["other"]
    df["x_id"] = np.arange(len(df))
    label_order = ROUGH_LABEL_ORDER
    state_to_y = {label: i for i, label in enumerate(label_order)}
    df["label_y"] = df["rough_label"].map(state_to_y)
    count_columns = ["ref", "alt", "other"]
    count_values = df[count_columns].to_numpy(dtype=float).ravel()
    finite_count_values = count_values[np.isfinite(count_values)]
    if len(finite_count_values) == 0:
        y_cap = 1.0
    else:
        y_cap = float(np.quantile(finite_count_values, 0.995))
        y_cap = max(y_cap, 1.0)
    n_clipped = int(np.sum(finite_count_values > y_cap))

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        row_heights=[0.42, 0.58],
        vertical_spacing=0.08,
        subplot_titles=("genome_wide_label", "genome_wide_count"),
    )

    for state in label_order:
        mask = df["rough_label"] == state
        if mask.any():
            fig.add_trace(
                go.Scattergl(
                    x=df.loc[mask, "x_id"],
                    y=df.loc[mask, "label_y"],
                    mode="markers",
                    marker=dict(color=STATE_COLORS.get(state, "#333333"), size=7),
                    name=f"{state} ({int(mask.sum())})",
                ),
                row=1,
                col=1,
            )

    count_traces = [
        ("ref", "#d62728"),
        ("alt", "#1f77b4"),
        ("other", "#7f7f7f"),
    ]
    for col, color in count_traces:
        fig.add_trace(
            go.Scattergl(
                x=df["x_id"],
                y=df[col],
                mode="lines+markers",
                line=dict(color=color, width=1.5),
                marker=dict(size=5),
                name=col,
            ),
            row=2,
            col=1,
        )

    chrom_starts = df.groupby("chromosome", sort=False)["x_id"].min()
    chrom_ends = df.groupby("chromosome", sort=False)["x_id"].max()
    for chrom, start in chrom_starts.iloc[1:].items():
        fig.add_vline(x=float(start), line_width=1, line_color="black", opacity=0.5)
    for chrom in chrom_starts.index:
        mid = (chrom_starts[chrom] + chrom_ends[chrom]) / 2
        fig.add_annotation(
            x=mid,
            y=1.08,
            xref="x",
            yref="paper",
            text=str(chrom),
            showarrow=False,
            font=dict(size=16),
        )

    fig.update_yaxes(
        tickmode="array",
        tickvals=list(range(len(label_order))),
        ticktext=label_order,
        title=dict(text="label", font=dict(size=18)),
        tickfont=dict(size=16),
        row=1,
        col=1,
    )
    fig.update_yaxes(
        title=dict(text="read count", font=dict(size=18)),
        tickfont=dict(size=16),
        range=[0, y_cap * 1.05],
        row=2,
        col=1,
    )
    fig.update_xaxes(title=dict(text="genome-wide bin index", font=dict(size=18)), tickfont=dict(size=16), row=2, col=1)
    fig.add_annotation(
        text=f"count y-axis capped at p99.5={y_cap:.1f}; clipped points={n_clipped}",
        x=0.5,
        y=0.48,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(size=16),
        bgcolor="rgba(255,255,255,0.75)",
    )
    fig.update_layout(
        template="plotly_white",
        title=title,
        title_font=dict(size=28),
        font=dict(size=18),
        height=720,
        width=1400,
        legend=dict(orientation="h", yanchor="bottom", y=-0.28, font=dict(size=16)),
    )
    _write_plotly_png(fig, output_path, width=1400, height=720)


def _read_mii_plot_table(path: Union[str, Path]) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.suffix == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def _combine_mii_state_and_count_tables(
    states_df: pd.DataFrame,
    counts_df: Optional[pd.DataFrame],
) -> pd.DataFrame:
    df = states_df.copy()
    if "chrom" in df.columns and "chromosome" not in df.columns:
        df["chromosome"] = df["chrom"]
    if "pos" in df.columns and "position" not in df.columns:
        df["position"] = df["pos"]

    if counts_df is None:
        return df

    counts = counts_df.copy()
    if "chrom" in counts.columns and "chromosome" not in counts.columns:
        counts["chromosome"] = counts["chrom"]
    if "pos" in counts.columns and "position" not in counts.columns:
        counts["position"] = counts["pos"]

    merge_cols = [col for col in ["chromosome", "position"] if col in df.columns and col in counts.columns]
    count_cols = [col for col in ["ref", "alt", "other", "start_pos", "end_pos", "n_snps"] if col in counts.columns]
    existing_count_cols = [col for col in ["ref", "alt", "other"] if col in df.columns]
    if len(merge_cols) < 2 or existing_count_cols:
        if (
            not existing_count_cols
            and len(df) == len(counts)
            and "chromosome" in counts.columns
            and set(df.get("chromosome", pd.Series(dtype=str)).astype(str)) <= {"all"}
        ):
            for col in ["chromosome", "position"] + count_cols:
                if col in counts.columns:
                    df[col] = counts[col].to_numpy()
        return df
    return df.merge(counts[merge_cols + count_cols], on=merge_cols, how="left")


def _load_mii_final_assignment_output(pipeline_output: Any) -> pd.DataFrame:
    """Load a final-state/count table from a dataframe, file, directory, or pipeline result dict."""
    if isinstance(pipeline_output, pd.DataFrame):
        return pipeline_output.copy()

    if isinstance(pipeline_output, SisterChromatidResults):
        return pd.DataFrame({
            "chromosome": pipeline_output.chromosome,
            "position": pipeline_output.positions,
            "state": pipeline_output.states,
        })

    if isinstance(pipeline_output, dict):
        output_files = pipeline_output.get("output_files", pipeline_output)
        if "states" in output_files:
            states_df = _read_mii_plot_table(output_files["states"])
            counts_path = output_files.get("read_counts") or output_files.get("counts")
            counts_df = _read_mii_plot_table(counts_path) if counts_path is not None else None
            return _combine_mii_state_and_count_tables(states_df, counts_df)
        if "results" in pipeline_output and isinstance(pipeline_output["results"], SisterChromatidResults):
            return _load_mii_final_assignment_output(pipeline_output["results"])

    if isinstance(pipeline_output, (str, Path)):
        path = Path(pipeline_output)
        if path.is_dir():
            state_files = sorted(path.glob("*_states.tsv"))
            if not state_files:
                raise FileNotFoundError(f"No *_states.tsv files found in {path}")
            states_df = pd.concat([_read_mii_plot_table(p) for p in state_files], ignore_index=True)
            count_files = sorted(path.glob("*_read_counts.tsv"))
            counts_df = None
            if count_files:
                counts_df = pd.concat([_read_mii_plot_table(p) for p in count_files], ignore_index=True)
            return _combine_mii_state_and_count_tables(states_df, counts_df)
        return _read_mii_plot_table(path)

    raise TypeError(
        "pipeline_output must be a pandas DataFrame, pipeline result dict, "
        "output directory, or tabular output file"
    )


def _normalize_mii_final_assignment_table(df: pd.DataFrame, require_counts: bool = True) -> pd.DataFrame:
    df = df.copy().reset_index(drop=True)
    if "chrom" in df.columns and "chromosome" not in df.columns:
        df["chromosome"] = df["chrom"]
    if "pos" in df.columns and "position" not in df.columns:
        df["position"] = df["pos"]
    if "label" in df.columns and "state" not in df.columns:
        df["state"] = df["label"]
    if "assigned_label" in df.columns and "state" not in df.columns:
        df["state"] = df["assigned_label"]

    missing = {"chromosome", "state"} - set(df.columns)
    if missing:
        raise ValueError(f"pipeline_output is missing required columns: {sorted(missing)}")

    if "position" not in df.columns:
        if {"start_pos", "end_pos"} <= set(df.columns):
            df["position"] = (df["start_pos"].astype(float) + df["end_pos"].astype(float)) / 2
        else:
            raise ValueError("pipeline_output must contain position, pos, or start_pos/end_pos columns")

    count_missing = {"ref", "alt", "other"} - set(df.columns)
    if require_counts and count_missing:
        raise ValueError(
            "Final label/count visualization requires ref, alt, and other columns. "
            f"Missing: {sorted(count_missing)}"
        )

    if "start_pos" not in df.columns:
        df["start_pos"] = df["position"]
    if "end_pos" not in df.columns:
        df["end_pos"] = df["position"]

    numeric_cols = ["position", "start_pos", "end_pos"]
    numeric_cols.extend([col for col in ["ref", "alt", "other"] if col in df.columns])
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["position", "start_pos", "end_pos"]).reset_index(drop=True)
    return df


def _parse_mii_plot_region(region: Optional[Union[str, List[Tuple[str, int]], Region, Tuple[Any, ...]]]):
    if region is None:
        return None
    if isinstance(region, Region):
        return region.r
    if isinstance(region, str):
        chrom, start, end = parse_ucsc_region(region)
        return [(chrom, start), (chrom, end)]
    if isinstance(region, tuple) and len(region) == 3:
        chrom, start, end = region
        return [(chrom, start), (chrom, end)]
    if isinstance(region, tuple):
        region = list(region)
    return Region(region).r


def _filter_mii_final_assignment_region(
    df: pd.DataFrame,
    region: Optional[Union[str, List[Tuple[str, int]], Region, Tuple[Any, ...]]],
) -> pd.DataFrame:
    parsed = _parse_mii_plot_region(region)
    if parsed is None:
        return df

    (start_chrom, start_pos), (end_chrom, end_pos) = parsed
    start_pos = -np.inf if start_pos is None else start_pos
    end_pos = np.inf if end_pos is None else end_pos

    if start_chrom == end_chrom:
        mask = (
            (df["chromosome"].astype(str) == str(start_chrom))
            & (df["end_pos"] >= start_pos)
            & (df["start_pos"] <= end_pos)
        )
        return df.loc[mask].reset_index(drop=True)

    chrom_order = list(dict.fromkeys(df["chromosome"].astype(str).tolist()))
    if str(start_chrom) not in chrom_order or str(end_chrom) not in chrom_order:
        raise ValueError(f"Region chromosomes {start_chrom}, {end_chrom} are not present in pipeline_output")
    start_idx = chrom_order.index(str(start_chrom))
    end_idx = chrom_order.index(str(end_chrom))
    if start_idx > end_idx:
        raise ValueError("Region start chromosome appears after end chromosome in pipeline_output")
    order_map = {chrom: i for i, chrom in enumerate(chrom_order)}
    chrom_idx = df["chromosome"].astype(str).map(order_map)
    mask = (chrom_idx > start_idx) & (chrom_idx < end_idx)
    mask |= (
        (df["chromosome"].astype(str) == str(start_chrom))
        & (df["end_pos"] >= start_pos)
    )
    mask |= (
        (df["chromosome"].astype(str) == str(end_chrom))
        & (df["start_pos"] <= end_pos)
    )
    return df.loc[mask].reset_index(drop=True)


def _prepare_mii_assignment_plot_table(
    pipeline_output: Any,
    region: Optional[Union[str, List[Tuple[str, int]], Region, Tuple[Any, ...]]] = None,
    require_counts: bool = True,
) -> Tuple[pd.DataFrame, bool, str]:
    df = _normalize_mii_final_assignment_table(
        _load_mii_final_assignment_output(pipeline_output),
        require_counts=require_counts,
    )
    df = _filter_mii_final_assignment_region(df, region)
    if len(df) == 0:
        raise ValueError("No pipeline_output rows overlap the requested region")

    df = _sort_chromosome_position_table(df)
    single_chrom = df["chromosome"].astype(str).nunique() == 1
    if single_chrom:
        df["x_plot"] = df["position"]
        x_title = "genomic position (bp)"
    else:
        df["x_plot"] = np.arange(len(df))
        x_title = "genome-wide bin index"
    return df, single_chrom, x_title


def _add_mii_chrom_boundaries(fig, df: pd.DataFrame, single_chrom: bool, annotation_y: float = 1.08):
    if single_chrom:
        return fig

    chrom_starts = df.groupby("chromosome", sort=False)["x_plot"].min()
    chrom_ends = df.groupby("chromosome", sort=False)["x_plot"].max()
    for start in chrom_starts.iloc[1:]:
        fig.add_vline(x=float(start), line_width=1, line_color="black", opacity=0.5)
    for chrom in chrom_starts.index:
        mid = (chrom_starts[chrom] + chrom_ends[chrom]) / 2
        fig.add_annotation(
            x=mid,
            y=annotation_y,
            xref="x",
            yref="paper",
            text=str(chrom),
            showarrow=False,
            font=dict(size=18),
        )
    return fig


def plot_mii_assignment_track(
    pipeline_output: Any,
    region: Optional[Union[str, List[Tuple[str, int]], Region, Tuple[Any, ...]]] = None,
    output_path: Optional[Union[str, Path]] = None,
    title: str = "mII_final_assignment",
    width: int = 1400,
    height: int = 360,
):
    """
    Plot the mII final assignment label track.

    Parameters
    ----------
    pipeline_output
        Pipeline output as a DataFrame, tabular file, output directory, pipeline
        result dictionary, or SisterChromatidResults. The resolved table must
        contain chromosome, state, and position or start_pos/end_pos columns.
    region
        Region to plot using hic_basic standard representations, such as
        "chr1:1000000-2000000", [(chrom, start), (chrom, end)], or Region.
        If None, all rows are shown.
    output_path
        Optional PNG path. If provided, the Plotly figure is written to this path.
    title
        Figure title.
    width, height
        PNG export size in pixels when output_path is provided.

    Returns
    -------
    plotly.graph_objects.Figure
        The generated label track figure.
    """
    import plotly.graph_objects as go

    df, single_chrom, x_title = _prepare_mii_assignment_plot_table(
        pipeline_output,
        region=region,
        require_counts=False,
    )
    state_order = STATE_ORDER + sorted(set(df["state"].dropna()) - set(STATE_ORDER))
    state_to_y = {state: i for i, state in enumerate(state_order)}
    df["label_y"] = df["state"].map(state_to_y)

    fig = go.Figure()
    hover_base = (
        "chrom=%{customdata[0]}<br>"
        "start=%{customdata[1]}<br>"
        "end=%{customdata[2]}<br>"
    )
    for state in state_order:
        mask = df["state"] == state
        if mask.any():
            fig.add_trace(
                go.Scattergl(
                    x=df.loc[mask, "x_plot"],
                    y=df.loc[mask, "label_y"],
                    mode="markers",
                    marker=dict(color=STATE_COLORS.get(state, "#333333"), size=8),
                    name=f"{state} ({int(mask.sum())})",
                    customdata=df.loc[mask, ["chromosome", "start_pos", "end_pos"]].to_numpy(),
                    hovertemplate=hover_base + f"state={state}<extra></extra>",
                )
            )

    _add_mii_chrom_boundaries(fig, df, single_chrom)
    fig.update_yaxes(
        tickmode="array",
        tickvals=list(range(len(state_order))),
        ticktext=state_order,
        title=dict(text="label", font=dict(size=20)),
        tickfont=dict(size=18),
        range=[-0.5, len(state_order) - 0.5],
    )
    fig.update_xaxes(title=dict(text=x_title, font=dict(size=20)), tickfont=dict(size=18))
    fig.update_xaxes(showline=True, linewidth=1.4, linecolor="black", mirror=True, showgrid=False, zeroline=False, ticks="outside")
    fig.update_yaxes(showline=True, linewidth=1.4, linecolor="black", mirror=True, showgrid=False, zeroline=False, ticks="outside")
    fig.update_layout(
        title=title,
        title_font=dict(size=28),
        font=dict(size=18),
        width=width,
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=-0.32, font=dict(size=16)),
    )

    if output_path is not None:
        _write_plotly_png(fig, output_path, width=width, height=height)
    return fig


def plot_mii_final_assignment(
    pipeline_output: Any,
    region: Optional[Union[str, List[Tuple[str, int]], Region, Tuple[Any, ...]]] = None,
    output_path: Optional[Union[str, Path]] = None,
    title: str = "Final assigned labels: genome_wide_label + genome_wide_count",
    width: int = 1400,
    height: int = 760,
):
    """
    Plot final mII HMM labels together with ref/alt/other counts.

    Parameters
    ----------
    pipeline_output
        Pipeline output as a DataFrame, tabular file, output directory, pipeline
        result dictionary. The resolved table must contain chromosome, position
        (or start_pos/end_pos), state, ref, alt, and other columns.
    region
        Region to plot using hic_basic standard representations, such as
        "chr1:1000000-2000000", [(chrom, start), (chrom, end)], or Region.
        If None, all rows are shown.
    output_path
        Optional PNG path. If provided, the Plotly figure is written to this path.
    title
        Figure title.
    width, height
        PNG export size in pixels when output_path is provided.

    Returns
    -------
    plotly.graph_objects.Figure
        The generated Plotly figure.
    """
    from plotly.subplots import make_subplots

    df, single_chrom, x_title = _prepare_mii_assignment_plot_table(
        pipeline_output,
        region=region,
        require_counts=True,
    )
    label_fig = plot_mii_assignment_track(
        df,
        region=None,
        title="final_label",
        width=width,
        height=max(int(height * 0.42), 300),
    )
    count_fig = plot_genome_phased_read_count_track(
        df,
        region=None,
        title="ref_alt_other_count",
        width=width,
        height=max(int(height * 0.58), 360),
    )
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        row_heights=[0.42, 0.58],
        vertical_spacing=0.08,
        subplot_titles=("final_label", "ref_alt_other_count"),
    )

    for trace in label_fig.data:
        fig.add_trace(trace, row=1, col=1)
    for trace in count_fig.data:
        fig.add_trace(trace, row=2, col=1)

    _add_mii_chrom_boundaries(fig, df, single_chrom)
    fig.update_yaxes(label_fig.layout.yaxis.to_plotly_json(), row=1, col=1)
    fig.update_yaxes(count_fig.layout.yaxis.to_plotly_json(), row=2, col=1)
    fig.update_xaxes(title=dict(text=x_title, font=dict(size=20)), tickfont=dict(size=18), row=2, col=1)
    fig.update_xaxes(showline=True, linewidth=1.4, linecolor="black", mirror=True, showgrid=False, zeroline=False, ticks="outside")
    fig.update_yaxes(showline=True, linewidth=1.4, linecolor="black", mirror=True, showgrid=False, zeroline=False, ticks="outside")

    for annotation in count_fig.layout.annotations:
        if annotation.text and "count y-axis capped" in annotation.text:
            ann = annotation.to_plotly_json()
            ann["y"] = 0.48
            ann["font"] = dict(size=17)
            fig.add_annotation(**ann)
    fig.update_layout(
        title=title,
        title_font=dict(size=30),
        font=dict(size=19),
        width=width,
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=-0.28, font=dict(size=17)),
    )

    if output_path is not None:
        _write_plotly_png(fig, output_path, width=width, height=height)
    return fig


def _plot_assigned_label_count_png(
    counts_df: pd.DataFrame,
    assigned_labels: Union[pd.Series, List[str]],
    output_path: Union[str, Path],
    title: str = "Final assigned labels: genome_wide_label + genome_wide_count",
) -> None:
    """Plot final assigned HMM labels through the public Plotly implementation."""
    df = counts_df.copy().reset_index(drop=True)
    df["state"] = pd.Series(assigned_labels).reset_index(drop=True)
    plot_mii_final_assignment(
        pipeline_output=df,
        region=None,
        output_path=output_path,
        title=title,
    )


def _plot_rough_decision_distributions_png(
    counts_df: pd.DataFrame,
    output_path: Union[str, Path],
    effective_depth_threshold: Union[pd.Series, float, int],
    other_threshold: float,
    balance_threshold: float,
) -> None:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    df = counts_df.copy().reset_index(drop=True)
    total = df["ref"] + df["alt"] + df["other"]
    total_log = np.log10(total + 1)
    safe_total = total.replace(0, np.nan)
    other_frac = (df["other"] / safe_total).fillna(0.0)
    ref_frac = (df["ref"] / safe_total).fillna(0.0)
    alt_frac = (df["alt"] / safe_total).fillna(0.0)
    ref_alt_delta = (ref_frac - alt_frac).abs()
    nonzero_other_frac = other_frac[other_frac > 0]
    if len(nonzero_other_frac) > 0:
        other_frac_eps = min(1e-4, float(nonzero_other_frac.min()) / 2)
    else:
        other_frac_eps = 1e-4
    other_frac_log = np.log10(other_frac + other_frac_eps)
    other_threshold_log = float(np.log10(other_threshold + other_frac_eps))

    if isinstance(effective_depth_threshold, pd.Series):
        depth_cutoff = effective_depth_threshold.reset_index(drop=True).astype(float)
    else:
        depth_cutoff = pd.Series(float(effective_depth_threshold), index=df.index)
    pass_depth = total >= depth_cutoff
    pass_other = pass_depth & (other_frac <= other_threshold)
    constant_depth_cutoff = depth_cutoff.nunique(dropna=False) == 1

    fig = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=(
            "1. log10(total depth + 1) (all bins)",
            "2. log10(other fraction + eps) (after depth filter)",
            "3. ref/alt delta (after depth+other filters)",
        ),
    )

    fig.add_trace(
        go.Histogram(
            x=total_log,
            nbinsx=50,
            marker=dict(color="#4c78a8"),
            name="total depth",
            showlegend=False,
        ),
        row=1,
        col=1,
    )
    if constant_depth_cutoff:
        cutoff_label = f"{float(depth_cutoff.iloc[0]):.0f}"
        fig.add_vline(
            x=float(np.log10(depth_cutoff.iloc[0] + 1)),
            line=dict(color="#d62728", width=3, dash="dash"),
            annotation_text=f"effective_depth_threshold={cutoff_label}",
            annotation_font_size=16,
            row=1,
            col=1,
        )
    else:
        median_cutoff = float(depth_cutoff.median())
        fig.add_vline(
            x=float(np.log10(median_cutoff + 1)),
            line=dict(color="#d62728", width=3, dash="dash"),
            annotation_text=f"median effective_depth_threshold={median_cutoff:.0f}",
            annotation_font_size=16,
            row=1,
            col=1,
        )

    fig.add_trace(
        go.Histogram(
            x=other_frac_log[pass_depth],
            nbinsx=50,
            marker=dict(color="#7f7f7f"),
            name="other fraction",
            showlegend=False,
        ),
        row=1,
        col=2,
    )
    fig.add_vline(
        x=other_threshold_log,
        line=dict(color="#d62728", width=3, dash="dash"),
        annotation_text=f"log10(other_threshold+eps)={other_threshold_log:.2f}",
        annotation_font_size=16,
        row=1,
        col=2,
    )

    fig.add_trace(
        go.Histogram(
            x=ref_alt_delta[pass_other],
            nbinsx=50,
            marker=dict(color="#2ca02c"),
            name="ref/alt delta",
            showlegend=False,
        ),
        row=1,
        col=3,
    )
    fig.add_vline(
        x=balance_threshold,
        line=dict(color="#d62728", width=3, dash="dash"),
        annotation_text=f"balance_threshold={balance_threshold:.2f}",
        annotation_font_size=16,
        row=1,
        col=3,
    )

    pass_depth_n = int(pass_depth.sum())
    pass_other_n = int(pass_other.sum())
    fig.add_annotation(
        text=(
            f"bins={len(df):,}; pass depth={pass_depth_n:,}; "
            f"pass depth+other={pass_other_n:,}; "
            f"other fraction eps={other_frac_eps:.2g}; zeros={int((other_frac == 0).sum()):,}"
        ),
        x=0.5,
        y=-0.18,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(size=18),
    )

    fig.update_xaxes(
        title=dict(text="log10(total depth + 1)", font=dict(size=17)),
        tickfont=dict(size=15),
        row=1,
        col=1,
    )
    fig.update_xaxes(
        title=dict(text="log10(other fraction + eps)", font=dict(size=17)),
        tickfont=dict(size=15),
        row=1,
        col=2,
    )
    fig.update_yaxes(title=dict(text="bin count", font=dict(size=18)), tickfont=dict(size=15))
    fig.update_layout(
        template="plotly_white",
        title="Rough classification decision distributions",
        title_font=dict(size=28),
        font=dict(size=18),
        bargap=0.05,
        height=560,
        width=1600,
        margin=dict(b=110),
    )
    _write_plotly_png(fig, output_path, width=1600, height=560)


def _plot_rough_segment_lengths_png(
    segment_table: pd.DataFrame,
    min_segment_bins: int,
    output_path: Union[str, Path],
) -> None:
    """Plot raw rough-label segment length distributions by state."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    n_cols = 2
    n_rows = int(np.ceil(len(ROUGH_LABEL_ORDER) / n_cols))
    fig = make_subplots(
        rows=n_rows,
        cols=2,
        subplot_titles=ROUGH_LABEL_ORDER,
        horizontal_spacing=0.08,
        vertical_spacing=0.14,
    )
    if not segment_table.empty:
        max_len = int(segment_table["length_bins"].max())
        nbinsx = max(10, min(60, max_len))
        for state_idx, state in enumerate(ROUGH_LABEL_ORDER):
            state_lengths = segment_table.loc[
                segment_table["label"] == state,
                "length_bins",
            ]
            row = state_idx // n_cols + 1
            col = state_idx % n_cols + 1
            fig.add_trace(
                go.Histogram(
                    x=state_lengths,
                    name=f"{state} ({len(state_lengths)})",
                    marker=dict(color=STATE_COLORS.get(state, "#333333")),
                    opacity=0.8,
                    nbinsx=nbinsx,
                    showlegend=False,
                ),
                row=row,
                col=col,
            )
            fig.add_vline(
                x=min_segment_bins,
                line=dict(color="#d62728", width=3, dash="dash"),
                annotation_text=f"threshold={min_segment_bins}",
                annotation_font_size=15,
                row=row,
                col=col,
            )
            fig.update_xaxes(title=dict(text="length (bins)", font=dict(size=15)), tickfont=dict(size=14), row=row, col=col)
            fig.update_yaxes(title=dict(text="segments", font=dict(size=15)), tickfont=dict(size=14), row=row, col=col)

    short_n = int((segment_table["length_bins"] < min_segment_bins).sum()) if not segment_table.empty else 0
    total_n = int(len(segment_table))
    short_pct = short_n / total_n if total_n > 0 else 0.0
    fig.add_annotation(
        text=(
            f"transition_min_segment_bins={min_segment_bins}; "
            f"raw segments={total_n:,}; filtered as gaps={short_n:,} ({short_pct:.1%})"
        ),
        x=0.5,
        y=-0.12,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(size=18),
    )
    fig.update_layout(
        template="plotly_white",
        title="Raw rough label segment length distribution",
        title_font=dict(size=28),
        font=dict(size=18),
        bargap=0.05,
        width=1300,
        height=980,
        margin=dict(b=120),
    )
    _write_plotly_png(fig, output_path, width=1300, height=980)


def _plot_raw_filtered_label_png(
    counts_df: pd.DataFrame,
    raw_labels: pd.Series,
    filtered_labels: pd.Series,
    output_path: Union[str, Path],
) -> None:
    """Plot genome-wide raw rough labels and transition-filtered labels."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    df = counts_df.copy().reset_index(drop=True)
    df["raw_label"] = raw_labels.reset_index(drop=True)
    df["filtered_label"] = filtered_labels.reset_index(drop=True)
    df["x_id"] = np.arange(len(df))
    label_order = ["MM", "PM", "PP", "U", AMBIG_LABEL, GAP_LABEL]
    label_to_y = {label: i for i, label in enumerate(label_order)}
    label_names = {
        "MM": "MM",
        "PM": "PM",
        "PP": "PP",
        "U": "U",
        AMBIG_LABEL: "AMBIG",
        GAP_LABEL: "gap",
    }

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        row_heights=[0.5, 0.5],
        vertical_spacing=0.12,
        subplot_titles=("genome_wide_raw_rough_label", "genome_wide_filtered_transition_label"),
    )

    for row_idx, label_col in enumerate(["raw_label", "filtered_label"], start=1):
        for label in label_order:
            mask = df[label_col] == label
            if not mask.any():
                continue
            color = "#bdbdbd" if label == GAP_LABEL else STATE_COLORS.get(label, "#333333")
            fig.add_trace(
                go.Scattergl(
                    x=df.loc[mask, "x_id"],
                    y=np.full(int(mask.sum()), label_to_y[label]),
                    mode="markers",
                    marker=dict(color=color, size=7),
                    name=f"{label_names[label]} ({int(mask.sum())})",
                    showlegend=row_idx == 2,
                ),
                row=row_idx,
                col=1,
            )

    chrom_starts = df.groupby("chromosome", sort=False)["x_id"].min()
    chrom_ends = df.groupby("chromosome", sort=False)["x_id"].max()
    for start in chrom_starts.iloc[1:]:
        fig.add_vline(x=float(start), line_width=1, line_color="black", opacity=0.45)
    for chrom in chrom_starts.index:
        mid = (chrom_starts[chrom] + chrom_ends[chrom]) / 2
        fig.add_annotation(
            x=mid,
            y=1.08,
            xref="x",
            yref="paper",
            text=str(chrom),
            showarrow=False,
            font=dict(size=15),
        )

    for row_idx in [1, 2]:
        fig.update_yaxes(
            tickmode="array",
            tickvals=list(range(len(label_order))),
            ticktext=[label_names[label] for label in label_order],
            title=dict(text="label", font=dict(size=18)),
            tickfont=dict(size=16),
            row=row_idx,
            col=1,
        )
    fig.update_xaxes(title=dict(text="genome-wide bin index", font=dict(size=18)), tickfont=dict(size=16), row=2, col=1)
    fig.update_layout(
        template="plotly_white",
        title="Raw rough labels vs transition-filtered labels",
        title_font=dict(size=28),
        font=dict(size=18),
        height=780,
        width=1400,
        legend=dict(orientation="h", yanchor="bottom", y=-0.22, font=dict(size=16)),
    )
    _write_plotly_png(fig, output_path, width=1400, height=780)


def _print_estimation_diagnostics(
    counts_df: pd.DataFrame,
    rough_labels: Optional[pd.Series],
    emission_source: str,
    emission_table: pd.DataFrame,
    transition_source: str,
    transition_param_dict: Dict[str, float],
    rough_segment_table: Optional[pd.DataFrame] = None,
    transition_min_segment_bins: int = 3,
) -> None:
    total_depth = counts_df[["ref", "alt", "other"]].sum(axis=1)
    print(
        "[mII estimate] units={:,}, mean_depth={:.3f}, median_depth={:.3f}".format(
            len(counts_df), float(total_depth.mean()), float(total_depth.median())
        )
    )
    if rough_labels is None:
        print("[mII estimate] rough labels skipped: emission and transition parameters were specified")
    else:
        label_counts = rough_labels.value_counts().reindex(ROUGH_LABEL_ORDER, fill_value=0)
        label_summary = ", ".join(
            f"{label}={int(count)} ({count / len(rough_labels):.1%})"
            for label, count in label_counts.items()
        )
        n_switches = int((rough_labels.shift(1) != rough_labels).iloc[1:].sum())
        print(f"[mII estimate] rough labels: {label_summary}; adjacent_label_switches={n_switches}")
        if rough_segment_table is not None:
            short_n = int((rough_segment_table["length_bins"] < transition_min_segment_bins).sum())
            print(
                "[mII estimate] raw rough segments={:,}; short segments (<{} bins) skipped for transition={:,}".format(
                    len(rough_segment_table), transition_min_segment_bins, short_n
                )
            )

    print(f"[mII estimate] emission source={emission_source}")
    for _, row in emission_table.iterrows():
        print(
            "[mII estimate] emission {label}: ref_ep={ref_ep:.4f}, alt_ep={alt_ep:.4f}, other_ep={other_ep:.4f}".format(
                label=row["label"],
                ref_ep=float(row["ref_ep"]),
                alt_ep=float(row["alt_ep"]),
                other_ep=float(row["other_ep"]),
            )
        )

    print(
        "[mII estimate] transition source={}; alpha={:.6f}, beta={:.6f}, gamma={:.6f}, rho={:.6f}".format(
            transition_source,
            transition_param_dict["alpha"],
            transition_param_dict["beta"],
            transition_param_dict["gamma"],
            transition_param_dict["rho"],
        )
    )


def estimate_mii_hmm_parameters(
    read_counts_df: pd.DataFrame,
    emission_params: Optional[Dict[str, Tuple[float, float, float]]] = None,
    transition_params: Union[str, Dict[str, float]] = "default",
    alpha: float = 0.01,
    beta: float = 0.02,
    gamma: float = 0.001,
    rho: float = 0.1,
    depth_threshold: int = 2,
    other_threshold: float = 0.4,
    balance_threshold: float = 0.25,
    binning_method: Optional[str] = None,
    binning_snp_number: int = 250,
    transition_param_unit: str = "step",
    transition_min_segment_bins: int = 3,
    transition_u_gamma_prior: Optional[float] = None,
    transition_u_gamma_factor: float = 0.1,
    transition_u_rho_prior: float = 0.32,
    u_emission_other_prior: float = 0.05,
    params_vis_prefix: Optional[Union[str, Path]] = None,
    verbose: int = 0,
) -> Dict[str, Any]:
    """
    Estimate or resolve mII HMM parameters from read-count data.

    Args:
        read_counts_df: DataFrame with chromosome, position, ref, alt, other.
        emission_params: Optional specified emission probabilities by label.
        transition_params: "default", "estimated", or a dict with alpha/beta/gamma/rho.
        alpha, beta, gamma, rho: Fixed transition parameters used for "default".
        depth_threshold: Minimum depth for rough classification.
        other_threshold: Maximum other fraction for rough classification.
        balance_threshold: Ref/alt balance threshold for rough classification.
        binning_method: None or "by_snp".
        binning_snp_number: Number of SNPs per bin when binning_method="by_snp".
        transition_param_unit: Unit for caller-specified transition parameters:
            "step" for per HMM observation, "mb" for per-Mb rates, or "snp" for
            per-SNP probabilities. Estimated parameters are always per-step.
        transition_min_segment_bins: Raw rough-label segments shorter than this
            many bins are marked as gaps and skipped for transition estimation.
        transition_u_gamma_prior: Optional fixed gamma prior used for U entry
            when binning_method="by_snp" and transition_params="estimated". If
            None, uses transition_u_gamma_factor * min(alpha, beta).
        transition_u_gamma_factor: Factor for deriving the by-SNP U-entry prior
            from estimated alpha/beta when transition_u_gamma_prior is None.
        transition_u_rho_prior: By-SNP prior for U -> each biological state.
        u_emission_other_prior: By-SNP fixed U-emission other probability
            used when emission_params is None. Ref and alt share the remaining
            probability equally.
        params_vis_prefix: Prefix for diagnostic PNG files. When verbose >= 1,
            writes rough_label_classification, rough_decision_distributions,
            rough_segment_lengths, rough_label_filtering, emission_matrix, and
            transition_matrix PNGs. During full decoding it also writes
            final_assigned_label_count PNGs.
        verbose: If >=1, print diagnostics and, when params_vis_prefix is set, write PNGs.

    Returns:
        Dictionary containing binned counts, rough labels, emission parameters,
        transition parameters, transition matrix, and diagnostic tables.
    """
    counts_df, original_counts_df = _prepare_counts_for_estimation(
        read_counts_df,
        binning_method=binning_method,
        binning_snp_number=binning_snp_number,
    )
    effective_depth_threshold = _effective_depth_threshold(
        counts_df,
        depth_threshold=depth_threshold,
        binning_method=binning_method,
    )
    median_distance_mb = _median_observation_distance_mb(counts_df)
    median_distance_snps = _median_observation_distance_snps(
        counts_df,
        binning_method=binning_method,
    )
    if transition_param_unit not in {"step", "mb", "snp"}:
        raise ValueError("transition_param_unit must be one of ['mb', 'snp', 'step']")
    if transition_min_segment_bins < 1:
        raise ValueError("transition_min_segment_bins must be >= 1")
    u_emission_prior = _u_emission_prior_from_other(u_emission_other_prior)

    transition_needs_rough = transition_params == "estimated"
    emission_needs_rough = emission_params is None
    rough_labels: Optional[pd.Series] = None
    rough_segment_table: Optional[pd.DataFrame] = None
    transition_rough_labels: Optional[pd.Series] = None
    filtered_rough_labels: Optional[pd.Series] = None

    if emission_needs_rough or transition_needs_rough:
        rough_labels = _rough_classification_with_thresholds(
            counts_df,
            effective_depth_threshold=effective_depth_threshold,
            other_threshold=other_threshold,
            balance_threshold=balance_threshold,
        ).reset_index(drop=True)
        rough_segment_table = _rough_label_segments(counts_df, rough_labels)

    if transition_needs_rough:
        filtered_rough_labels = _filter_short_rough_segments(
            rough_labels,
            rough_segment_table,
            min_segment_bins=transition_min_segment_bins,
        )
        filtered_rough_labels = filtered_rough_labels.mask(filtered_rough_labels == AMBIG_LABEL, GAP_LABEL)
        if binning_method == "by_snp":
            filtered_rough_labels = filtered_rough_labels.mask(filtered_rough_labels == "U", GAP_LABEL)
        transition_rough_labels = _labels_with_chromosome_gap_separators(
            counts_df,
            filtered_rough_labels,
        )

    if emission_params is None:
        resolved_emission_params = estimate_emissions_from_counts(counts_df, rough_labels)
        emission_source = "estimated"
        if binning_method == "by_snp":
            resolved_emission_params["U"] = u_emission_prior
            emission_source = "estimated_with_by_snp_u_prior"
    else:
        missing = set(STATE_ORDER) - set(emission_params)
        if missing:
            raise ValueError(f"emission_params is missing labels: {sorted(missing)}")
        resolved_emission_params = {
            state: tuple(float(x) for x in emission_params[state])
            for state in STATE_ORDER
        }
        emission_source = "specified"

    transition_param_dict, transition_source, transition_scale_info = _transition_params_from_source(
        transition_params=transition_params,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        rho=rho,
        rough_labels_for_transition=transition_rough_labels,
        transition_param_unit=transition_param_unit,
        median_distance_mb=median_distance_mb,
        median_distance_snps=median_distance_snps,
    )
    if transition_needs_rough and rough_segment_table is not None:
        short_segments = rough_segment_table["length_bins"] < transition_min_segment_bins
        transition_scale_info.update({
            "transition_min_segment_bins": transition_min_segment_bins,
            "raw_rough_segments": int(len(rough_segment_table)),
            "short_segments_skipped": int(short_segments.sum()),
            "short_segment_bins_skipped": int(rough_segment_table.loc[short_segments, "length_bins"].sum()),
        })
        if binning_method == "by_snp":
            _apply_by_snp_u_transition_prior(
                transition_param_dict=transition_param_dict,
                transition_scale_info=transition_scale_info,
                transition_u_gamma_prior=transition_u_gamma_prior,
                transition_u_gamma_factor=transition_u_gamma_factor,
                transition_u_rho_prior=transition_u_rho_prior,
            )
    hmm = SisterChromatidHMM(
        alpha=transition_param_dict["alpha"],
        beta=transition_param_dict["beta"],
        gamma=transition_param_dict["gamma"],
        rho=transition_param_dict["rho"],
        emission_params=resolved_emission_params,
    )
    emission_table = _emission_stats_table(
        counts_df=counts_df,
        rough_labels=rough_labels,
        emission_params=resolved_emission_params,
        source=emission_source,
    )

    if verbose >= 1:
        _print_estimation_diagnostics(
            counts_df=counts_df,
            rough_labels=rough_labels,
            emission_source=emission_source,
            emission_table=emission_table,
            transition_source=transition_source,
            transition_param_dict=transition_param_dict,
            rough_segment_table=rough_segment_table,
            transition_min_segment_bins=transition_min_segment_bins,
        )

    if verbose >= 1 and params_vis_prefix is not None:
        both_specified = emission_source == "specified" and transition_source == "specified"
        if both_specified:
            _blank_diagnostic_png(
                "Both emission and transition parameters were specified; rough label classification was skipped.",
                _diagnostic_path(params_vis_prefix, "rough_label_classification"),
            )
        else:
            _plot_rough_label_count_png(
                counts_df=counts_df,
                rough_labels=rough_labels,
                output_path=_diagnostic_path(params_vis_prefix, "rough_label_classification"),
                title="Rough label classification: genome_wide_label + genome_wide_count",
            )
        _plot_rough_decision_distributions_png(
            counts_df=counts_df,
            output_path=_diagnostic_path(params_vis_prefix, "rough_decision_distributions"),
            effective_depth_threshold=effective_depth_threshold,
            other_threshold=other_threshold,
            balance_threshold=balance_threshold,
        )
        if rough_segment_table is None:
            _blank_diagnostic_png(
                "Both emission and transition parameters were specified; raw rough label segment diagnostics were skipped.",
                _diagnostic_path(params_vis_prefix, "rough_segment_lengths"),
            )
            _blank_diagnostic_png(
                "Both emission and transition parameters were specified; rough label filtering diagnostics were skipped.",
                _diagnostic_path(params_vis_prefix, "rough_label_filtering"),
            )
        else:
            _plot_rough_segment_lengths_png(
                segment_table=rough_segment_table,
                min_segment_bins=transition_min_segment_bins,
                output_path=_diagnostic_path(params_vis_prefix, "rough_segment_lengths"),
            )
            if filtered_rough_labels is None:
                _blank_diagnostic_png(
                    "Transition parameters were not estimated; rough label filtering was skipped.",
                    _diagnostic_path(params_vis_prefix, "rough_label_filtering"),
                )
            else:
                _plot_raw_filtered_label_png(
                    counts_df=counts_df,
                    raw_labels=rough_labels,
                    filtered_labels=filtered_rough_labels,
                    output_path=_diagnostic_path(params_vis_prefix, "rough_label_filtering"),
                )
        if emission_source == "estimated_with_by_snp_u_prior":
            emission_title = "Emission matrix estimated from rough labels with fixed by-SNP U prior"
        elif emission_source == "estimated":
            emission_title = "Emission matrix estimated from rough labels"
        else:
            emission_title = "Emission matrix specified by caller"
        _plot_table_png(
            emission_table,
            emission_title,
            _diagnostic_path(params_vis_prefix, "emission_matrix"),
        )
        _plot_transition_matrix_png(
            transition_matrix=hmm.transition_matrix,
            transition_param_dict=transition_param_dict,
            source=transition_source,
            output_path=_diagnostic_path(params_vis_prefix, "transition_matrix"),
        )

    return {
        "counts_df": counts_df,
        "original_counts_df": original_counts_df,
        "rough_labels": rough_labels,
        "transition_rough_labels": transition_rough_labels,
        "filtered_rough_labels": filtered_rough_labels,
        "rough_label_segments": rough_segment_table,
        "emission_params": resolved_emission_params,
        "emission_source": emission_source,
        "emission_table": emission_table,
        "transition_params": transition_param_dict,
        "transition_source": transition_source,
        "transition_param_unit": transition_param_unit,
        "transition_scale_info": transition_scale_info,
        "transition_matrix": hmm.transition_matrix,
        "effective_depth_threshold": effective_depth_threshold,
        "binning_method": binning_method,
        "binning_snp_number": binning_snp_number,
        "transition_min_segment_bins": transition_min_segment_bins,
        "transition_u_gamma_prior": transition_u_gamma_prior,
        "transition_u_gamma_factor": transition_u_gamma_factor,
        "transition_u_rho_prior": transition_u_rho_prior,
        "u_emission_other_prior": u_emission_other_prior,
        "u_emission_prior_applied": emission_params is None and binning_method == "by_snp",
    }


def convert_unit_coordinates_to_snp(unit_results: SisterChromatidResults, 
                                     original_read_counts: pd.DataFrame,
                                     n_snps: int = 10) -> SisterChromatidResults:
    """
    Convert unit-level state assignments back to original SNP coordinates.
    
    Args:
        unit_results: Results from HMM at unit level
        original_read_counts: Original SNP-level read counts
        n_snps: Number of SNPs per unit
    
    Returns:
        SisterChromatidResults with SNP-level states
    """
    # Normalize column names
    if 'pos' in original_read_counts.columns and 'position' not in original_read_counts.columns:
        original_read_counts = original_read_counts.copy()
        original_read_counts['position'] = original_read_counts['pos']
    
    # Expand states from units to SNPs
    snp_states = []
    snp_positions = []
    
    for i, state in enumerate(unit_results.states):
        # Get SNPs in this unit
        start_idx = i * n_snps
        end_idx = min((i + 1) * n_snps, len(original_read_counts))
        
        for j in range(start_idx, end_idx):
            if j < len(original_read_counts):
                snp_states.append(state)
                snp_positions.append(original_read_counts.iloc[j]['position'])
    
    # Recalculate blocks at SNP level
    snp_blocks = []
    n_snps_total = len(snp_states)
    if n_snps_total > 0:
        i = 0
        while i < n_snps_total:
            current_state = snp_states[i]
            start_pos = snp_positions[i]
            start_idx = i
            
            while i < n_snps_total and snp_states[i] == current_state:
                i += 1
            
            end_pos = snp_positions[i - 1]
            num_snps = i - start_idx
            
            block = SisterChromatidBlock(
                chromosome=unit_results.chromosome,
                start_pos=start_pos,
                end_pos=end_pos,
                state=current_state,
                confidence=1.0,
                num_snps=num_snps
            )
            snp_blocks.append(block)
    
    # Create new results object
    snp_results = SisterChromatidResults(
        sample_id=unit_results.sample_id,
        chromosome=unit_results.chromosome,
        states=snp_states,
        positions=snp_positions,
        blocks=snp_blocks,
        statistics=unit_results.statistics
    )
    
    return snp_results


def extract_read_counts(
    bam_file: str,
    snp_file: str,
    output_file: Optional[str] = None,
    min_mapq: int = 20,
    min_baseq: int = 20
) -> pd.DataFrame:
    """
    Extract read counts at SNP positions from BAM file.

    Args:
        bam_file: Path to BAM/BCF file
        snp_file: Path to SNP file (VCF or TSV with positions)
        output_file: Optional output file path
        min_mapq: Minimum mapping quality
        min_baseq: Minimum base quality

    Returns:
        DataFrame with columns: chromosome, position, ref, alt, other
    """
    import pysam

    logger.info(f"Extracting read counts from {bam_file}")

    # Load SNP positions
    if snp_file.endswith('.vcf') or snp_file.endswith('.vcf.gz'):
        snps = []
        with pysam.VariantFile(snp_file) as vcf:
            for record in vcf:
                snps.append({
                    'chromosome': record.chrom,
                    'position': record.pos,
                    'ref': record.ref,
                    'alt': record.alts[0] if record.alts else '.',
                })
        snp_df = pd.DataFrame(snps)
    else:
        snp_df = pd.read_csv(snp_file, sep='\t')

    # Count reads at each SNP
    counts = []
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for _, row in snp_df.iterrows():
            chrom = row['chromosome']
            pos = row['position']

            ref_count = 0
            alt_count = 0
            other_count = 0

            try:
                for pileupcolumn in bam.pileup(
                    chrom,
                    start=pos - 1,
                    end=pos,
                    stepper='samtools',
                    min_mapping_quality=min_mapq
                ):
                    if pileupcolumn.pos != pos - 1:
                        continue

                    for pileupread in pileupcolumn.pileups:
                        if pileupread.is_del or pileupread.is_refskip:
                            continue

                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        qual = pileupread.alignment.query_qualities[pileupread.query_position]

                        if qual < min_baseq:
                            continue

                        ref_allele = row.get('ref', '')
                        alt_allele = row.get('alt', '')

                        if len(ref_allele) == 1 and len(alt_allele) == 1:
                            if base.upper() == ref_allele.upper():
                                ref_count += 1
                            elif base.upper() == alt_allele.upper():
                                alt_count += 1
                            else:
                                other_count += 1
            except ValueError:
                # Chromosome not found in BAM
                pass

            counts.append({
                'chromosome': chrom,
                'position': pos,
                'ref': ref_count,
                'alt': alt_count,
                'other': other_count,
            })

    counts_df = pd.DataFrame(counts)

    if output_file:
        counts_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved read counts to {output_file}")

    return counts_df


def _process_sample_mii(
    sample_id: str,
    bam_file: str,
    snp_file: str,
    output_dir: str,
    cache_dir: str,
    chromosome: Optional[str] = None,
    transition_params: Union[str, Dict[str, float]] = "default",  # "default" | "estimated" | dict
    emission_params: Optional[Dict[str, Tuple[float, float, float]]] = None,
    alpha: float = 0.01,
    beta: float = 0.02,
    gamma: float = 0.001,
    rho: float = 0.1,
    depth_threshold: int = 2,
    other_threshold: float = 0.4,
    balance_threshold: float = 0.25,
    binning_method: Optional[str] = None,
    binning_snp_number: int = 250,
    transition_param_unit: str = "step",
    transition_min_segment_bins: int = 3,
    transition_u_gamma_prior: Optional[float] = None,
    transition_u_gamma_factor: float = 0.1,
    transition_u_rho_prior: float = 0.32,
    u_emission_other_prior: float = 0.05,
    estimate_only: bool = False,
    params_vis_prefix: Optional[Union[str, Path]] = None,
    verbose: int = 0,
) -> Dict[str, Any]:
    """
    Process a single sample for sister chromatid inference.

    Args:
        sample_id: Sample identifier
        bam_file: Path to BAM file
        snp_file: Path to SNP file
        output_dir: Output directory
        cache_dir: Cache directory
        chromosome: Optional chromosome to analyze
        transition_params: Transition parameter source ("default" | "estimated")
        alpha, beta, gamma, rho: Transition parameters (used when transition_params="default")
        depth_threshold: Minimum depth for confident call (default: 2)
        other_threshold: Maximum other fraction for confident call (default: 0.4)
        balance_threshold: Threshold for classifying as balanced (default: 0.25)
                          PP: ref_frac > 0.75, MM: alt_frac > 0.75
        binning_method: Binning method (None for single-SNP, "by_snp" for merging by SNP count)
        binning_snp_number: Number of SNPs per unit (default: 250)

    Returns:
        Dictionary with results and statistics
    """
    cache = CacheManager(cache_dir)

    # Compute cache key
    params = {
        "sample_id": sample_id,
        "bam_file": bam_file,
        "snp_file": snp_file,
        "chromosome": chromosome,
        "alpha": alpha,
        "beta": beta,
        "gamma": gamma,
        "rho": rho,
        "emission_params": emission_params,
        "transition_params": transition_params,
        "depth_threshold": depth_threshold,
        "other_threshold": other_threshold,
        "balance_threshold": balance_threshold,
        "binning_method": binning_method,
        "binning_snp_number": binning_snp_number,
        "transition_param_unit": transition_param_unit,
        "transition_min_segment_bins": transition_min_segment_bins,
        "transition_u_gamma_prior": transition_u_gamma_prior,
        "transition_u_gamma_factor": transition_u_gamma_factor,
        "transition_u_rho_prior": transition_u_rho_prior,
        "u_emission_other_prior": u_emission_other_prior,
        "estimate_only": estimate_only,
    }
    key = cache.compute_key(params)

    result = {
        "sample_id": sample_id,
        "success": False,
        "steps_completed": [],
        "output_files": {},
        "error": None,
        "results": None,
    }

    try:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Step 1: Extract read counts
        counts_file = output_path / f"{sample_id}_read_counts.tsv"
        if chromosome:
            counts_file = output_path / f"{sample_id}_{chromosome}_read_counts.tsv"

        if not cache.exists(key, "extract_counts"):
            logger.info(f"Step 1/4: Extracting read counts for {sample_id}")
            counts_df = extract_read_counts(bam_file, snp_file, str(counts_file))
            cache.save(counts_df, key, "extract_counts")
            cache.set_flag(key, "extract_counts", {"n_snps": len(counts_df)})
        else:
            logger.info("Step 1/4: Using cached read counts")
            counts_df = cache.load(key, "extract_counts")

        if chromosome:
            counts_df = counts_df[counts_df['chromosome'] == chromosome]

        if len(counts_df) == 0:
            raise ValueError("No SNPs found after filtering")

        result["steps_completed"].append("extract_counts")
        result["output_files"]["read_counts"] = str(counts_file)

        # Step 2-3: Binning, rough classification, and parameter estimation
        logger.info("Step 2/4: Estimating or resolving HMM parameters")
        params_dict = estimate_mii_hmm_parameters(
            counts_df,
            emission_params=emission_params,
            transition_params=transition_params,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            rho=rho,
            depth_threshold=depth_threshold,
            other_threshold=other_threshold,
            balance_threshold=balance_threshold,
            binning_method=binning_method,
            binning_snp_number=binning_snp_number,
            transition_param_unit=transition_param_unit,
            transition_min_segment_bins=transition_min_segment_bins,
            transition_u_gamma_prior=transition_u_gamma_prior,
            transition_u_gamma_factor=transition_u_gamma_factor,
            transition_u_rho_prior=transition_u_rho_prior,
            u_emission_other_prior=u_emission_other_prior,
            params_vis_prefix=params_vis_prefix,
            verbose=verbose,
        )
        cache.set_flag(key, "estimate_params", {
            "emission_source": params_dict["emission_source"],
            "transition_source": params_dict["transition_source"],
                "binning_method": binning_method,
                "binning_snp_number": binning_snp_number,
                "transition_param_unit": transition_param_unit,
                "transition_min_segment_bins": transition_min_segment_bins,
                "transition_u_gamma_prior": transition_u_gamma_prior,
                "transition_u_gamma_factor": transition_u_gamma_factor,
                "transition_u_rho_prior": transition_u_rho_prior,
                "u_emission_other_prior": u_emission_other_prior,
            })

        counts_df = params_dict["counts_df"]
        original_counts_df = params_dict["original_counts_df"]
        resolved_emission_params = params_dict["emission_params"]
        transition_param_dict = params_dict["transition_params"]
        alpha_est = transition_param_dict["alpha"]
        beta_est = transition_param_dict["beta"]
        gamma_est = transition_param_dict["gamma"]
        rho_est = transition_param_dict["rho"]

        result["steps_completed"].append("estimate_params")

        if binning_method == "by_snp":
            result["steps_completed"].append("binning")
        if params_dict.get("rough_labels") is not None:
            result["steps_completed"].append("rough_classification")

        params_file = output_path / f"{sample_id}_hmm_params.json"
        if chromosome:
            params_file = output_path / f"{sample_id}_{chromosome}_hmm_params.json"

        import json
        json_payload = {
            "emission_params": {k: list(v) for k, v in resolved_emission_params.items()},
            "emission_source": params_dict["emission_source"],
            "transition_params": transition_param_dict,
            "transition_source": params_dict["transition_source"],
            "transition_param_unit": params_dict["transition_param_unit"],
            "transition_scale_info": params_dict["transition_scale_info"],
            "u_emission_other_prior": params_dict["u_emission_other_prior"],
            "u_emission_prior_applied": params_dict["u_emission_prior_applied"],
            "binning_method": binning_method,
            "binning_snp_number": binning_snp_number,
        }
        with open(params_file, "w") as f:
            json.dump(json_payload, f, indent=2)
        result["output_files"]["hmm_params"] = str(params_file)

        if estimate_only:
            result["steps_completed"].append("estimate_only")
            result["results"] = params_dict
            result["success"] = True
            logger.info(f"Completed parameter estimation for sample {sample_id}")
            return result

        # Step 4: HMM inference
        if not cache.exists(key, "hmm_inference"):
            logger.info("Step 3/4: Running HMM inference")
            hmm = SisterChromatidHMM(
                alpha=alpha_est,
                beta=beta_est,
                gamma=gamma_est,
                rho=rho_est,
                emission_params=resolved_emission_params
            )

            positions = counts_df['position'].values
            results = hmm.decode_with_blocks(
                read_counts=counts_df,
                chromosome=chromosome or "all",
                positions=positions,
                sample_id=sample_id
            )
            if verbose >= 1 and params_vis_prefix is not None:
                suffix = "final_assigned_label_count"
                if chromosome:
                    suffix = f"{suffix}.{chromosome}"
                _plot_assigned_label_count_png(
                    counts_df=counts_df,
                    assigned_labels=results.states,
                    output_path=_diagnostic_path(params_vis_prefix, suffix),
                )

            # Convert unit-level coordinates back to SNP coordinates if binning was used
            if binning_method == "by_snp":
                logger.info(f"  Converting {len(results.states)} units back to {len(original_counts_df)} SNPs")
                results = convert_unit_coordinates_to_snp(
                    results,
                    original_counts_df,
                    n_snps=binning_snp_number
                )

            cache.save(results, key, "hmm_inference")
            cache.set_flag(key, "hmm_inference", {
                "n_snps": len(results.states),
                "n_blocks": len(results.blocks),
            })
        else:
            logger.info("Step 3/4: Using cached HMM inference")
            results = cache.load(key, "hmm_inference")

        result["steps_completed"].append("hmm_inference")
        result["results"] = results

        # Save output files
        blocks_file = output_path / f"{sample_id}_blocks.tsv"
        if chromosome:
            blocks_file = output_path / f"{sample_id}_{chromosome}_blocks.tsv"

        blocks_df = pd.DataFrame([
            {
                'chromosome': b.chromosome,
                'start_pos': b.start_pos,
                'end_pos': b.end_pos,
                'state': b.state,
                'confidence': b.confidence,
                'num_snps': b.num_snps,
            }
            for b in results.blocks
        ])
        blocks_df.to_csv(blocks_file, sep='\t', index=False)
        result["output_files"]["blocks"] = str(blocks_file)

        states_file = output_path / f"{sample_id}_states.tsv"
        if chromosome:
            states_file = output_path / f"{sample_id}_{chromosome}_states.tsv"

        states_df = pd.DataFrame({
            'chromosome': chromosome or "all",
            'position': results.positions,
            'state': results.states,
        })
        states_df.to_csv(states_file, sep='\t', index=False)
        result["output_files"]["states"] = str(states_file)

        stats_file = output_path / f"{sample_id}_statistics.json"
        if chromosome:
            stats_file = output_path / f"{sample_id}_{chromosome}_statistics.json"

        import json
        with open(stats_file, 'w') as f:
            json.dump(results.statistics, f, indent=2)
        result["output_files"]["statistics"] = str(stats_file)

        result["success"] = True
        logger.info(f"Completed processing sample {sample_id}")

    except Exception as e:
        logger.error(f"Error processing sample {sample_id}: {e}")
        result["error"] = str(e)

    return result


def detect_sister_chromatid_states_end2end(
    bam_file: str,
    snp_file: str,
    output_dir: str,
    sample_id: Optional[str] = None,
    chromosomes: Optional[List[str]] = None,
    cache_dir: Optional[str] = None,
    transition_params: Union[str, Dict[str, float]] = "default",  # "default" | "estimated" | dict
    emission_params: Optional[Dict[str, Tuple[float, float, float]]] = None,
    alpha: float = 0.01,
    beta: float = 0.02,
    gamma: float = 0.001,
    rho: float = 0.1,
    depth_threshold: int = 2,
    other_threshold: float = 0.4,
    balance_threshold: float = 0.25,
    binning_method: Optional[str] = None,
    binning_snp_number: int = 250,
    transition_param_unit: str = "step",
    estimate_use: Union[str, List[str]] = "all",
    transition_min_segment_bins: int = 3,
    transition_u_gamma_prior: Optional[float] = None,
    transition_u_gamma_factor: float = 0.1,
    transition_u_rho_prior: float = 0.32,
    u_emission_other_prior: float = 0.05,
    estimate_only: bool = False,
    params_vis_prefix: Optional[Union[str, Path]] = None,
    verbose: int = 0,
) -> Union[SisterChromatidResults, Dict[str, SisterChromatidResults], Dict[str, Any]]:
    """
    End-to-end sister chromatid state inference from BAM file.

    This is the main entry point for the mII crossover detection pipeline.
    It processes BAM files through all steps:
    1. Extract read counts at SNP positions
    2. Binning (optional): merge adjacent SNPs into units
    3. Rough classification for parameter initialization
    4. Estimate emission and transition parameters
    5. HMM inference with Viterbi decoding
    6. Group states into contiguous blocks
    7. Convert unit-level coordinates back to SNP coordinates (if binning was used)

    Args:
        bam_file: Path to input BAM file
        snp_file: Path to SNP file (VCF or TSV)
        output_dir: Output directory for results
        sample_id: Sample ID (default: derived from bam_file name)
        chromosomes: Optional list of chromosomes to analyze (default: all)
        cache_dir: Cache directory (default: output_dir/.cache)
        transition_params: Transition parameter source ("default" | "estimated") or alpha/beta/gamma/rho dict
        emission_params: Optional specified emission probabilities by PP/PM/MM/U label
        alpha: PP<->PM transition probability (default: 0.01, used when transition_params="default")
        beta: PM->PP/MM transition probability (default: 0.02, used when transition_params="default")
        gamma: probability of entering U state (default: 0.001, used when transition_params="default")
        rho: probability of leaving U state (default: 0.1, used when transition_params="default")
        depth_threshold: Minimum depth for confident call (default: 2)
        other_threshold: Maximum other fraction for confident call (default: 0.4)
        balance_threshold: Threshold for classifying as balanced (default: 0.25)
                          PP: ref_frac > 0.75, MM: alt_frac > 0.75
        binning_method: Binning method (None for single-SNP mode, "by_snp" for merging by SNP count)
        binning_snp_number: Number of SNPs per unit (default: 250)
        transition_param_unit: Unit for specified transition parameters:
            "step", "mb", or "snp". Estimated parameters are always per-step.
        estimate_use: Chromosomes used for parameter estimation when chromosomes
            is provided. Use "all" for all extracted chromosomes, or a list of
            chromosome names.
        transition_min_segment_bins: Raw rough-label segments shorter than this
            many bins are marked as gaps and skipped for transition estimation.
        transition_u_gamma_prior: Optional fixed gamma prior used for U entry
            when binning_method="by_snp" and transition_params="estimated". If
            None, uses transition_u_gamma_factor * min(alpha, beta).
        transition_u_gamma_factor: Factor for deriving the by-SNP U-entry prior
            from estimated alpha/beta when transition_u_gamma_prior is None.
        transition_u_rho_prior: By-SNP prior for U -> each biological state.
        u_emission_other_prior: By-SNP fixed U-emission other probability
            used when emission_params is not specified.
        estimate_only: If True, stop after parameter estimation and return parameter diagnostics
        params_vis_prefix: Prefix for parameter diagnostic PNG files when verbose >= 1.
            Writes rough_label_classification, rough_decision_distributions,
            rough_segment_lengths, rough_label_filtering, emission_matrix, and
            transition_matrix PNGs. During full decoding it also writes
            final_assigned_label_count PNGs.
        verbose: If >=1, print parameter estimation diagnostics

    Returns:
        SisterChromatidResults object (single chromosome) or dictionary (multiple chromosomes)

    Example:
        >>> # Default parameters (fixed transition matrix)
        >>> results = detect_sister_chromatid_states_end2end(
        ...     bam_file="mII_cell.bam",
        ...     snp_file="snps.vcf",
        ...     output_dir="output/"
        ... )

        >>> # Estimated transition parameters from data
        >>> results = detect_sister_chromatid_states_end2end(
        ...     bam_file="mII_cell.bam",
        ...     snp_file="snps.vcf",
        ...     output_dir="output/",
        ...     transition_params="estimated"
        ... )

        >>> # Custom binning (250 SNPs per unit)
        >>> results = detect_sister_chromatid_states_end2end(
        ...     bam_file="mII_cell.bam",
        ...     snp_file="snps.vcf",
        ...     output_dir="output/",
        ...     binning_method="by_snp",
        ...     binning_snp_number=250
        ... )
    """
    if sample_id is None:
        sample_id = Path(bam_file).stem

    if not Path(bam_file).exists():
        raise FileNotFoundError(f"BAM file not found: {bam_file}")
    if not Path(snp_file).exists():
        raise FileNotFoundError(f"SNP file not found: {snp_file}")

    if cache_dir is None:
        cache_dir = str(Path(output_dir) / ".cache")

    logger.info(f"Starting sister chromatid inference pipeline")
    logger.info(f"  Sample: {sample_id}")
    logger.info(f"  BAM file: {bam_file}")
    logger.info(f"  SNPs: {snp_file}")
    logger.info(f"  Output: {output_dir}")

    if chromosomes is None:
        # Process all chromosomes together
        result = _process_sample_mii(
            sample_id=sample_id,
            bam_file=bam_file,
            snp_file=snp_file,
            output_dir=output_dir,
            cache_dir=cache_dir,
            chromosome=None,
            transition_params=transition_params,
            emission_params=emission_params,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            rho=rho,
            depth_threshold=depth_threshold,
            other_threshold=other_threshold,
            balance_threshold=balance_threshold,
            binning_method=binning_method,
            binning_snp_number=binning_snp_number,
            transition_param_unit=transition_param_unit,
            transition_min_segment_bins=transition_min_segment_bins,
            transition_u_gamma_prior=transition_u_gamma_prior,
            transition_u_gamma_factor=transition_u_gamma_factor,
            transition_u_rho_prior=transition_u_rho_prior,
            u_emission_other_prior=u_emission_other_prior,
            estimate_only=estimate_only,
            params_vis_prefix=params_vis_prefix,
            verbose=verbose,
        )

        if not result["success"]:
            raise RuntimeError(f"Pipeline failed: {result.get('error', 'Unknown error')}")

        return result["results"]

    else:
        # Estimate parameters once from the requested genome-wide subset, then
        # decode each chromosome independently to avoid cross-chromosome HMM
        # transitions.
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        counts_file = output_path / f"{sample_id}_read_counts.tsv"
        counts_df_all = extract_read_counts(bam_file, snp_file, str(counts_file))
        if len(counts_df_all) == 0:
            raise ValueError("No SNPs found")

        if estimate_use == "all":
            estimate_counts_df = counts_df_all
        elif isinstance(estimate_use, list):
            estimate_counts_df = counts_df_all[counts_df_all["chromosome"].isin(estimate_use)]
            if len(estimate_counts_df) == 0:
                raise ValueError(f"No SNPs found for estimate_use={estimate_use}")
        else:
            raise ValueError("estimate_use must be 'all' or a list of chromosome names")

        params_dict = estimate_mii_hmm_parameters(
            estimate_counts_df,
            emission_params=emission_params,
            transition_params=transition_params,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            rho=rho,
            depth_threshold=depth_threshold,
            other_threshold=other_threshold,
            balance_threshold=balance_threshold,
            binning_method=binning_method,
            binning_snp_number=binning_snp_number,
            transition_param_unit=transition_param_unit,
            transition_min_segment_bins=transition_min_segment_bins,
            transition_u_gamma_prior=transition_u_gamma_prior,
            transition_u_gamma_factor=transition_u_gamma_factor,
            transition_u_rho_prior=transition_u_rho_prior,
            u_emission_other_prior=u_emission_other_prior,
            params_vis_prefix=params_vis_prefix,
            verbose=verbose,
        )

        params_file = output_path / f"{sample_id}_hmm_params.json"
        import json
        with open(params_file, "w") as f:
            json.dump({
                "emission_params": {k: list(v) for k, v in params_dict["emission_params"].items()},
                "emission_source": params_dict["emission_source"],
                "transition_params": params_dict["transition_params"],
                "transition_source": params_dict["transition_source"],
                "transition_param_unit": params_dict["transition_param_unit"],
                "transition_scale_info": params_dict["transition_scale_info"],
                "estimate_use": estimate_use,
                "binning_method": binning_method,
                "binning_snp_number": binning_snp_number,
                "transition_min_segment_bins": transition_min_segment_bins,
                "transition_u_gamma_prior": transition_u_gamma_prior,
                "transition_u_gamma_factor": transition_u_gamma_factor,
                "transition_u_rho_prior": transition_u_rho_prior,
                "u_emission_other_prior": params_dict["u_emission_other_prior"],
                "u_emission_prior_applied": params_dict["u_emission_prior_applied"],
            }, f, indent=2)

        if estimate_only:
            return params_dict

        all_results = {}
        transition_param_dict = params_dict["transition_params"]
        resolved_emission_params = params_dict["emission_params"]

        for chrom in chromosomes:
            logger.info(f"Processing chromosome {chrom}")
            chrom_counts = counts_df_all[counts_df_all["chromosome"] == chrom].copy()
            if len(chrom_counts) == 0:
                logger.error(f"  ✗ {chrom}: no SNPs found")
                continue

            counts_df, original_counts_df = _prepare_counts_for_estimation(
                chrom_counts,
                binning_method=binning_method,
                binning_snp_number=binning_snp_number,
            )
            hmm = SisterChromatidHMM(
                alpha=transition_param_dict["alpha"],
                beta=transition_param_dict["beta"],
                gamma=transition_param_dict["gamma"],
                rho=transition_param_dict["rho"],
                emission_params=resolved_emission_params,
            )
            result_obj = hmm.decode_with_blocks(
                read_counts=counts_df,
                chromosome=chrom,
                positions=counts_df["position"].values,
                sample_id=sample_id,
            )
            if verbose >= 1 and params_vis_prefix is not None:
                _plot_assigned_label_count_png(
                    counts_df=counts_df,
                    assigned_labels=result_obj.states,
                    output_path=_diagnostic_path(params_vis_prefix, f"final_assigned_label_count.{chrom}"),
                )
            if binning_method == "by_snp":
                result_obj = convert_unit_coordinates_to_snp(
                    result_obj,
                    original_counts_df,
                    n_snps=binning_snp_number,
                )

            blocks_file = output_path / f"{sample_id}_{chrom}_blocks.tsv"
            blocks_df = pd.DataFrame([
                {
                    "chromosome": b.chromosome,
                    "start_pos": b.start_pos,
                    "end_pos": b.end_pos,
                    "state": b.state,
                    "confidence": b.confidence,
                    "num_snps": b.num_snps,
                }
                for b in result_obj.blocks
            ])
            blocks_df.to_csv(blocks_file, sep="\t", index=False)

            states_file = output_path / f"{sample_id}_{chrom}_states.tsv"
            states_df = pd.DataFrame({
                "chromosome": chrom,
                "position": result_obj.positions,
                "state": result_obj.states,
            })
            states_df.to_csv(states_file, sep="\t", index=False)

            stats_file = output_path / f"{sample_id}_{chrom}_statistics.json"
            with open(stats_file, "w") as f:
                json.dump(result_obj.statistics, f, indent=2)

            all_results[chrom] = result_obj
            logger.info(f"  ✓ {chrom}: {len(result_obj.blocks)} blocks")

        return all_results


def save_results_mii(
    results: Union[SisterChromatidResults, Dict[str, SisterChromatidResults]],
    output_dir: str,
    sample_id: str
) -> Dict[str, str]:
    """
    Save sister chromatid inference results to files.

    Args:
        results: Results object or dictionary
        output_dir: Output directory
        sample_id: Sample identifier

    Returns:
        Dictionary mapping output type to file path
    """
    from pathlib import Path

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    output_files = {}

    if isinstance(results, dict):
        for chrom, result in results.items():
            chrom_output = save_results_mii(result, output_dir, f"{sample_id}_{chrom}")
            output_files.update(chrom_output)
    else:
        # Save blocks
        blocks_file = output_path / f"{sample_id}_blocks.tsv"
        blocks_df = pd.DataFrame([
            {
                'chromosome': b.chromosome,
                'start_pos': b.start_pos,
                'end_pos': b.end_pos,
                'state': b.state,
                'confidence': b.confidence,
                'num_snps': b.num_snps,
            }
            for b in results.blocks
        ])
        blocks_df.to_csv(blocks_file, sep='\t', index=False)
        output_files['blocks'] = str(blocks_file)

        # Save states
        states_file = output_path / f"{sample_id}_states.tsv"
        states_df = pd.DataFrame({
            'chromosome': results.chromosome,
            'position': results.positions,
            'state': results.states,
        })
        states_df.to_csv(states_file, sep='\t', index=False)
        output_files['states'] = str(states_file)

        # Save statistics
        import json
        stats_file = output_path / f"{sample_id}_statistics.json"
        with open(stats_file, 'w') as f:
            json.dump(results.statistics, f, indent=2)
        output_files['statistics'] = str(stats_file)

    return output_files
