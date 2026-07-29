"""Figures for the SIGMAR1 / amyloid analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .blast import SIGNIFICANCE_THRESHOLD, Hit  # noqa: E402
from .interpro import Signature  # noqa: E402
from .stringdb import Interaction  # noqa: E402

BLUE = "#1f77b4"
ORANGE = "#e8a33d"
RED = "#b4413c"
GREY = "#b8c4cc"


def plot_blast_null_model(
    hits: list[Hit], null: dict[str, float], path: Path
) -> Path:
    """Observed best-hit bitscores against the shuffled-sequence null."""
    top = hits[:12]
    labels = [hit.subject_name.replace("_HUMAN", "") for hit in top]
    scores = [hit.bitscore for hit in top]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(labels, scores, color=GREY, edgecolor="white")

    ax.axhline(
        null["percentile_95"],
        color=RED,
        linestyle="--",
        linewidth=1.4,
        label=f"95th percentile of shuffled null ({null['percentile_95']:.1f} bits)",
    )
    ax.axhline(
        null["max_best_bitscore"],
        color=RED,
        linestyle=":",
        linewidth=1.4,
        label=f"maximum of shuffled null ({null['max_best_bitscore']:.1f} bits)",
    )

    # Leave headroom above the null lines so the legend never sits on top of them.
    ceiling = max([*scores, null["max_best_bitscore"]])
    ax.set_ylim(0, ceiling * 1.45)

    ax.set_ylabel("BLAST bitscore")
    ax.set_title(
        "No SIGMAR1 hit exceeds what shuffled sequences achieve by chance",
        fontsize=11,
    )
    ax.tick_params(axis="x", rotation=60)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def plot_evalue_distribution(hits: list[Hit], path: Path) -> Path:
    """Where every hit falls relative to the significance threshold."""
    fig, ax = plt.subplots(figsize=(7, 4))

    evalues = [hit.evalue for hit in hits]
    ax.scatter(
        range(1, len(evalues) + 1),
        evalues,
        s=40,
        color=BLUE,
        zorder=3,
    )
    ax.axhline(
        SIGNIFICANCE_THRESHOLD,
        color=RED,
        linestyle="--",
        linewidth=1.4,
        label=f"significance threshold (E = {SIGNIFICANCE_THRESHOLD})",
    )

    ax.set_yscale("log")
    ax.set_xlabel("Hit rank")
    ax.set_ylabel("E-value (log scale)")
    ax.set_title("Every hit falls above the significance threshold", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#eef1f3", zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


TOPOLOGY_COLOURS = {"cytoplasmic": "#c3d5e4", "non-cytoplasmic": "#e6dcc6"}


def plot_domain_architecture(
    signatures: list[Signature], protein_length: int, path: Path
) -> Path:
    """Draw predicted topology plus one track per family-assigning analysis.

    Each analysis gets its own row, so overlapping matches from PANTHER and Pfam
    stay legible instead of printing on top of one another.
    """
    from .interpro import summarise

    summary = summarise(signatures)
    families = [s for s in signatures if s.interpro_accession not in {"-", ""}]
    analyses = sorted({s.analysis for s in families})

    row_height = 0.5
    fig, ax = plt.subplots(figsize=(9, 2.2 + 0.55 * len(analyses)))

    backbone_y = len(analyses) * row_height + 0.35

    # Topology bands along the backbone.
    for segment in summary["topology"]:
        ax.add_patch(
            plt.Rectangle(
                (segment["start"], backbone_y),
                segment["stop"] - segment["start"] + 1,
                0.3,
                color=TOPOLOGY_COLOURS.get(segment["side"], "#e9edf0"),
                zorder=1,
            )
        )

    # Merged transmembrane segments on top of the backbone.
    for segment in summary["transmembrane_segments"]:
        width = segment["stop"] - segment["start"] + 1
        ax.add_patch(
            plt.Rectangle(
                (segment["start"], backbone_y - 0.06),
                width,
                0.42,
                color=ORANGE,
                zorder=3,
            )
        )
        ax.text(
            segment["start"] + width / 2,
            backbone_y + 0.42,
            f"TM {segment['start']}–{segment['stop']}",
            ha="center",
            fontsize=8,
            color="#7a5518",
        )

    # One row per analysis.
    for index, analysis in enumerate(analyses):
        y = index * row_height
        for signature in (s for s in families if s.analysis == analysis):
            ax.add_patch(
                plt.Rectangle(
                    (signature.start, y),
                    signature.length,
                    0.3,
                    color=BLUE,
                    alpha=0.55,
                    zorder=2,
                )
            )
            ax.text(
                signature.start + signature.length / 2,
                y + 0.15,
                f"{signature.interpro_accession}  {signature.start}–{signature.stop}",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                zorder=4,
            )
        ax.text(-4, y + 0.15, analysis, ha="right", va="center", fontsize=8.5)

    ax.set_xlim(-protein_length * 0.16, protein_length + 6)
    ax.set_ylim(-0.3, backbone_y + 0.85)
    ax.set_yticks([])
    ax.set_xlabel("Residue")
    ax.set_title(
        f"SIGMAR1 ({protein_length} aa): {summary['interpro_families'][0]['description']}",
        fontsize=11,
    )

    for side, colour in TOPOLOGY_COLOURS.items():
        ax.add_patch(plt.Rectangle((0, 0), 0, 0, color=colour, label=side))
    ax.add_patch(plt.Rectangle((0, 0), 0, 0, color=ORANGE, label="transmembrane"))
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="lower right")

    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.tight_layout()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def plot_string_partners(
    interactions: list[Interaction], path: Path, top: int = 20
) -> Path:
    """Rank STRING partners by confidence, flagging the evidence class."""
    selected = interactions[:top][::-1]
    labels = [i.partner_b for i in selected]
    scores = [i.combined_score for i in selected]
    colours = [BLUE if i.has_physical_evidence else GREY for i in selected]

    fig, ax = plt.subplots(figsize=(7, max(4, len(selected) * 0.3)))
    ax.barh(labels, scores, color=colours)

    ax.axvline(0.9, color=RED, linestyle="--", linewidth=1.2)
    ax.text(0.898, -0.7, "highest confidence", ha="right", fontsize=8, color=RED)

    ax.set_xlim(0, 1.02)
    ax.set_xlabel("STRING combined score")
    ax.set_title("SIGMAR1 functional interaction partners", fontsize=11)
    ax.barh([], [], color=BLUE, label="experimental / curated evidence")
    ax.barh([], [], color=GREY, label="predicted or text-mined only")
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path
