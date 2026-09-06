#!/usr/bin/env python3
"""
run_pipeline.py

Runs the full analysis pipeline (Parts 2-4) against cell_counts.db and
writes the results to outputs/ as CSV files plus a boxplot figure, so
results can be reviewed without needing the dashboard running.

Assumes load_data.py has already been run (see Makefile `make pipeline`,
which runs both in order).
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from analysis import (
    get_avg_bcell_male_melanoma_responders_baseline,
    get_baseline_subset,
    get_connection,
    get_frequency_table,
    get_responder_comparison_data,
    compute_response_statistics,
    summarize_baseline_subset,
)

OUTPUT_DIR = Path(__file__).parent / "outputs"


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    conn = get_connection()

    # Part 2
    freq = get_frequency_table(conn)
    freq.to_csv(OUTPUT_DIR / "part2_frequency_table.csv", index=False)
    print(f"[Part 2] wrote {len(freq)} rows -> outputs/part2_frequency_table.csv")

    # Part 3
    comp = get_responder_comparison_data(conn)
    stats_df = compute_response_statistics(comp)
    stats_df.to_csv(OUTPUT_DIR / "part3_response_statistics.csv", index=False)
    print(f"[Part 3] wrote statistics -> outputs/part3_response_statistics.csv")

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        data=comp,
        x="population",
        y="percentage",
        hue="response",
        ax=ax,
    )
    ax.set_title(
        "Relative Frequency by Population: Responders vs Non-Responders\n"
        "(Melanoma, Miraclib, PBMC samples)"
    )
    ax.set_ylabel("Relative frequency (%)")
    ax.set_xlabel("Cell population")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "part3_boxplot.png", dpi=150)
    plt.close(fig)
    print("[Part 3] wrote boxplot -> outputs/part3_boxplot.png")

    # Part 4
    baseline = get_baseline_subset(conn)
    summary = summarize_baseline_subset(baseline)
    summary["samples_per_project"].to_csv(
        OUTPUT_DIR / "part4_samples_per_project.csv", index=False
    )
    summary["subjects_by_response"].to_csv(
        OUTPUT_DIR / "part4_subjects_by_response.csv", index=False
    )
    summary["subjects_by_sex"].to_csv(
        OUTPUT_DIR / "part4_subjects_by_sex.csv", index=False
    )
    print("[Part 4] wrote breakdown CSVs -> outputs/part4_*.csv")

    avg_bcell = get_avg_bcell_male_melanoma_responders_baseline(conn)
    with open(OUTPUT_DIR / "part4_avg_bcell_male_responders_baseline.txt", "w") as f:
        f.write(
            "Average B cell count, melanoma males (all sample/treatment types), "
            f"responders, time=0: {avg_bcell:.2f}\n"
        )
    print(f"[Part 4] avg B cell count (melanoma males, responders, t=0): {avg_bcell:.2f}")

    conn.close()
    print("\nPipeline complete. All outputs written to outputs/")


if __name__ == "__main__":
    main()
