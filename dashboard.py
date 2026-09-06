#!/usr/bin/env python3
"""
dashboard.py

Interactive Streamlit dashboard presenting the Part 2-4 results computed
from cell_counts.db (created by load_data.py).

Run with:
    streamlit run dashboard.py
"""

from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

import load_data
from analysis import (
    get_avg_bcell_male_melanoma_responders_baseline,
    get_baseline_subset,
    get_connection,
    get_frequency_table,
    get_responder_comparison_data,
    compute_response_statistics,
    summarize_baseline_subset,
)

DB_PATH = Path(__file__).parent / "cell_counts.db"

st.set_page_config(
    page_title="Loblaw Bio | Immune Cell Population Dashboard",
    layout="wide",
)


def require_db():
    """
    Build the database on first load if it doesn't exist yet. This makes
    the dashboard self-sufficient on platforms like Streamlit Community
    Cloud, which only run `streamlit run dashboard.py` and don't invoke
    `make pipeline` first. Locally, `make pipeline` still builds it ahead
    of time so this is typically a no-op.
    """
    if not DB_PATH.exists():
        with st.spinner("First run: building cell_counts.db from cell-count.csv..."):
            load_data.main()


def main():
    require_db()
    conn = get_connection()
    try:
        st.title("Immune Cell Population Dashboard")
        st.caption(
            "Cell count analysis for Loblaw Bio's miraclib clinical trial data."
        )

        tab2, tab3, tab4 = st.tabs(
            ["Part 2 - Frequency Overview", "Part 3 - Responder Analysis", "Part 4 - Baseline Subset"]
        )

        # ---------------- Part 2 ----------------
        with tab2:
            st.header("Relative frequency of each cell population per sample")
            freq = get_frequency_table(conn)

            samples = sorted(freq["sample"].unique())
            col1, col2 = st.columns([1, 3])
            with col1:
                selected = st.multiselect(
                    "Filter by sample (leave empty to show all)",
                    samples,
                    default=[],
                )
            display_df = freq[freq["sample"].isin(selected)] if selected else freq
            st.dataframe(display_df, use_container_width=True, height=400)
            st.caption(f"Showing {len(display_df):,} of {len(freq):,} rows.")

            st.download_button(
                "Download full frequency table (CSV)",
                freq.to_csv(index=False),
                file_name="part2_frequency_table.csv",
                mime="text/csv",
            )

        # ---------------- Part 3 ----------------
        with tab3:
            st.header("Responders vs Non-Responders")
            st.caption(
                "Melanoma patients treated with miraclib, PBMC samples only. "
                "Statistical test: Mann-Whitney U (non-parametric, robust to "
                "small samples and skew)."
            )

            comp = get_responder_comparison_data(conn)

            if comp.empty:
                st.warning("No matching samples found for this comparison.")
            else:
                fig = px.box(
                    comp,
                    x="population",
                    y="percentage",
                    color="response",
                    points="outliers",
                    labels={"percentage": "Relative frequency (%)", "population": "Cell population"},
                    title="Relative frequency by population: responders vs non-responders",
                    color_discrete_map={"yes": "#2E86AB", "no": "#E07A5F"},
                )
                fig.update_layout(boxmode="group")
                st.plotly_chart(fig, use_container_width=True)

                stats_df = compute_response_statistics(comp)
                st.subheader("Statistical test results (Mann-Whitney U)")
                st.dataframe(
                    stats_df.style.format(
                        {
                            "median_responders_pct": "{:.2f}",
                            "median_non_responders_pct": "{:.2f}",
                            "u_statistic": "{:.1f}",
                            "p_value": "{:.4f}",
                            "effect_size_rank_biserial": "{:.3f}",
                        }
                    ),
                    use_container_width=True,
                )

                sig = stats_df[stats_df["significant_p<0.05"]]
                if len(sig) > 0:
                    pops = ", ".join(sig["population"].tolist())
                    st.success(
                        f"**Significant difference (p < 0.05):** {pops}. "
                        "These populations differ between responders and "
                        "non-responders and are candidates for predicting "
                        "response to miraclib."
                    )
                else:
                    st.info("No population reached significance at p < 0.05.")

        # ---------------- Part 4 ----------------
        with tab4:
            st.header("Baseline subset: melanoma, PBMC, miraclib, time = 0")

            baseline = get_baseline_subset(conn)
            summary = summarize_baseline_subset(baseline)

            c1, c2, c3 = st.columns(3)
            with c1:
                st.subheader("Samples per project")
                st.dataframe(summary["samples_per_project"], use_container_width=True)
            with c2:
                st.subheader("Subjects by response")
                st.dataframe(summary["subjects_by_response"], use_container_width=True)
            with c3:
                st.subheader("Subjects by sex")
                st.dataframe(summary["subjects_by_sex"], use_container_width=True)

            st.divider()
            st.subheader("Melanoma males, all sample/treatment types, responders at time = 0")
            avg_bcell = get_avg_bcell_male_melanoma_responders_baseline(conn)
            st.metric("Average B cell count", f"{avg_bcell:.2f}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
