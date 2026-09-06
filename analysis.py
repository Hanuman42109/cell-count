"""
analysis.py

Reusable query/analysis functions for Parts 2-4 of the project.
Both run_pipeline.py (static output generation) and dashboard.py
(interactive Streamlit app) import from this module so the logic
lives in exactly one place.
"""

import sqlite3
from pathlib import Path

import pandas as pd
from scipy import stats

DB_PATH = Path(__file__).parent / "cell_counts.db"


def get_connection(db_path: Path = DB_PATH) -> sqlite3.Connection:
    return sqlite3.connect(db_path)


# ---------------------------------------------------------------------------
# Part 2: relative frequency table
# ---------------------------------------------------------------------------

FREQUENCY_QUERY = """
SELECT
    s.sample_id AS sample,
    SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS total_count,
    p.name AS population,
    cc.count AS count,
    100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS percentage
FROM cell_counts cc
JOIN populations p ON cc.population_id = p.population_id
JOIN samples s ON cc.sample_id = s.sample_id
ORDER BY s.sample_id, p.name;
"""


def get_frequency_table(conn: sqlite3.Connection) -> pd.DataFrame:
    """Part 2: one row per (sample, population) with count and % of sample total."""
    df = pd.read_sql_query(FREQUENCY_QUERY, conn)
    return df


# ---------------------------------------------------------------------------
# Part 3: responders vs non-responders (melanoma, miraclib, PBMC)
# ---------------------------------------------------------------------------

RESPONSE_COMPARISON_QUERY = """
SELECT
    s.sample_id AS sample,
    sub.subject_id AS subject,
    sub.response AS response,
    p.name AS population,
    100.0 * cc.count / tot.total_count AS percentage
FROM cell_counts cc
JOIN populations p ON cc.population_id = p.population_id
JOIN samples s ON cc.sample_id = s.sample_id
JOIN subjects sub ON s.subject_id = sub.subject_id
JOIN (
    SELECT sample_id, SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
) tot ON tot.sample_id = cc.sample_id
WHERE sub.condition = 'melanoma'
  AND sub.treatment = 'miraclib'
  AND s.sample_type = 'PBMC'
  AND sub.response IS NOT NULL
ORDER BY population, response, sample;
"""


def get_responder_comparison_data(conn: sqlite3.Connection) -> pd.DataFrame:
    """Part 3 input data: melanoma + miraclib + PBMC samples with response label."""
    return pd.read_sql_query(RESPONSE_COMPARISON_QUERY, conn)


def compute_response_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Mann-Whitney U test per population comparing responder vs non-responder
    relative frequencies. Returns one row per population with the test
    statistic, p-value, group medians, and a rank-biserial effect size.
    """
    results = []
    for population, sub in df.groupby("population"):
        yes = sub.loc[sub.response == "yes", "percentage"]
        no = sub.loc[sub.response == "no", "percentage"]

        if len(yes) < 2 or len(no) < 2:
            continue

        u_stat, p_value = stats.mannwhitneyu(yes, no, alternative="two-sided")
        # rank-biserial effect size derived from U
        n1, n2 = len(yes), len(no)
        effect_size = 1 - (2 * u_stat) / (n1 * n2)

        results.append(
            {
                "population": population,
                "n_responders": n1,
                "n_non_responders": n2,
                "median_responders_pct": yes.median(),
                "median_non_responders_pct": no.median(),
                "u_statistic": u_stat,
                "p_value": p_value,
                "effect_size_rank_biserial": effect_size,
                "significant_p<0.05": p_value < 0.05,
            }
        )

    result_df = pd.DataFrame(results).sort_values("p_value")
    return result_df


# ---------------------------------------------------------------------------
# Part 4: baseline melanoma / miraclib / PBMC subset
# ---------------------------------------------------------------------------

BASELINE_SUBSET_QUERY = """
SELECT
    s.sample_id AS sample,
    sub.subject_id AS subject,
    p.name AS project,
    sub.response AS response,
    sub.sex AS sex
FROM samples s
JOIN subjects sub ON s.subject_id = sub.subject_id
JOIN projects p ON sub.project_id = p.project_id
WHERE sub.condition = 'melanoma'
  AND sub.treatment = 'miraclib'
  AND s.sample_type = 'PBMC'
  AND s.time_from_treatment_start = 0;
"""


def get_baseline_subset(conn: sqlite3.Connection) -> pd.DataFrame:
    """Part 4.1: melanoma + PBMC + miraclib samples at baseline (time = 0)."""
    return pd.read_sql_query(BASELINE_SUBSET_QUERY, conn)


def summarize_baseline_subset(df: pd.DataFrame) -> dict:
    """Part 4.2: breakdowns by project, response, and sex (subject-level counts)."""
    by_project = df.groupby("project")["sample"].nunique().rename("n_samples")

    subj_level = df.drop_duplicates("subject")
    by_response = subj_level["response"].value_counts().rename("n_subjects")
    by_sex = subj_level["sex"].value_counts().rename("n_subjects")

    return {
        "samples_per_project": by_project.reset_index(),
        "subjects_by_response": by_response.reset_index().rename(
            columns={"index": "response"}
        ),
        "subjects_by_sex": by_sex.reset_index().rename(columns={"index": "sex"}),
    }


MALE_MELANOMA_BASELINE_BCELL_QUERY = """
SELECT cc.count AS b_cell_count
FROM cell_counts cc
JOIN populations p ON cc.population_id = p.population_id
JOIN samples s ON cc.sample_id = s.sample_id
JOIN subjects sub ON s.subject_id = sub.subject_id
WHERE sub.condition = 'melanoma'
  AND sub.sex = 'M'
  AND sub.response = 'yes'
  AND s.time_from_treatment_start = 0
  AND p.name = 'b_cell';
"""


def get_avg_bcell_male_melanoma_responders_baseline(conn: sqlite3.Connection) -> float:
    """
    Part 4.3: Considering melanoma males of ALL sample and treatment types,
    the average B cell count for responders at time = 0.
    """
    df = pd.read_sql_query(MALE_MELANOMA_BASELINE_BCELL_QUERY, conn)
    return round(df["b_cell_count"].mean(), 2)


if __name__ == "__main__":
    conn = get_connection()
    freq = get_frequency_table(conn)
    print("Part 2 sample rows:")
    print(freq.head())

    comp = get_responder_comparison_data(conn)
    stats_df = compute_response_statistics(comp)
    print("\nPart 3 statistics:")
    print(stats_df)

    baseline = get_baseline_subset(conn)
    summary = summarize_baseline_subset(baseline)
    print("\nPart 4 breakdowns:")
    for k, v in summary.items():
        print(k)
        print(v)

    avg_bcell = get_avg_bcell_male_melanoma_responders_baseline(conn)
    print(f"\nAvg B cell count, melanoma males, responders, t=0: {avg_bcell:.2f}")

    conn.close()
