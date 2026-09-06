#!/usr/bin/env python3
"""
load_data.py

Initializes a SQLite database (cell_counts.db) with a normalized relational
schema and loads all rows from cell-count.csv into it.

Usage:
    python load_data.py

No command-line arguments required. Running this script will:
  1. Create (or recreate) cell_counts.db in the repository root.
  2. Create the schema (projects, subjects, samples, populations, cell_counts).
  3. Parse cell-count.csv and populate all tables.

Schema design
-------------
The raw CSV is "wide": one row per sample, with one column per cell
population (b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte) plus subject-
and sample-level metadata all mixed together. That flat structure repeats
subject metadata (project, condition, age, sex, treatment, response) on
every sample row, and hard-codes the five cell populations as columns,
which makes it awkward to add a sixth population or to run per-population
SQL analytics later.

Instead we normalize into five tables:

  projects(project_id PK, name)
      One row per project (prj1, prj2, prj3, ...).

  subjects(subject_id PK, project_id FK, condition, age, sex, treatment,
           response)
      One row per patient/subject. These attributes are constant across
      all of a subject's samples in the source data, so they belong on
      the subject, not repeated per sample.

  samples(sample_id PK, subject_id FK, sample_type, time_from_treatment_start)
      One row per biological sample (a subject contributes multiple
      samples over time / across sample types).

  populations(population_id PK, name)
      One row per immune cell population (b_cell, cd8_t_cell, ...). Kept
      as a lookup table rather than fixed columns so new populations can
      be added without altering table schemas or existing queries.

  cell_counts(sample_id FK, population_id FK, count)
      One row per (sample, population) pair -- the long/tidy form of the
      wide CSV. Composite primary key (sample_id, population_id).

This is a standard 3NF "star-ish" design: subjects/samples/populations are
dimension tables, cell_counts is the fact table. See README.md for the
rationale behind this design and how it scales to many more projects,
samples, and analysis types.
"""

import csv
import sqlite3
from pathlib import Path

DB_PATH = Path(__file__).parent / "cell_counts.db"
CSV_PATH = Path(__file__).parent / "cell-count.csv"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
CREATE TABLE projects (
    project_id  INTEGER PRIMARY KEY AUTOINCREMENT,
    name        TEXT NOT NULL UNIQUE
);

CREATE TABLE subjects (
    subject_id  TEXT PRIMARY KEY,
    project_id  INTEGER NOT NULL REFERENCES projects(project_id),
    condition   TEXT NOT NULL,
    age         INTEGER NOT NULL,
    sex         TEXT NOT NULL CHECK (sex IN ('M', 'F')),
    treatment   TEXT NOT NULL,
    response    TEXT CHECK (response IN ('yes', 'no') OR response IS NULL)
);

CREATE TABLE samples (
    sample_id                   TEXT PRIMARY KEY,
    subject_id                  TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type                 TEXT NOT NULL,
    time_from_treatment_start   INTEGER NOT NULL
);

CREATE TABLE populations (
    population_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    name            TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id       TEXT NOT NULL REFERENCES samples(sample_id),
    population_id   INTEGER NOT NULL REFERENCES populations(population_id),
    count           INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id)
);

CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_samples_subject ON samples(subject_id);
CREATE INDEX idx_cellcounts_sample ON cell_counts(sample_id);
CREATE INDEX idx_cellcounts_population ON cell_counts(population_id);
"""


def build_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA)


def load_csv(conn: sqlite3.Connection, csv_path: Path) -> None:
    cur = conn.cursor()

    # populations lookup table
    pop_ids = {}
    for name in POPULATIONS:
        cur.execute("INSERT INTO populations (name) VALUES (?)", (name,))
        pop_ids[name] = cur.lastrowid

    project_ids = {}
    seen_subjects = set()
    seen_samples = set()

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            project_name = row["project"]
            if project_name not in project_ids:
                cur.execute(
                    "INSERT INTO projects (name) VALUES (?)", (project_name,)
                )
                project_ids[project_name] = cur.lastrowid
            project_id = project_ids[project_name]

            subject_id = row["subject"]
            if subject_id not in seen_subjects:
                seen_subjects.add(subject_id)
                response = row["response"] if row["response"] else None
                cur.execute(
                    """INSERT INTO subjects
                       (subject_id, project_id, condition, age, sex,
                        treatment, response)
                       VALUES (?, ?, ?, ?, ?, ?, ?)""",
                    (
                        subject_id,
                        project_id,
                        row["condition"],
                        int(row["age"]),
                        row["sex"],
                        row["treatment"],
                        response,
                    ),
                )

            sample_id = row["sample"]
            if sample_id not in seen_samples:
                seen_samples.add(sample_id)
                cur.execute(
                    """INSERT INTO samples
                       (sample_id, subject_id, sample_type,
                        time_from_treatment_start)
                       VALUES (?, ?, ?, ?)""",
                    (
                        sample_id,
                        subject_id,
                        row["sample_type"],
                        int(row["time_from_treatment_start"]),
                    ),
                )

            for pop_name in POPULATIONS:
                cur.execute(
                    """INSERT INTO cell_counts (sample_id, population_id, count)
                       VALUES (?, ?, ?)""",
                    (sample_id, pop_ids[pop_name], int(row[pop_name])),
                )

    conn.commit()


def main() -> None:
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"Could not find {CSV_PATH}")

    if DB_PATH.exists():
        DB_PATH.unlink()

    conn = sqlite3.connect(DB_PATH)
    try:
        conn.execute("PRAGMA foreign_keys = ON;")
        build_schema(conn)
        load_csv(conn, CSV_PATH)

        n_projects = conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0]
        n_subjects = conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0]
        n_samples = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        n_counts = conn.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0]
        print(f"Created {DB_PATH.name}")
        print(f"  projects:     {n_projects}")
        print(f"  subjects:     {n_subjects}")
        print(f"  samples:      {n_samples}")
        print(f"  cell_counts:  {n_counts}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
