#! /usr/bin/env python3

import argparse
import multiprocessing
import subprocess
from typing import NamedTuple, Optional
import sqlite3
import textwrap
import numpy as np

class RunResult(NamedTuple):
    q: int
    m: int
    n: int
    number_of_possible_columns: int
    number_of_different_columns: int
    count_configurations: int
    count_recursive_calls: int
    cpu_preparation_time: float
    cpu_calculation_time: float
    distance_counts: list[int]


def run_smart_code_counter(q: int, m: int, n: int) -> Optional[RunResult]:
    args = ["./smart_code_counter", str(q), str(m), str(n)]
    result = subprocess.run(args, capture_output=True)
    if result.returncode != 0:
        return None

    lines = result.stdout.decode('ascii').splitlines()

    header_line = lines[0]
    data_lines = lines[1:]

    # Extract header items

    header_fields_start = header_line.find(":") + 1

    header_fields = header_line[header_fields_start:].split()
    hdict = {}
    while len(header_fields) != 0:
        key = header_fields.pop(0)
        value_string = header_fields.pop(0)
        value_type = float if "." in value_string else int
        value = value_type(value_string)
        hdict[key] = value

    number_of_possible_columns = hdict["number_of_possible_columns"]
    number_of_different_columns = hdict["number_of_different_columns"]
    count_configurations = hdict["count_configurations"]
    count_recursive_calls = hdict["count_recursive_calls"]
    cpu_preparation_time = hdict["cpu_preparation_time"]
    cpu_calculation_time = hdict["cpu_calculation_time"]

    distance_counts = [int(data_line.split()[-1]) for data_line in data_lines]

    return RunResult(
        q, m, n,
        number_of_possible_columns, number_of_different_columns,
        count_configurations, count_recursive_calls,
        cpu_preparation_time, cpu_calculation_time,
        distance_counts
    )


def setup_tables(con):

    query = "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='runs'"
    (count, ) = con.execute(query).fetchone()
    if count != 0:
        return

    # Create tables.

    create_runs_table_query = textwrap.dedent("""\
        CREATE TABLE runs (
            --
            q                            INTEGER NOT NULL,
            m                            INTEGER NOT NULL,
            n                            INTEGER NOT NULL,
            --
            number_of_possible_columns   INTEGER NOT NULL,
            number_of_different_columns  INTEGER NOT NULL,
            count_configurations         INTEGER NOT NULL,
            count_recursive_calls        INTEGER NOT NULL,
            cpu_preparation_time         REAL    NOT NULL,
            cpu_calculation_time         REAL    NOT NULL,
            --
            PRIMARY KEY (q, m, n)
        );
    """)

    con.execute(create_runs_table_query)

    create_counts_table_query = textwrap.dedent("""\
        CREATE TABLE counts (
            --
            q      INTEGER NOT NULL,
            m      INTEGER NOT NULL,
            n      INTEGER NOT NULL,
            d      INTEGER NOT NULL,
            --
            count  TEXT    NOT NULL, -- Can get very big, so we store count values as decimal strings.
            --
            PRIMARY KEY (q, m, n, d)
        );
    """)

    con.execute(create_counts_table_query)


class RunJudge:

    def __init__(self, con, max_n: int, max_runtime: float):
        self.con = con
        self.max_n = max_n
        self.max_runtime = max_runtime

    def __call__(self, q: int, m: int, n: int) -> bool:

        if n > self.max_n:
            return False

        predicted_runtime = self.predict_runtime(q, m, n)
        if predicted_runtime > self.max_runtime:
            print(f"Predicted runtime for q={q}, m={m}, n={n} is {predicted_runtime} seconds; stopping now.")
            return False

        return True

    def predict_runtime(self, q: int, m: int, n: int) -> float:

        max_entries = 10

        query = "SELECT n, total_time FROM (SELECT n, cpu_preparation_time + cpu_calculation_time FROM runs WHERE q = ? AND m = ? ORDER BY n DEC LIMIT ?) ORDER BY n;"
        results = self.con.execute(query, (q, m, max_entries)).fetchall()

        n_values = np.array([f1 for (f1, f2) in results])
        t_values = np.array([f2 for (f1, f2) in results])

        if len(n_values) == 0:
            return 0.0 # Cannot give a good estimate.

        max_time_seen = np.max(t_values)

        if len(n_values) == 1:
            return max_time_seen

        # If we have at least 2 points, we can make a reasonable fit.

        # For two points, we make a linear fit.
        # For three points o more, we make a quadratic fit.
        fit_degree = 2 if len(n_values) <= 2 else 3

        # Note that we fit the log of the times.
        fit = np.polynomial.Polynomial.fit(n_values, np.log(t_values), deg=fit_degree)

        log_t_predicted = fit(n)

        t_predicted = np.exp(log_t_predicted)

        return max(max_time_seen, t_predicted)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--pool-size", default=None, type=int)
    parser.add_argument("--max-q", default=10, type=int)
    parser.add_argument("--max-m", default=6, type=int)
    parser.add_argument("--max-n", default=1000, type=int)
    parser.add_argument("--max-runtime", default=3600.0, type=float)

    args = parser.parse_args()

    db_filename = "code_counter_data.sqlite3"

    con = sqlite3.connect(db_filename)
    setup_tables(con)

    pool = multiprocessing.Pool(args.pool_size)
    queue = multiprocessing.Queue()

    num_runs_active = 0

    run_judge = RunJudge(con, args.max_n, args.max_runtime)

    for q in range(2, args.max_q + 1):
        for m in range(2, args.max_m + 1):
            query = "SELECT MAX(n) FROM runs WHERE q=? AND m=?;"
            (last_n, ) = con.execute(query, (q, m)).fetchone()
            next_n = 1 if last_n is None else last_n + 1
            if run_judge(q, m, next_n):
                pool.apply_async(run_smart_code_counter, (q, m, next_n), callback=queue.put)
                num_runs_active += 1

    # Process results as they come in.

    while num_runs_active != 0:

        run_result = queue.get()
        num_runs_active -= 1

        print(f"[{num_runs_active}] Processing results of job q={run_result.q} m={run_result.m} n={run_result.n} (prep: {run_result.cpu_preparation_time:.6f} ; calc: {run_result.cpu_calculation_time:.6f})")


        query = "INSERT INTO runs(q, m, n, number_of_possible_columns, number_of_different_columns, count_configurations, count_recursive_calls, cpu_preparation_time, cpu_calculation_time) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);"
        con.execute(query, (run_result.q, run_result.m, run_result.n,
                            run_result.number_of_possible_columns, run_result.number_of_different_columns,
                            run_result.count_configurations, run_result.count_recursive_calls,
                            run_result.cpu_preparation_time, run_result.cpu_calculation_time))

        query = "INSERT INTO counts(q, m, n, d, count) VALUES (?, ?, ?, ?, ?);"
        for (d, count) in enumerate(run_result.distance_counts):
            con.execute(query, (run_result.q, run_result.m, run_result.n, d, str(count)))
        con.commit()

        # Determine if we want to initiate a run with the next value of n.
        if run_judge(run_result.q, run_result.m, run_result.n + 1):
            pool.apply_async(run_smart_code_counter, (run_result.q, run_result.m, run_result.n + 1), callback=queue.put)
            num_runs_active += 1

    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
