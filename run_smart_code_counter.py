#! /usr/bin/env python3

import argparse
import multiprocessing
import subprocess
import time
from typing import NamedTuple, Optional
import sqlite3

class RunResult(NamedTuple):
    q: int
    m: int
    n: int
    wallclock_time: float
    nr_of_poss_columns: int
    nr_of_diff_columns: int
    count_of_configurations: int
    count_of_recursive_calls: int
    distance_counts: list[int]


def run_smart_code_counter(q: int, m: int, n: int) -> Optional[RunResult]:
    args = ["./smart-code-counter", str(q), str(m), str(n)]
    t1 = time.monotonic()
    result = subprocess.run(args, capture_output=True)
    t2 = time.monotonic()
    if result.returncode != 0:
        return None

    wallclock_time = (t2 - t1)

    lines = result.stdout.decode('ascii').splitlines()

    distance_counts = [int(lines[i].split()[-1]) for i in range(0, n + 1)]

    last_line_items = lines[-1].split()
    nr_of_poss_columns = int(last_line_items[11])
    nr_of_diff_columns = int(last_line_items[13])
    count_of_configurations = int(last_line_items[15])
    count_of_recursive_calls = int(last_line_items[17])

    return RunResult(q, m, n, wallclock_time, nr_of_poss_columns, nr_of_diff_columns, count_of_configurations, count_of_recursive_calls, distance_counts)


def setup_tables(con):

    query = "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='runs'"
    (count, ) = con.execute(query).fetchone()
    if count != 0:
        return

    # Create tables.

    query = "CREATE TABLE runs (q INTEGER NOT NULL, m INTEGER NOT NULL, n INTEGER NOT NULL, wallclock_time FLOAT NOT NULL, nr_of_poss_columns INTEGER NOT NULL, nr_of_diff_columns INTEGER NOT NULL, count_of_configurations INTEGER NOT NULL, count_of_recursive_calls INTEGER NOT NULL, PRIMARY KEY (q, m, n));"
    con.execute(query)

    query = "CREATE TABLE counts (q INTEGER NOT NULL, m INTEGER NOT NULL, n INTEGER NOT NULL, d INTEGER NOT NULL, count TEXT, PRIMARY KEY (q, m, n, d));"
    con.execute(query)


def main():

    db_filename = "code_counter.sqlite3"

    con = sqlite3.connect(db_filename)
    setup_tables(con)

    pool = multiprocessing.Pool()
    queue = multiprocessing.Queue()

    num_jobs = 0
    for q in (2, 3, 4, 5, 6):
        for m in (2, 3, 4):
            for n in range(1, 21):
                query = "SELECT COUNT(*) FROM runs WHERE q=? AND m=? and n=?;"
                (count, ) = con.execute(query, (q, m, n)).fetchone()
                if count != 0:
                    continue

                pool.apply_async(run_smart_code_counter, (q, m, n), callback=queue.put)
                num_jobs += 1

    # Process results as they come in.

    print(f"Number of jobs submitted: {num_jobs}")

    while num_jobs != 0:

        run_result = queue.get()

        print(f"[{num_jobs}] Processing results of job q={run_result.q} m={run_result.m} n={run_result.n}")

        query = "INSERT INTO runs(q, m, n, wallclock_time, nr_of_poss_columns, nr_of_diff_columns, count_of_configurations, count_of_recursive_calls) VALUES (?, ?, ?, ?, ?, ?, ?, ?);"
        con.execute(query, (run_result.q, run_result.m, run_result.n, run_result.wallclock_time, run_result.nr_of_poss_columns, run_result.nr_of_diff_columns, run_result.count_of_configurations, run_result.count_of_recursive_calls))

        query = "INSERT INTO counts(q, m, n, d, count) VALUES (?, ?, ?, ?, ?);"
        for (d, count) in enumerate(run_result.distance_counts):
            #print(run_result.q, run_result.m, run_result.n, d, str(count))
            con.execute(query, (run_result.q, run_result.m, run_result.n, d, str(count)))
        con.commit()

        num_jobs -= 1

    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
