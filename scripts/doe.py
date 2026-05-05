import csv
import random
from itertools import product

RUNTIME_SCHED = [
    ("starpu", "lws"),
    ("starpu", "dmda"),
    ("parsec", "lfq"),
    ("parsec", "ll"),
]

ALGORITHMS = ["potrf", "geqrf"]
BLOCK_SIZE = 320
N_VALUES = [1920, 3840, 7680, 11520, 15360, 19200]
THREADS = 24
REPS = 5
SEED = 42
FIELDS = ["runtime", "scheduler", "algorithm", "n", "b", "threads", "rep"]


def build_design():
    rows = []
    for (runtime, scheduler), algorithm, n_value, rep in product(
        RUNTIME_SCHED, ALGORITHMS, N_VALUES, range(1, REPS + 1)
    ):
        rows.append(
            {
                "runtime": runtime,
                "scheduler": scheduler,
                "algorithm": algorithm,
                "n": n_value,
                "b": BLOCK_SIZE,
                "threads": THREADS,
                "rep": rep,
            }
        )

    random.Random(SEED).shuffle(rows)
    return rows


def main():
    out_path = "experiments.csv"

    rows = build_design()
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in FIELDS})


if __name__ == "__main__":
    main()
