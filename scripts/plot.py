import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

STYLE = {
    ("starpu", "lws"): ("#1f77b4", "o", "-"),
    ("starpu", "dmda"): ("#1f77b4", "s", "--"),
    ("parsec", "lfq"): ("#d62728", "^", "-"),
    ("parsec", "ll"): ("#d62728", "D", "--"),
}


def main():
    csv_path = "data/chameleon-runtimes_781234/results.csv"
    out_dir = "plots"

    df = pd.read_csv(csv_path)
    df = df[df["exit_status"] == 0]

    algorithms = sorted(df["algorithm"].unique())
    fig, axes = plt.subplots(
        1, len(algorithms), figsize=(6.5 * len(algorithms), 4.5), sharey=False
    )
    if len(algorithms) == 1:
        axes = [axes]

    for ax, algo in zip(axes, algorithms):
        sub = df[df["algorithm"] == algo]
        agg = (
            sub.groupby(["runtime", "scheduler", "n"])["gflops"]
            .agg(
                median="median",
                q25=lambda s: s.quantile(0.25),
                q75=lambda s: s.quantile(0.75),
            )
            .reset_index()
        )
        for (runtime, sched), grp in agg.groupby(["runtime", "scheduler"]):
            grp = grp.sort_values("n")
            color, marker, linestyle = STYLE.get((runtime, sched), ("gray", "x", "-"))
            ax.plot(
                grp["n"],
                grp["median"],
                color=color,
                marker=marker,
                linestyle=linestyle,
                linewidth=1.8,
                markersize=6,
                label=f"{runtime}/{sched}",
            )
            ax.fill_between(
                grp["n"],
                grp["q25"],
                grp["q75"],
                color=color,
                alpha=0.18,
                linewidth=0,
            )

        ax.set_xticks(sorted(sub["n"].unique()))
        ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
        ax.set_xlabel("N (matriz N×N)")
        ax.set_ylabel("GFLOPS")
        ax.set_title(f"{algo}  (b=320, t=24, 5 reps)")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)
        ax.legend(loc="lower right", framealpha=0.9)

    fig.suptitle(
        "Chameleon CPU-only na cei (2× Xeon Silver 4116)",
        y=1.02,
    )
    fig.tight_layout()

    out_path = f"{out_dir}/gflops_vs_n.png"
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    out_pdf = f"{out_dir}/gflops_vs_n.pdf"
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"saved: {out_path}")
    print(f"saved: {out_pdf}")


if __name__ == "__main__":
    main()
