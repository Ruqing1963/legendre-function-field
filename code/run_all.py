"""
Reproduce everything: tests, data, tables, figures.

    python code/run_all.py            # full run (about 10-20 minutes on a laptop)
    python code/run_all.py --quick    # tests + tables/figures from the committed CSVs only

Then compile the paper:
    cd paper && pdflatex Legendre_FF_v6.tex && pdflatex Legendre_FF_v6.tex
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS_FULL = ["test_fflib.py", "thresholds.py", "interval_counts.py", "twisted_identity.py",
              "exhaustive_small.py", "make_tables.py", "make_figures.py"]
STEPS_QUICK = ["test_fflib.py", "make_tables.py", "make_figures.py"]


def main():
    steps = STEPS_QUICK if "--quick" in sys.argv else STEPS_FULL
    for s in steps:
        print(f"==> {s}", flush=True)
        subprocess.run([sys.executable, os.path.join(HERE, s)], check=True)
    print("all steps finished")


if __name__ == "__main__":
    main()
