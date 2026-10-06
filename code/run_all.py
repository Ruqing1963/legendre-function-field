"""
Reproduce everything: tests, data, tables, figures.

    python code/run_all.py            # full run (about 2-3 hours on an 8-core laptop)
    python code/run_all.py --quick    # tests + tables/figures from the committed CSVs only

The full run includes the open-range grid (about 30 minutes) and the single largest case
q = 5, d = 8 (about 1 hour), run as  python code/open_regime.py --only 5,8 --full --max-tests 2e10

Then compile the paper:
    cd paper && pdflatex Legendre_FF_v7.tex && pdflatex Legendre_FF_v7.tex
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS_FULL = [["test_fflib.py"], ["test_fqlib.py"], ["thresholds.py"], ["interval_counts.py"],
              ["twisted_identity.py"], ["exhaustive_small.py"], ["open_regime.py"],
              ["open_regime.py", "--only", "5,8", "--full", "--max-tests", "2e10"],
              ["make_tables.py"], ["make_figures.py"]]
STEPS_QUICK = [["test_fflib.py"], ["test_fqlib.py"], ["make_tables.py"], ["make_figures.py"]]


def main():
    steps = STEPS_QUICK if "--quick" in sys.argv else STEPS_FULL
    for s in steps:
        print(f"==> {' '.join(s)}", flush=True)
        subprocess.run([sys.executable, os.path.join(HERE, s[0])] + s[1:], check=True)
    print("all steps finished")


if __name__ == "__main__":
    main()
