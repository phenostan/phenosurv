"""Copy the Stan models from phenostan into src/phenosurv/stan/.

phenostan is the single source of truth for the models; phenosurv ships a
pinned copy so it works offline. Run through pixi:

    pixi run sync-models                 # latest phenostan main
    pixi run sync-models --ref v0.1.0    # a tag, branch or commit
    pixi run sync-models --source ../phenostan --ref my-branch

The .stan files in src/phenosurv/stan/ are replaced to mirror phenostan's
src/models/, and PHENOSTAN_VERSION records exactly where they came from.
"""

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

DEFAULT_SOURCE = "https://github.com/phenostan/phenostan.git"
MODELS_SUBDIR = Path("src") / "models"
STAN_DIR = Path(__file__).resolve().parent.parent / "src" / "phenosurv" / "stan"


def git(*args, cwd):
    result = subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True, text=True)
    return result.stdout.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", default=DEFAULT_SOURCE, help="phenostan git URL or local path")
    parser.add_argument("--ref", default="main", help="tag, branch or commit to copy from (default: main)")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        checkout = Path(tmp) / "phenostan"
        subprocess.run(["git", "clone", "--quiet", args.source, str(checkout)], check=True)
        git("checkout", "--quiet", args.ref, cwd=checkout)
        commit = git("rev-parse", "HEAD", cwd=checkout)

        stan_files = sorted((checkout / MODELS_SUBDIR).glob("*.stan"))
        if not stan_files:
            raise SystemExit(f"No .stan files found in {MODELS_SUBDIR} at {args.ref}")

        STAN_DIR.mkdir(parents=True, exist_ok=True)
        for old in STAN_DIR.glob("*.stan"):
            old.unlink()
        for stan_file in stan_files:
            shutil.copy2(stan_file, STAN_DIR / stan_file.name)

    (STAN_DIR / "PHENOSTAN_VERSION").write_text(
        f"source: {args.source}\nref: {args.ref}\ncommit: {commit}\n"
    )

    print(f"Copied {len(stan_files)} model(s) from phenostan {args.ref} ({commit[:7]}):")
    for stan_file in stan_files:
        print(f"  src/phenosurv/stan/{stan_file.name}")


if __name__ == "__main__":
    main()
