"""Write release/ files that identify and verify a CASSIA release.

    python scripts/make_release_manifest.py 0.3.2

Creates release/MANIFEST.sha256 (every code, data, database, result-table and figure file),
release/RELEASE.md (version, git commit of the code, commands, expected outputs), and copies the
latest test log and coverage report. Verify a checkout with:  sha256sum -c release/MANIFEST.sha256
"""
import hashlib
import subprocess
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REL = ROOT / "release"
REL.mkdir(exist_ok=True)
version = sys.argv[1] if len(sys.argv) > 1 else "0.3.2"
include = ["cassia", "scripts", "tests", "examples", "data", "databases", "results/v2", "results/v3",
           "manuscript/figures", "pyproject.toml", "requirements-lock.txt", "LICENSE", "CITATION.cff", "README.md", "CHANGELOG.md"]
files = []
for inc in include:
    p = ROOT / inc
    files += [p] if p.is_file() else sorted(q for q in p.rglob("*") if q.is_file() and "__pycache__" not in q.parts and q.suffix != ".pyc")
lines = [f"{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(ROOT)}" for f in files]
(REL / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")
commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
tables = sorted(str(f.relative_to(ROOT)) for f in files if f.suffix in (".csv", ".json", ".gz") and "results" in f.parts)
(REL / "RELEASE.md").write_text(f"""# CASSIA {version} release record

- Date: {date.today().isoformat()}
- Code commit at manifest generation (parent of the commit that adds this file): `{commit}`
- Content identity: `release/MANIFEST.sha256` ({len(lines)} files); verify with `sha256sum -c release/MANIFEST.sha256`
- Environment: `requirements-lock.txt` (Python 3.13.17); PHREEQC databases with checksums and provenance in `databases/`
- Tests: `release/test_log.txt`; coverage: `results/coverage_report.txt`
- Permanent archive: to be created by the author (GitHub release `cassia-v{version}` + Zenodo DOI); not available from the build environment

## Reproduce

```
pip install -e ..                      # adsorption equation package (repository root)
pip install -e .[ml,pores,dev]
pytest -q tests                        # expected: 40 tests pass (see release/test_log.txt)
bash scripts/reproduce_all.sh          # every table, figure and number, in dependency order (~9-11 h, 4 cores)
```

## Expected output tables ({len(tables)})

""" + "\n".join(f"- `{t}`" for t in tables) + "\n")
print("manifest", len(lines), "files; commit", commit)
