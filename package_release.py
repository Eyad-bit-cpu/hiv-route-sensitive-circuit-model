#!/usr/bin/env python
"""Create the canonical Model v1.2.1 source + results release ZIP.

The release intentionally includes the version-controlled source, frozen parameter/configuration
files, tests, documentation, and canonical regenerated results used by the manuscript.
Transient caches and compiled Python artifacts are excluded.
"""

import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT.parent / "hiv_circuit_sim_v1.2.1_release.zip"


def create_zip():
    """Create the canonical cross-platform release ZIP using forward-slash paths."""
    include_dirs = [
        "src", "scripts", "configs", "parameters", "data", "docs", "tests", "results"
    ]
    include_files = [
        "requirements.txt", "environment.yml", "README.md", "LICENSE",
        "pytest.ini", ".gitignore", "package_release.py",
        "REPAIR_REPORT_v1.2.md", "RECOVERY_PROTOCOL_CORRECTION_v1.2.1.md",
        "AUDIT_HARDENING_v1.2.1.md", "REGENERATED_RESULTS_v1.2.1.md",
        "FIGURE8_PATCH_NOTE.md", "VERIFICATION.txt"
    ]

    file_count = 0
    with zipfile.ZipFile(OUTPUT, "w", zipfile.ZIP_DEFLATED) as zf:
        for dir_name in include_dirs:
            dir_path = ROOT / dir_name
            if not dir_path.is_dir():
                continue
            for file_path in sorted(dir_path.rglob("*")):
                if not file_path.is_file():
                    continue
                rel_path = file_path.relative_to(ROOT)
                if any(part in {"__pycache__", ".pytest_cache"} for part in rel_path.parts):
                    continue
                if file_path.suffix in {".pyc", ".pyo"}:
                    continue
                zf.write(file_path, str(rel_path).replace("\\", "/"))
                file_count += 1

        for file_name in include_files:
            file_path = ROOT / file_name
            if file_path.is_file():
                zf.write(file_path, file_name.replace("\\", "/"))
                file_count += 1

    print(f"Created canonical release ZIP: {OUTPUT}")
    print(f"Total files: {file_count}")
    print(f"Total size: {OUTPUT.stat().st_size / (1024*1024):.2f} MB")


if __name__ == "__main__":
    create_zip()
