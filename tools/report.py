import os
import sys
import shutil
import tempfile
import datetime
from pathlib import Path
from zipfile import ZipFile

from tools.python3 import pushd

SEARCH_PATHS = list(set([os.getcwd()] + sys.path))


def find_phare_dirs():
    phare_dirs = []
    for path in SEARCH_PATHS:
        dot_phare_dir = Path(path) / ".phare"
        if dot_phare_dir.exists():
            phare_dirs += [dot_phare_dir]
    return phare_dirs


def main():
    phare_dirs = find_phare_dirs()
    assert find_phare_dirs
    cwd = Path(os.getcwd())
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmp_dir = Path(tmpdirname)
        for dir_idx, phare_dir in enumerate(phare_dirs):
            with pushd(phare_dir.parent):
                shutil.copytree(phare_dir.stem, tmp_dir / phare_dir.stem)
                shutil.move(
                    tmp_dir / phare_dir.stem, f"{tmpdirname}/{phare_dir.stem}_{dir_idx}"
                )

        zip_name = "PHARE_REPORT.zip"
        final_report_zip = cwd / zip_name
        if final_report_zip.exists():
            """move existing to subdirectory with creation timestamp in name"""
            reports_dir = cwd / ".phare_reports"
            reports_dir.mkdir(exist_ok=True, parents=True)
            timestamp = datetime.datetime.fromtimestamp(
                final_report_zip.stat().st_ctime
            ).isoformat()
            shutil.move(
                final_report_zip,
                reports_dir / f"{final_report_zip.stem}.{timestamp}.zip",
            )

        tmp_report_zip = tmp_dir / zip_name
        with pushd(tmpdirname):
            """shutil.make_archive doesn't like multiple directory inputs"""
            with ZipFile(tmp_report_zip, "w") as zip_object:
                for phare_dir in Path(tmpdirname).glob(".phare_*"):
                    for item in phare_dir.glob("**/*"):
                        zip_object.write(item, arcname=item.relative_to(tmpdirname))
            shutil.move(tmp_report_zip, cwd)


if __name__ == "__main__":
    main()
