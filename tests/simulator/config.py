project_root = "@PHARE_PROJECT_DIR@"

at = "@"  # try to work without cmake config

if project_root == f"{at}PHARE_PROJECT_DIR{at}":
    from pathlib import Path

    project_root = Path(__file__).resolve().parent.parent.parent
