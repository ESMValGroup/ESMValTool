from __future__ import annotations

import argparse
from pathlib import Path

from rf_shap.pipeline import run_pipeline


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Diagnose model-vs-observation input/output with tree models + SHAP")
    parser.add_argument("-c", "--config", type=Path, default=Path("configs/default.yaml"))
    args = parser.parse_args(argv)
    result = run_pipeline(args.config)
    if isinstance(result, list):
        print(f"Completed {len(result)} folds")
        for i, item in enumerate(result, 1):
            print(f"fold {i}: {item.metrics}")
    else:
        print(result.metrics)


if __name__ == "__main__":
    main()
