from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import yaml


def _merge(base: dict, override: Mapping[str, Any] | None) -> dict:
    if not override:
        return dict(base)
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _merge(out[key], value)
        else:
            out[key] = value
    return out


@dataclass
class PipelineConfig:
    raw: dict = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "PipelineConfig":
        return cls(raw=dict(data))

    @classmethod
    def from_yaml(cls, path: str | Path) -> "PipelineConfig":
        with open(path, "r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
        return cls.from_dict(data)

    def get(self, *keys, default=None):
        cur: Any = self.raw
        for key in keys:
            if not isinstance(cur, dict) or key not in cur:
                return default
            cur = cur[key]
        return cur

    def updated(self, override: Mapping[str, Any]) -> "PipelineConfig":
        return PipelineConfig(_merge(self.raw, override))

    @property
    def seed(self) -> int:
        return int(self.raw.get("seed", 42))

    @property
    def output_dir(self) -> Path:
        return Path(self.raw.get("output_dir", "outputs"))


def load_config(path: str | Path | None = None, override: Mapping[str, Any] | None = None) -> PipelineConfig:
    if path is None:
        cfg = PipelineConfig()
    else:
        cfg = PipelineConfig.from_yaml(path)
    if override:
        cfg = cfg.updated(override)
    return cfg
