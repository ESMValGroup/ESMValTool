"""Type definitions for CMORizers."""

from typing import TypedDict


class DatasetInfo(TypedDict):
    """Dataset description."""

    tier: int
    source: str
    last_access: str
    info: str
