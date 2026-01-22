"""Type definitions for CMORizers."""

import datetime
from typing import TypedDict


class DatasetInfo(TypedDict):
    """Dataset description."""

    tier: int
    source: str
    last_access: datetime.datetime
    info: str
