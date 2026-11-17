"""Tests of pharesee.amr. Importing this package points the on-disk index cache to a
temporary directory (unless PHARESEE_CACHE_DIR is already set), so tests never write
into the user's ~/.cache/pharesee."""

import atexit
import os
import shutil
import tempfile

if not os.environ.get("PHARESEE_CACHE_DIR"):
    _cache = tempfile.mkdtemp(prefix="pharesee-test-cache-")
    os.environ["PHARESEE_CACHE_DIR"] = _cache
    atexit.register(shutil.rmtree, _cache, ignore_errors=True)
