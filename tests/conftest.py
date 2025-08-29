import logging
import os
from pathlib import Path

import vcr

# set vcr logging level
logging.basicConfig()
logger = logging.getLogger("vcr")

# set default location for vcr cassettes
test_dir = Path(__file__).parent
test_data_dir = test_dir / "data" / "cassettes"

# initialize vcr
vcr.default_vcr = vcr.VCR(
    cassette_library_dir=str(test_data_dir),
    filter_headers=["Authorization"],
    filter_post_data_parameters=["Authorization"],
    record_mode=os.environ.get("VCR_RECORD_MODE", "once"),
)
vcr.use_cassette = vcr.default_vcr.use_cassette
