import logging
import os

import pytest
import vcr


# set vcr logging level
logging.basicConfig()
logger = logging.getLogger("vcr")
logger.setLevel(logging.DEBUG)

# set default location for vcr cassettes
test_dir = os.path.dirname(__file__)
test_data_dir = os.path.join(test_dir, "data", "cassettes")

# initialize vcr
vcr.default_vcr = vcr.VCR(
    cassette_library_dir=test_data_dir,
    filter_headers=["Authorization"],
    filter_post_data_parameters=['Authorization'],
    record_mode=os.environ.get("VCR_RECORD_MODE", "once"),
    )
vcr.use_cassette = vcr.default_vcr.use_cassette
