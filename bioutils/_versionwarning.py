"""emits a warning when imported under Python < 3.6

This module may be used by other biocommons packages

"""

import logging
import sys
import warnings

version_warning = ("biocommons tools are tested and supported only on Python >= 3.6"
                   " (https://github.com/biocommons/org/wiki/Migrating-to-Python-3.6)")

_logger = logging.getLogger(__package__)

if sys.version_info < (3, 6):
    _logger.warning(version_warning)
    
