import pkg_resources
import re
import warnings

import sys


try:
    __version__ = pkg_resources.get_distribution(__name__).version
    if re.match('^\d+\.\d+\.\d+$', __version__) is not None:
        _is_released_version = True  # pragma: no cover
except pkg_resources.DistributionNotFound as e:  # pragma: no cover
    warnings.warn("can't get __version__ because %s package isn't installed" %
                  __package__, Warning)
    __version__ = None


if sys.version_info < (3,5):
    warnings.warn("Support for Python <3.6 will be dropped on 2019-03-31. See https://github.com/biocommons/org/wiki/Migrating-to-Python-3.6")
