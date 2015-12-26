import pkg_resources
import re
import warnings

try:
    __version__ = pkg_resources.get_distribution(__name__).version
    if re.match('^\d+\.\d+\.\d+$', __version__) is not None:
        _is_released_version = True
except pkg_resources.DistributionNotFound as e:
    warnings.warn("can't get __version__ because %s package isn't installed" % __package__, Warning)
    __version__ = None
