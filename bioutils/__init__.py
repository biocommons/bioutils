import pkg_resources
import re
import warnings


try:
    __version__ = pkg_resources.get_distribution(__name__).version
    if re.match(r"^\d+\.\d+\.\d+$", __version__) is not None:
        _is_released_version = True  # pragma: no cover
except pkg_resources.DistributionNotFound as e:  # pragma: no cover
    warnings.warn("can't get __version__ because %s package isn't installed" %
                  __package__, Warning)
    __version__ = None


from ._versionwarning import warnings
