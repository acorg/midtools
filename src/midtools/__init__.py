from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("midtools")
except PackageNotFoundError:
    __version__ = "unknown"
