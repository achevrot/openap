import configparser
from pathlib import Path

from appdirs import user_cache_dir, user_config_dir

from .extra import aero
from .extra import nav
from .extra import filters
from .extra import statistics

from .thrust import Thrust
from .drag import Drag
from .fuel import FuelFlow
from .emission import Emission
from .kinematic import WRAP
from .phase import FlightPhase


__version__ = importlib_metadata.version("traffic")  # type: ignore
__all__ = ["config_dir", "config_file"]

config_dir = Path(user_config_dir("openap"))
config_file = config_dir / "openap.conf"

if not config_dir.exists():
    config_template = (Path(__file__).parent / "openap.conf").read_text()
    config_dir.mkdir(parents=True)
    config_file.write_text(config_template)

config = configparser.ConfigParser()
config.read(config_file.as_posix())
