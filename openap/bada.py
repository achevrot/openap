from . import config, config_file
import xml.etree.ElementTree as ET
from pathlib import Path


class Bada(object):
    def __init__(self, aircraft, **kwargs):

        self.bada_dir = Path(config.get("bada", "bada_path", fallback=""))
        self.file_path = self.bada_dir / f"{aircraft}/{aircraft}.xml"
        self.tree = ET.parse(self.file_path)

    def get_param(self, path: str, findall: str = None) -> str:
        if findall is not None:
            return self.tree.find(path).findall(findall)
        else:
            return self.tree.find(path)
