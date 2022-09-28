from typing import List
from .. import config, config_file
import xml.etree.ElementTree as ET
from pathlib import Path
from .acm_params import acm_params
from pitot import Q_
import numpy as np


class BadaReader(object):
    def __init__(self, aircraft, **kwargs):

        self.bada_dir = Path(config.get("bada", "bada_path", fallback=""))
        self.file_path = self.bada_dir / f"{aircraft}/{aircraft}.xml"
        self.tree = ET.parse(self.file_path)

    def get_param(self, path: str) -> Q_ | List[Q_] | str:

        element = self.tree.find(path)

        if (unit := acm_params[path]) is not None:
            param = Q_(float(element.text), unit)
        else:
            children = [child.text for child in element]
            if children:
                param = Q_(np.array(children).astype(float))
            else:
                try:
                    param = Q_(float(element.text))
                except ValueError:
                    param = element.text

        return param
