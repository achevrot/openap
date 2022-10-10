from functools import lru_cache
from pathlib import Path
from typing import List

import numpy as np
from lxml import etree
from pitot import Q_

from .. import config, config_file
from .acm_params import acm_params


class BadaReader(object):
    INSTANCES: dict[str, "BadaReader"] = {}

    @classmethod
    def getAircraft(cls, aircraft: str) -> "BadaReader":
        if instance := cls.INSTANCES.get(aircraft, None):
            return instance

        instance = BadaReader()
        instance.bada_dir = Path(config.get("bada", "bada_path", fallback=""))
        instance.file_path = instance.bada_dir / f"{aircraft}/{aircraft}.xml"
        instance.tree = etree.parse(instance.file_path)
        return instance

    @lru_cache()
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
