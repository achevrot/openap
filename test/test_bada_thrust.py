#%%

import pytest
import pandas as pd
from pathlib import Path
from openap import config, config_file
from openap.bada4 import Thrust
from pitot import Q_


@pytest.fixture
def thrust_table():
    def _method(
        aircraft: str = "A320-232", flight_phase: str = "CLIMB"
    ) -> pd.DataFrame:

        bada_dir = Path(config.get("bada", "bada_path", fallback=""))
        file_path = bada_dir / f"{aircraft}/{aircraft}_ISA.PTD"

        if flight_phase == "CLIMB":
            return pd.read_fwf(
                file_path, infer_nrows=20, skiprows=list(range(6)) + [8]
            )[:27]
        if flight_phase == "DESC":
            return pd.read_fwf(
                file_path, infer_nrows=20, skiprows=list(range(6)) + [8]
            )[93:120]
        if flight_phase == "CRUISE":
            return pd.read_fwf(
                file_path, infer_nrows=20, skiprows=list(range(6)) + [8]
            )[186:213]

    return _method


def test_climb(
    thrust_table: pd.DataFrame,
) -> None:
    th = Thrust("A320-232")
    for _, row in thrust_table().iterrows():
        res = th.climb(
            Q_(float(row.TAS), "kt"),
            Q_(float(row.FL + "00"), "ft"),
            roc=float(row.ROCD),
        ).m
        assert res == pytest.approx(float(row.Thrust), abs=1)


def test_cruise(
    thrust_table: pd.DataFrame,
) -> None:
    th = Thrust("A320-232")
    for _, row in thrust_table(flight_phase="CRUISE").iterrows():
        res = th.cruise(
            Q_(float(row.TAS), "kt"),
            Q_(float(row.FL + "00"), "ft"),
            roc=float(row.ROCD),
            mass=float(row.mass),
        ).m
        assert res == pytest.approx(float(row.Thrust), abs=5)


def test_desc(
    thrust_table: pd.DataFrame,
) -> None:
    th = Thrust("A320-232")
    for _, row in thrust_table(flight_phase="DESC").iterrows():
        res = th.descent_idle(
            Q_(float(row.TAS), "kt"),
            Q_(float(row.FL + "00"), "ft"),
            # roc=float(row.ROCD),
            #  mass=float(row.mass),
        ).m
        assert res == pytest.approx(float(row.Thrust), abs=5)


def test_thrust(
    thrust_table: pd.DataFrame,
) -> None:
    th = Thrust("A320-232")
    res = th.descent_idle(
        Q_(float(126.06), "kt"),
        Q_(float("5" + "00"), "ft"),
        # roc=float(row.ROCD),
        #  mass=float(row.mass),
    ).m
    assert res == pytest.approx(float(35457), abs=1)


# %%
