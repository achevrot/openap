import functools
import math
from numbers import Number
from re import A
from typing import Any, List
from openap.bada4.bada_reader import BadaReader
from pitot import Q_
from pitot.wrapper import default_units

import pitot.isa as isa
import pitot.aero as aero
from openap.bada4.drag import Drag

from openap.bada4.thrust import Thrust


class FuelFlow(object):
    """Fuel flow model based on ICAO emission databank."""

    def __init__(self, aircraft_type, eng=None, **kwargs):
        """Initialize FuelFlow object.

        Args:
            ac (string): ICAO aircraft type (for example: A320).
            eng (string): Engine type (for example: CFM56-5A3).
                Leave empty to use the default engine specified
                by in the aircraft database.
            polydeg (int): Order of the polynomials for fuel flow model (2 or 3), defaults to 2.

        """

        self.h_to = Q_(400, "ft")
        self.h_ic = Q_(2000, "ft")
        self.h_ap = Q_(8000, "ft")
        self.h_ld = Q_(3000, "ft")
        self.Cvmin_to = 1.13
        self.Cvmin = 1.23

        self.aircraft_type = aircraft_type
        self.thrust = Thrust(self.aircraft_type, eng, **kwargs)
        self.drag = Drag(self.aircraft_type, **kwargs)
        self.bada_reader = BadaReader(self.aircraft_type)

    def horner(
        self, coeffs: List[Number], x: Number, smallest_order: int = 0
    ) -> Number:
        return (
            functools.reduce(lambda acc, coef: acc * x + coef, reversed(coeffs), 0)
            * x**smallest_order
        )

    @default_units(alt="m")
    def get_dp(self, alt: Any = 0):

        if alt < Q_(11000, "m"):
            # bada (2.2-18)
            dP = (
                isa.pressure(Q_(0, "m"))
                * (isa.temperature(alt) / isa.temperature(Q_(0, "m")))
                ** (-isa.G_0 / (isa.BETA_T * isa.R))
                / isa.pressure(Q_(0, "m"))
            )
        else:
            # bada (2.2-20)
            p_trop = isa.pressure(Q_(0, "m")) * (
                isa.temperature(isa.H_TROP) / isa.temperature(Q_(0, "m"))
            ) ** (-isa.G_0 / (isa.BETA_T * isa.R))
            dP = (
                p_trop
                * math.exp(
                    -(isa.G_0 / (isa.R * isa.temperature(isa.H_TROP)))
                    * (alt - isa.H_TROP)
                )
                / isa.pressure(Q_(0, "m"))
            )
        return dP

    def get_phase(self, alt: float, tas: float, path_angle: int = 0) -> str:
        """Simplified version of phase calculation

        Args:
            alt (float): _description_
            tas (float): _description_
            path_angle (int, optional): _description_. Defaults to 0.

        Returns:
            str: _description_
        """

        if path_angle >= 0:
            if alt < self.h_to:
                phase = "TAKEOFF"
            elif alt < self.h_ic:
                phase = "CLIMB"
            else:
                phase = "CRUISE"

        else:
            if alt < self.h_ld:
                phase = "LANDING"
            elif alt < self.h_ap:
                phase = "DESCENT"
            else:
                phase = "CRUISE"

        return phase

    def at_thrust(self, acthr, alt=0, limit=True):
        """Compute the fuel flow at a given total thrust.

        Args:
            acthr (int or ndarray): The total net thrust of the aircraft (unit: N).
            alt (int or ndarray): Aircraft altitude (unit: ft).

        Returns:
            float: Fuel flow (unit: kg/s).

        """
        n_eng = self.aircraft["engine"]["number"]
        engthr = acthr / n_eng

        maxthr = self.thrust.takeoff(tas=0, alt=0)
        ratio = acthr / maxthr

        if limit:
            ratio = self.np.where(ratio < 0.07, 0.07, ratio)
            ratio = self.np.where(ratio > 1, 1, ratio)

        ff_sl = self.polyfuel(ratio)
        ff_corr_alt = self.engine["fuel_ch"] * (engthr / 1000) * (alt * 0.3048)
        ff_eng = ff_sl + ff_corr_alt

        fuelflow = ff_eng * n_eng

        return fuelflow

    def takeoff(self, tas, alt=None, throttle=1):
        """Compute the fuel flow at takeoff.

        The net thrust is first estimated based on the maximum thrust model
        and throttle setting. Then FuelFlow.at_thrust() is called to compted
        the thrust.

        Args:
            tas (int or ndarray): Aircraft true airspeed (unit: kt).
            alt (int or ndarray): Altitude of airport (unit: ft). Defaults to sea-level.
            throttle (float or ndarray): The throttle setting, between 0 and 1.
                Defaults to 1, which is at full thrust.

        Returns:
            float: Fuel flow (unit: kg/s).

        """
        Tmax = self.thrust.takeoff(tas=tas, alt=alt)
        fuelflow = throttle * self.at_thrust(Tmax)
        return fuelflow

    @default_units(mass="kg", tas="kts", alt="ft", path_angle="degree", temp="K")
    def enroute(self, mass, tas, alt, path_angle=0, limit=True, temp=None):
        """Compute the fuel flow during climb, cruise, or descent.

        The net thrust is first estimated based on the dynamic equation.
        Then FuelFlow.at_thrust() is called to compted the thrust. Assuming
        no flap deflection and no landing gear extended.

        Args:
            mass (int or ndarray): Aircraft mass (unit: kg).
            tas (int or ndarray): Aircraft true airspeed (unit: kt).
            alt (int or ndarray): Aircraft altitude (unit: ft).
            path_angle (float or ndarray): Flight path angle (unit: degrees).

        Returns:
            float: Fuel flow (unit: kg/s).

        """

        phase = self.get_phase(alt, tas, path_angle)
        mach = aero.tas2mach(tas, alt)
        mref = self.bada_reader.get_param("PFM/MREF")
        w_mref = mref * isa.G_0
        lhv = self.bada_reader.get_param("PFM/LHV")
        dP = self.get_dp(alt)
        temp = temp if temp is not None else isa.temperature(alt)
        dT = temp / isa.temperature(Q_(0, "m"))

        if phase == "CLIMB":
            ct = self.thrust.climb_ct(tas, alt)
            f_coef = self.bada_reader.get_param("PFM/TFM/CF")
            d_coef = []
            for m_coef in [f_coef[i : i + 5] for i in range(0, len(f_coef), 5)]:
                d_coef.append(self.horner(m_coef, ct.m))

            cf = self.horner(d_coef, mach)

        elif phase == "CRUISE":
            d = self.drag.clean(mass=mass, tas=tas, alt=alt, path_angle=path_angle)
            ct = d / (w_mref * dP)
            f_coef = self.bada_reader.get_param("PFM/TFM/CF")
            d_coef = []
            for m_coef in [f_coef[i : i + 5] for i in range(0, len(f_coef), 5)]:
                d_coef.append(self.horner(m_coef, ct.m))

            cf = self.horner(d_coef, mach)

        else:
            fi_coef = self.bada_reader.get_param("PFM/TFM/LIDL/CF")
            d_coef = []
            for m_coef in [fi_coef[i : i + 3] for i in range(0, len(fi_coef), 3)]:
                d_coef.append(self.horner(m_coef, dP))

            cf = self.horner(d_coef, mach) * math.sqrt(dT) / dP

        fuelflow = dP * math.sqrt(dT) * w_mref * isa.sound_speed(Q_(0, "ft")) * cf / lhv
        return fuelflow

    def plot_model(self, plot=True):
        """Plot the engine fuel model, or return the pyplot object.

        Args:
            plot (bool): Display the plot or return an object.

        Returns:
            None or pyplot object.

        """
        import matplotlib.pyplot as plt

        x = [0.07, 0.3, 0.85, 1.0]
        y = [
            self.engine["ff_idl"],
            self.engine["ff_app"],
            self.engine["ff_co"],
            self.engine["ff_to"],
        ]
        plt.scatter(x, y, color="k")

        xx = self.np.linspace(0, 1, 50)
        yy = self.polyfuel(xx)
        plt.plot(xx, yy, "--", color="gray")

        if plot:
            plt.show()
        else:
            return plt
