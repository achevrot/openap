import functools
from math import cos, exp, pow
from numbers import Number
from typing import List
import math


# from . import BadaReader
# from . import Drag
# from . import acm_params

from pitot.wrapper import default_units
from pitot import Q_, ureg
import pitot.isa as isa
import pitot.aero as aero

from openap.bada4.bada_reader import BadaReader
from openap.bada4.drag import Drag


class Thrust(object):
    def __init__(self, aircraft_type, eng=None, **kwargs):

        self.bada_reader = BadaReader.getAircraft(aircraft_type)
        self.drag = Drag(aircraft_type, **kwargs)
        self.t0 = isa.temperature(Q_(0, "m"))

    def horner(
        self, coeffs: List[Number], x: Number, smallest_order: int = 0
    ) -> Number:
        return functools.reduce(
            lambda acc, coef: acc * x + coef, reversed(coeffs), 0
        ) * pow(x, smallest_order)

    @default_units(alt="m")
    def get_dp(self, alt):

        t_0 = self.t0
        p_0 = isa.pressure(Q_(0, "m"))
        t_trop = isa.temperature(isa.H_TROP)

        if alt < isa.H_TROP:
            # bada (2.2-18)
            dP = (
                p_0
                * (isa.temperature(alt) / t_0) ** (-isa.G_0 / (isa.BETA_T * isa.R))
                / p_0
            )
        else:
            # bada (2.2-20)
            p_trop = p_0 * (t_trop / t_0) ** (-isa.G_0 / (isa.BETA_T * isa.R))
            dP = p_trop * exp(-(isa.G_0 / (isa.R * t_trop)) * (alt - isa.H_TROP)) / p_0
        return dP

    @default_units(tas="kts", alt="ft")
    def takeoff(self, tas, alt=Q_(0, "m"), temp=None):
        temp = temp if temp is not None else isa.temperature(alt)
        return self.climb(tas, alt, temp)

    @default_units(tas="kts", alt="ft")
    def cruise(self, tas, alt, roc=None, temp=None, mass=None):
        mass = (
            Q_(mass, "kg")
            if mass is not None
            else self.bada_reader.get_param("ALM/DLM/MTOW")
        )
        D = self.drag.clean(mass, tas, alt)
        return D

    @default_units(tas="kts", alt="m")
    def climb(self, tas, alt, roc=Q_(0, "m / s"), temp=None):
        temp = temp if temp is not None else isa.temperature(alt)
        dP = self.get_dp(alt)
        ct = self.climb_ct(tas, alt, dP, temp)
        return self.get_thrust(ct, dP)

    @default_units(tas="kts", alt="ft")
    def descent_idle(self, tas, alt):
        mach = aero.tas2mach(Q_(tas, "kts"), Q_(alt, "m"))
        dP = self.get_dp(alt)

        if alt > Q_(1000, "ft"):

            ti_coeff = self.bada_reader.get_param("PFM/TFM/LIDL/CT")
            d_coef = []
            for m_coef in [ti_coeff[i : i + 4] for i in range(0, len(ti_coeff), 4)]:
                d_coef.append(self.horner(m_coef, dP, -1))
            ct = self.horner(d_coef, mach)
        else:
            # TODO : Find how to calculate the landing thrust
            throttle = self.bada_reader.get_param("PFM/TFM/throttle/low")
            ct = self.climb_ct(tas, alt, dP, throttle=throttle)
        return self.get_thrust(ct, dP)

    @default_units(tas="kts", alt="m")
    def compute_throttle(self, tas, alt, dP, temp=None):
        temp = temp if temp is not None else isa.temperature(alt)
        kink = self.bada_reader.get_param("PFM/TFM/MCMB/kink")
        mach = aero.tas2mach(Q_(tas, "kts"), Q_(alt, "m"))
        dT = temp / self.t0

        if (temp - isa.temperature(Q_(alt, "m"))) > kink:
            # temp-rated area: eq (3.3-6)
            ratio = dT * (1 + mach * mach * ((isa.GAMMA - 1) / 2))
            th_coeff = self.bada_reader.get_param("PFM/TFM/MCMB/temp_rating")

            d_coef = []
            for m_coef in [th_coeff[i : i + 5] for i in range(0, len(th_coeff), 5)]:
                d_coef.append(self.horner(m_coef, mach))

            throttle_temp = self.horner(d_coef[:5], ratio)
            throttle_pres = self.horner(d_coef[5:], dP)
            throttle = throttle_temp + throttle_pres
        else:
            # flat-rated area: eq (3.3-5)
            th_coeff = self.bada_reader.get_param("PFM/TFM/MCMB/flat_rating")
            d_coef = []
            for m_coef in [th_coeff[i : i + 6] for i in range(0, len(th_coeff), 6)]:
                d_coef.append(self.horner(m_coef, mach))

            throttle = self.horner(d_coef, dP)
        return throttle

    @default_units(tas="kts", alt="ft", temp="K")
    def climb_ct(self, tas, alt, dP, temp=None, throttle=None):
        temp = temp if temp is not None else isa.temperature(alt)
        mach = aero.tas2mach(tas, alt)
        throttle = (
            throttle
            if throttle is not None
            else self.compute_throttle(tas, alt, dP, temp)
        )

        # non-idle rating: eq (3.3-3)
        ct_coeff = self.bada_reader.get_param("PFM/TFM/CT")
        d_coef = []
        for m_coef in [ct_coeff[i : i + 6] for i in range(0, len(ct_coeff), 6)]:
            d_coef.append(self.horner(m_coef, mach))
        ct = self.horner(d_coef, throttle)
        return ct

    @default_units(tas="kts", alt="ft", temp="K")
    def idle_ct(self, tas, alt, roc=Q_(0, "m / s"), temp=None):
        temp = temp if temp is not None else isa.temperature(alt)
        # TODO
        print("to implement")
        return

    def get_thrust(self, ct, dP):

        mref = self.bada_reader.get_param("PFM/MREF")
        w_mref = mref * isa.G_0
        T = dP * w_mref * ct

        return T.to("N")
