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
        self.t0 = isa.temperature(Q_(0, "m")).m

    def horner(
        self, coeffs: List[Number], x: Number, smallest_order: int = 0
    ) -> Number:
        return functools.reduce(
            lambda acc, coef: acc * x + coef, reversed(coeffs), 0
        ) * pow(x, smallest_order)

    def get_dp(self, alt):

        t_0 = self.t0
        p_0 = isa.pressure(Q_(0, "m")).m
        t_trop = isa.temperature(isa.H_TROP).m

        if alt < isa.H_TROP.m:
            # bada (2.2-18)
            dP = (
                p_0
                * pow(
                    (isa.temperature(Q_(alt, "m")).m / t_0),
                    (-isa.G_0.m / (isa.BETA_T.m * isa.R.m)),
                )
                / p_0
            )
        else:
            # bada (2.2-20)
            p_trop = p_0 * pow((t_trop / t_0), (-isa.G_0.m / (isa.BETA_T.m * isa.R.m)))
            dP = (
                p_trop
                * exp(-(isa.G_0.m / (isa.R.m * t_trop)) * (alt - isa.H_TROP.m))
                / p_0
            )
        return dP

    @default_units(tas="kts", alt="ft")
    def takeoff(self, tas, alt=Q_(0, "m"), temp=None):
        temp = temp if temp is not None else self.t0
        return self.climb(tas, alt, temp)

    @default_units(tas="kts", alt="ft")
    def cruise(self, tas, alt, temp=None, mass=None):
        mass = (
            Q_(mass, "kg")
            if mass is not None
            else self.bada_reader.get_param("ALM/DLM/MTOW")
        )
        D = self.drag.clean(mass, tas, alt)
        return D

    @default_units(tas="kts", alt="ft")
    def climb(self, tas, alt, roc=Q_(0, "m / s"), temp=None):
        temp = temp if temp is not None else self.t0
        dP = self.get_dp(alt.m)
        ct = self.climb_ct(tas.m, alt.m, dP, temp)
        return self.get_thrust(ct, alt, dP)

    @default_units(tas="kts", alt="ft")
    def descent_idle(self, tas, alt):
        mach = aero.tas2mach(Q_(tas, "kts"), Q_(alt, "m")).m
        dP = self.get_dp(alt.m)
        ti_coeff = self.bada_reader.get_param("PFM/TFM/LIDL/CT").m
        d_coef = []
        for m_coef in [ti_coeff[i : i + 4] for i in range(0, len(ti_coeff), 4)]:
            d_coef.append(self.horner(m_coef, dP, -1))
        ct = self.horner(d_coef, mach)
        return self.get_thrust(ct, alt, dP)

    def compute_throttle(self, tas, alt, dP, temp=None):
        temp = temp if temp is not None else self.t0
        kink = self.bada_reader.get_param("PFM/TFM/MCMB/kink").m
        mach = aero.tas2mach(Q_(tas, "kts"), Q_(alt, "m")).m
        dT = temp / self.t0

        if (temp - isa.temperature(Q_(alt, "m")).m) > kink:
            # temp-rated area: eq (3.3-6)
            ratio = dT * (1 + mach * mach * ((isa.GAMMA.m - 1) / 2))
            th_coeff = self.bada_reader.get_param("PFM/TFM/MCMB/temp_rating").m

            d_coef = []
            for m_coef in [th_coeff[i : i + 5] for i in range(0, len(th_coeff), 5)]:
                d_coef.append(self.horner(m_coef, mach))

            throttle_temp = self.horner(d_coef[:5], ratio)
            throttle_pres = self.horner(d_coef[5:], dP)
            throttle = throttle_temp + throttle_pres
        else:
            # flat-rated area: eq (3.3-5)
            th_coeff = self.bada_reader.get_param("PFM/TFM/MCMB/flat_rating").m
            d_coef = []
            for m_coef in [th_coeff[i : i + 6] for i in range(0, len(th_coeff), 6)]:
                d_coef.append(self.horner(m_coef, mach))

            throttle = self.horner(d_coef, dP)
        return throttle

    def climb_ct(self, tas, alt, dP, temp=None):
        temp = temp if temp is not None else self.t0
        mach = aero.tas2mach(Q_(tas, "kts"), Q_(alt, "m")).m
        throttle = self.compute_throttle(tas, alt, dP, temp)

        # non-idle rating: eq (3.3-3)
        ct_coeff = self.bada_reader.get_param("PFM/TFM/CT").m
        d_coef = []
        for m_coef in [ct_coeff[i : i + 6] for i in range(0, len(ct_coeff), 6)]:
            d_coef.append(self.horner(m_coef, mach))
        ct = self.horner(d_coef, throttle)
        return ct

    @default_units(tas="kts", alt="ft", temp="K")
    def idle_ct(self, tas, alt, roc=Q_(0, "m / s"), temp=None):
        temp = temp if temp is not None else self.t0
        # TODO
        print("to implement")
        return

    def get_thrust(self, ct, alt, dP):

        mref = self.bada_reader.get_param("PFM/MREF").m
        w_mref = mref * isa.G_0.m
        T = dP * w_mref * ct

        return Q_(T, "N")
