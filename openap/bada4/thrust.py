import functools

from numbers import Number
from typing import List
import math
from typing import Any, Tuple

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

        self.aircraft_type = aircraft_type
        self.bada_reader = BadaReader(self.aircraft_type)
        self.drag = Drag(aircraft_type, **kwargs)

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

    @default_units(tas="kts", alt="ft", temp="K")
    def takeoff(self, tas, alt=None, temp=None):
        return self.climb(tas, alt=None, temp=None)

    @default_units(tas="kts", alt="ft", temp="K", mass="kg")
    def cruise(self, tas, alt, temp=None, mass=None):
        mass = (
            mass if mass is not None else self.bada_reader.get_param("ACM/ALM/DLM/MTOW")
        )
        drag = Drag(self.aircraft_type)
        D = drag.clean(mass, tas, alt)
        return D

    @default_units(tas="kts", alt="ft", roc="m / s", temp="K")
    def climb(self, tas, alt, roc=0, temp=None):
        temp = temp if temp is not None else isa.temperature(alt)
        ct = self.climb_ct(tas, alt, temp)

        return self.get_thrust(ct, alt)

    @default_units(tas="kts", alt="ft", temp="K")
    def descent_idle(self, tas, alt):
        mach = aero.tas2mach(tas, alt)
        dP = self.get_dp(alt)
        ti_coeff = self.bada_reader.get_param("PFM/TFM/LIDL/CT")
        d_coef = []
        for m_coef in [ti_coeff[i : i + 4] for i in range(0, len(ti_coeff), 4)]:
            d_coef.append(self.horner(m_coef, dP, -1))
        ct = self.horner(d_coef, mach)

        return self.get_thrust(ct, alt)

    @default_units(tas="kts", alt="ft", temp="K")
    def compute_throttle(self, tas, alt, temp=None):
        kink = self.bada_reader.get_param("PFM/TFM/MCMB/kink")
        mach = aero.tas2mach(tas, alt)
        temp = temp if temp is not None else isa.temperature(alt)
        dT = temp / isa.temperature(Q_(0, "m"))
        dP = self.get_dp(alt)

        if (temp - isa.temperature(alt)) > kink:
            # temp-rated area: eq (3.3-6)
            ratio = dT * (1 + mach**2 * ((isa.GAMMA - 1) / 2))
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
    def climb_ct(self, tas, alt, temp=None):

        mach = aero.tas2mach(tas, alt)
        temp = temp if temp is not None else isa.temperature(alt)
        throttle = self.compute_throttle(tas, alt, temp)

        # non-idle rating: eq (3.3-3)
        ct_coeff = self.bada_reader.get_param("PFM/TFM/CT")
        d_coef = []
        for m_coef in [ct_coeff[i : i + 6] for i in range(0, len(ct_coeff), 6)]:
            d_coef.append(self.horner(m_coef, mach))
        ct = self.horner(d_coef, throttle)
        return ct

    @default_units(tas="kts", alt="ft", roc="m / s", temp="K")
    def idle_ct(self, tas, alt, roc=0, temp=None):
        # TODO
        print("to implement")
        return

    @default_units(alt="ft")
    def get_thrust(self, ct, alt):

        dP = self.get_dp(alt)

        mref = self.bada_reader.get_param("PFM/MREF")
        w_mref = mref * isa.G_0
        T = dP * w_mref * ct

        return T.to("N")
