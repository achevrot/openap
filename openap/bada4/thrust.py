import functools
import importlib
from numbers import Number
from typing import List
import math
from . import BadaReader


class Thrust(object):
    def __init__(self, aircraft_type, eng=None, **kwargs):

        if not hasattr(self, "aero"):
            self.aero = importlib.import_module("openap").aero

        self.bada_reader = BadaReader(aircraft_type)

    def horner(coeffs: List[Number], x: Number, smallest_order: int = 0) -> Number:
        return (
            functools.reduce(lambda acc, coef: acc * x + coef, reversed(coeffs), 0)
            * x**smallest_order
        )

    def get_dp(self, alt=0):
        h = alt * self.aero.ft
        if h < 11000:
            # bada (2.2-18)
            dP = (
                self.aero.pressure(0)
                * (self.aero.temperature(h) / self.aero.temperature(0))
                ** (-self.aero.g0 / (self.aero.beta * self.aero.R))
                / self.aero.pressure(0)
            )
        else:
            # bada (2.2-20)
            p_trop = self.aero.pressure(0) * (
                self.aero.temperature(11000) / self.aero.temperature(0)
            ) ** (-self.aero.g0 / (self.aero.beta * self.aero.R))
            dP = (
                p_trop
                * math.exp(
                    -(self.aero.g0 / (self.aero.R * self.aero.temperature(11000)))
                    * (h - 11000)
                )
                / self.aero.pressure(0)
            )

        return dP

    def climb(self, tas, alt, roc=0, temp=None):

        roc = roc
        kink = int(self.bada_reader.get_param("PFM/TFM/MCMB/kink"))
        h = alt * self.aero.ft
        mach = self.aero.tas2mach(tas * self.aero.kts, h)

        if temp is not None:
            dT = temp / self.aero.temperature(0)
        else:
            dT = self.aero.temperature(h) / self.aero.temperature(0)

        dP = self.get_dp(alt)

        if dT > kink:
            # temp-rated area: eq (3.3-6)
            ratio = dT * (1 + mach**2 * ((self.aero.gamma - 1) / 2))
            th_coeff = self.bada_reader.get_param("PFM/TFM/MCMB/temp_rating", "c")
            th_coeff = [float(elem) for elem in th_coeff]

            d_coef = []
            for m_coef in [th_coeff[i : i + 5] for i in range(0, len(th_coeff), 5)]:
                d_coef.append(self.poly(m_coef, mach))

            throttle_temp = self.poly(d_coef[:5], ratio)
            throttle_pres = self.poly(d_coef[5:], dP)
            throttle = throttle_temp + throttle_pres

        else:
            # flat-rated area: eq (3.3-5)
            th_coeff = self.bada_reader.get_param("PFM/TFM/MCMB/flat_rating", "b")
            th_coeff = [float(elem) for elem in th_coeff]

            d_coef = []
            for m_coef in [th_coeff[i : i + 6] for i in range(0, len(th_coeff), 6)]:
                d_coef.append(self.poly(m_coef, mach))

            throttle = self.poly(d_coef, dP)

        # flat-rated area: eq (3.3-3)
        ct_coeff = self.bada_reader.get_param("PFM/TFM/CT", "a")
        ct_coeff = [float(elem) for elem in ct_coeff]

        d_coef = []
        for m_coef in [ct_coeff[i : i + 6] for i in range(0, len(ct_coeff), 6)]:
            d_coef.append(self.poly(m_coef, mach))

        ct = self.poly(d_coef, throttle)

        return ct, throttle

    def get_thrust(self, ct, alt):

        dP = self.get_dp(alt)

        mref = float(self.bada_reader.get_param("PFM/MREF"))
        w_mref = mref * self.aero.g0

        T = dP * w_mref * ct

        return T
