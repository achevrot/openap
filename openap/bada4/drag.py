from termios import CR0
from . import BadaReader
import math
from typing import Any, Tuple

from pitot.wrapper import default_units
from pitot import Q_, ureg
import pitot.isa as isa
import pitot.aero as aero

import numpy as np

from functools import reduce

class Drag(object):
    def __init__(self, aircraft_type, eng=None, **kwargs):

        self.bada_reader = BadaReader(aircraft_type)


    @default_units(s="m^2")
    def get_drag(self, mach, dP, s, cd):
        result = 1/2 * dP * isa.P_0 * isa.GAMMA * s * mach**2 * cd
        return result.to('N')
    
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
    
    @default_units(mass='kg', tas='kts', alt="m", s="m^2", path_angle='rad')
    def get_cl(self, mass, tas, alt, dP, s, path_angle=0):
        dP = self.get_dp(alt)
        s = self.bada_reader.get_param("AFCM/S")
        mach = aero.tas2mach(tas, alt)
        cl = 2 * mass * isa.G_0 / (dP * isa.P_0 * isa.GAMMA * s * mach**2 * math.cos(path_angle))
        return cl
    
    @default_units(alt="m")
    def clean(self, mass, tas, alt, path_angle=0):
        
        def calc_c(acc, c):
            return acc + c[0]/((1-mach**2)**c[1])

        
        mach = aero.tas2mach(tas, alt)
        dP = self.get_dp(alt)
        s = self.bada_reader.get_param("AFCM/S")
        cl = self.get_cl(mass, tas, alt, dP, s, path_angle)
        mach_max = self.bada_reader.get_param("AFCM/Configuration/LGUP/DPM_clean/M_max")
        
        if mach > mach_max:
            result = 2
        else:
            d_c = self.bada_reader.get_param("AFCM/Configuration/LGUP/DPM_clean/CD_clean")
            l0 = zip(d_c[0:5], np.concatenate(([0], np.arange(0.5, 2.5, 0.5))))
            c0 = reduce(calc_c, l0, 0)
            l2 = zip(d_c[5:10], np.concatenate(([0], np.arange(1.5, 6.5, 1.5))))
            c2 = reduce(calc_c, l2, 0)
            l6 = zip(d_c[10:15], np.concatenate(([0], np.arange(7,9.0,0.5))))
            c6 = reduce(calc_c, l6, 0)
            scalar = self.bada_reader.get_param("AFCM/Configuration/LGUP/DPM_clean/scalar")
            cd = scalar * (c0 + (c2 * cl**2) + (c6 * cl**6))
        
        return self.get_drag(mach, dP, s, cd)
    
    @default_units(alt="m")
    def nonclean(self, mass, tas, alt, flap_angle, path_angle=0, landing_gear=False):
        lg = 'LGDN' if landing_gear else 'LGUP'
        # TODO : Check different high-lift configs
        d1, d2, d3 = self.bada_reader.get_param(f"AFCM/Configuration/{lg}/DPM_nonclean/CD_nonclean")
        cl = self.get_cl(mass, tas, alt, path_angle)
        cd = d1 + d2 * cl + d3 * cl**2
                
        return cd
    
    