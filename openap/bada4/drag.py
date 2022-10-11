from functools import reduce
from math import cos, exp, pow
import numpy as np
import pitot.aero as aero
import pitot.isa as isa
from pitot import Q_
from pitot.wrapper import default_units

from . import BadaReader


class Drag(object):
    
    def __init__(self, aircraft_type, eng=None, **kwargs):

        self.bada_reader = BadaReader.getAircraft(aircraft_type)
        self.dp_dict = {}

    def get_drag(self, mach, dP, s, cd):
        result = 1/2 * dP * isa.P_0.m * isa.GAMMA.m * s * mach*mach * cd
        return Q_(result, 'N')
    
    def get_dp(self, alt=Q_(0, "m / s")):
        
        if f"{alt}" in self.dp_dict:
            return self.dp_dict[f"{alt}"]

        t_0 = isa.temperature(Q_(0, "m")).m
        p_0 = isa.pressure(Q_(0, "m")).m
        t_trop = isa.temperature(isa.H_TROP).m

        if alt < isa.H_TROP.m:
            # bada (2.2-18)
            dP = (
                p_0
                * pow(
                    (isa.temperature(Q_(alt, "m")).m / t_0), (-isa.G_0.m / (isa.BETA_T.m * isa.R.m))
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
            
        self.dp_dict["alt"] = dP
        return dP
    
    def get_cl(self, mass, tas, alt, dP, s, path_angle=0):

        s = self.bada_reader.get_param("AFCM/S").m
        mach = aero.tas2mach(Q_(tas, 'kts'), Q_(alt, 'm')).m
        cl = 2 * mass * isa.G_0.m / (dP * isa.P_0.m * isa.GAMMA.m * s * mach*mach * cos(path_angle))
        return cl
    
    @default_units(mass='kg', tas='kts', alt="m", path_angle='degree')
    def clean(self, mass, tas, alt, path_angle=Q_(0, 'degree')):
        return self.get_clean(mass.m, tas.m, alt.m, path_angle.m)

    def get_clean(self, mass, tas, alt, path_angle):
    
        def calc_c(acc, c):
            return acc + c[0]/(pow((1-mach*mach),c[1]))

        mach = aero.tas2mach(Q_(tas, 'kts'), Q_(alt, 'm')).m
        dP = self.get_dp(alt)
        s = self.bada_reader.get_param("AFCM/S").m
        cl = self.get_cl(mass, tas, alt, dP, s, path_angle)
        mach_max = self.bada_reader.get_param("AFCM/Configuration/LGUP/DPM_clean/M_max").m
        
        if mach > mach_max:
            result = 2
        else:
            d_c = self.bada_reader.get_param("AFCM/Configuration/LGUP/DPM_clean/CD_clean").m
            l0 = zip(d_c[0:5], np.concatenate(([0], np.arange(0.5, 2.5, 0.5))))
            c0 = reduce(calc_c, l0, 0)
            l2 = zip(d_c[5:10], np.concatenate(([0], np.arange(1.5, 6.5, 1.5))))
            c2 = reduce(calc_c, l2, 0)
            l6 = zip(d_c[10:15], np.concatenate(([0], np.arange(7,9.0,0.5))))
            c6 = reduce(calc_c, l6, 0)
            scalar = self.bada_reader.get_param("AFCM/Configuration/LGUP/DPM_clean/scalar").m
            cd = scalar * (c0 + (c2 * cl*cl) + (c6 * pow(cl,6)))
        
        return self.get_drag(mach, dP, s, cd)
    
    @default_units(alt="m")
    def nonclean(self, mass, tas, alt, flap_angle, path_angle=Q_(0, 'degree'), landing_gear=False):
        lg = 'LGDN' if landing_gear else 'LGUP'
        # TODO : Check different high-lift configs
        d1, d2, d3 = self.bada_reader.get_param(f"AFCM/Configuration/{lg}/DPM_nonclean/CD_nonclean")
        cl = self.get_cl(mass, tas, alt, path_angle)
        cd = d1 + d2 * cl + d3 * cl*cl
                
        return cd
    
    