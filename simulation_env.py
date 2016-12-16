# -*- coding: cp1252 -*-

from greenampt import *
from pestievent import *
from metalevent import *
import numpy as np
from hydroplots import *
from leach_hydrology import *

""" Contaminant characteristics & conditions """
Kd_Smeto = 4.625  # S-metolachor; Koc: 185 ml/g x foc (mean): 0.025 == 4.625 mL/g
# Copper log(Kd) range: 0.1 - 7.0, max-mean = 5.5
Kd_copper_min = 10**2.7  # [mL/g] = [L/Kg] (Allison and Allison, 2005 - EPA/600/R-05/074)
Kd_copper_max = 10**3.6  # [mL/g] = [L/Kg] (Allison and Allison, 2005 - EPA/600/R-05/074)

# Zinc log(Kd) range: 1.5 - 6.9
Kd_zinc_min = 10**2.9  # [mL/g] = [L/Kg] (Allison and Allison, 2005 - EPA/600/R-05/074)
Kd_zinc_max = 10**3.9  # [mL/g] = [L/Kg] (Allison and Allison, 2005 - EPA/600/R-05/074)

# Mass leached at 6min, 12min, 30min, 30min @ 135, 55, 55, 30mm/h
sol_sterile = np.array([11.29, 11.63, 306.80, 21.08])  # [ug]
sol_untreat = np.array([0.5, 1.405, 37.0, 1])

""" Soil Characteristics """
pb_vine = 1.04  # bulk density (g/cm^3) # Rouffach (Martine Trautmann, sampled pre-event)
pb_crop = 0.99  # bulk density (g/cm^3) # Alteck (Martine Trautmann, sampled pre-event)
porosity_crop = 0.61  # Crop soil
porosity_vine = 0.55  # Vinyard soil
kSat_vine = 0.293  # cm/min (17.6 cm/h - Crop Soil) - Rouffach (Martine Trautmann, sampled pre-event)
kSat_crop = 0.225  # cm/min (13.5 cm/h - Crop Soil) - Alteck (Martine Trautmann, sampled pre-event)
ov_1 = 0.25   # Initial water content m3. m-3
ov_2 = 0.30   # Initial water content m3. m-3
ovSat_crop = 0.45  # Saturated water content (assumed)
ovSat_vine = 0.45  # Saturated water content (assumed)
psi_crop = 110  # soil suction Alteck
# psi_vine = xxxx  # soil suction Rouffach

""" Microcosm """
d = (1.493 * 2)  # Diameter of falcon tube (cm)
area = ((d / 2) ** 2) * 3.1416  # (cmÂ²)
zl = soil_height = 3.0  # Mixing layer depth in cm

""" Soil modelled """
pb = pb_crop
ov = ov_1
ovSat = ovSat_crop
psi = psi_crop

""" Leached calculation """
water_data = leachsim()

cum_time_30min = water_data[:, 0]
cum_inf_135mmh = water_data[:, 4]
cum_inf_55mmh = water_data[:, 5]
cum_inf_30mmh = water_data[:, 6]
cum_leach_135mmh = water_data[:, 7]
cum_leach_55mmh = water_data[:, 8]
cum_leach_30mmh = water_data[:, 9]

data = stackdata6(cum_time_30min,
                 cum_inf_135mmh, cum_inf_55mmh, cum_inf_30mmh,
                 cum_leach_135mmh, cum_leach_55mmh, cum_leach_30mmh)

hydroplot(data)

""" Copper simulation 1 - Retardation """

""" Copper Lab Data"""
#  Mass initial (untreated)
mass_ini_sterile = (1627 + 1107) / float(2)  # all intensities:{0d, 10d)
mass_ini_untreated = (1184 + 1177) / float(2)  # all intesities:{0d, 10d)

cum_copper_Alteck = metal_ret(Kd_copper_min, Kd_copper_max,
                              pb,
                              ovSat,
                              water_data,
                              area, soil_height,
                              mass_ini_sterile, mass_ini_untreated)

""" Observed Cupper """
cu_sol_sterile = np.array([11.29, 11.63, 306.80, 21.08])
cu_sol_untreat = np.array([0.5, 1.405, 37.0, 1])

metalplot(cum_copper_Alteck, cu_sol_sterile, cu_sol_untreat, 'Cumulative Cu')

""" Zinc simulation 1 - Retardation """
mass_ini_sterile = (3106 + 2594) / float(2)  # all intensities:{0d, 10d)
mass_ini_untreated = (2636 + 2586) / float(2)
cum_zinc_Alteck = metal_ret(Kd_zinc_min, Kd_zinc_max,
                            pb,
                            ovSat,
                            water_data,
                            area, soil_height,
                            mass_ini_sterile, mass_ini_untreated)

""" Observed Cupper Zinc"""
zn_sol_sterile = np.array([17, 17, 404.25, 30.13])
zn_sol_untreat = np.array([2.19, 2, 36.85, 3.09])

metalplot(cum_zinc_Alteck, zn_sol_sterile, zn_sol_untreat, 'Cumulative Zn')
