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
Kd_zinc_min = 10**3.1  # [mL/g] = [L/Kg] (Allison and Allison, 2005 - EPA/600/R-05/074)
Kd_zinc_max = 10**5.1  # [mL/g] = [L/Kg] (Allison and Allison, 2005 - EPA/600/R-05/074)

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

hydroplot2(data)

""" Copper simulation 1 - Retardation """
cum_copper_Alteck = copper_ret(Kd_copper_min, Kd_copper_max,
                               pb,
                               ovSat,
                               water_data,
                               area, soil_height)

copperplot(cum_copper_Alteck)

