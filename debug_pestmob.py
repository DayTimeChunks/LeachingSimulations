from leach_hydrology import *
from pestmob import *

d = (14.93 * 2)  # Diameter of falcon tube (mm)
area = ((d / 2) ** 2) * 3.1416  # (mm2)
zl = soil_height = 30  # Mixing layer depth in mm

''' Hydrological controlling parameters - all in mm/h'''
# Alteck (Martine Trautmann, sampled pre-event)
pb_crop = 0.99  # bulk density (g/cm^3)
porosity_crop = 0.61  # Crop soil
kSat_crop = 2.24  # mm/min = 0.224 cm/min (13.45 cm/h - Crop Soil)
kSat_crop2 = kSat_crop
ovSat = 0.45  # Saturated water content (assumed)
ov_1 = 0.25   # Initial water content m3. m-3
ov_2 = 0.20   # Initial water content m3. m-3
ovSat_crop = 0.45  # Saturated water content (assumed)
psi_crop = 1100  # soil suction Alteck cm
#  (Lefrancq, 2014: 61.7 cm , p.160; 110 cm, p.189)

water_data1 = leachsim(dtGA=1,
                       ov=0.30,
                       ovSat=ovSat_crop,
                       kSat=kSat_crop,
                       psi=psi_crop,
                       soil_height=30)




