
from hydroplots import *
from leach_hydrology import *

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



water_data2 = leachsim(dtGA=1,
                       ov=0.30,
                       ovSat=ovSat_crop,
                       kSat=kSat_crop2,
                       psi=psi_crop,
                       soil_height=30)
"""

dtGA = 1  # timestep size, 1 minute
intensity = 2.25  # mm/min
StormD = 30  # Storm duration in min

wfs = psi_crop * (ovSat - ov_1)

# Initialization
cum_time = 0.0
cum_runoff = 0.0
cum_precip = 0.0
cum_leach = 0.0
cum_inf = 0.0

cum_time_dt = []
cum_precip_dt = []
cum_runoff_dt = []
cum_inf_dt = []  # F2
cum_leach_dt = []

runoff_dt = []
infil_dt = []
leach_dt = []

# Start of rainfall event X
while cum_time < StormD:
    cum_time += dtGA
    cum_precip += intensity*dtGA

    green_output = greenampt(ksat=kSat_crop2,
                             wfs=wfs,
                             dtga=dtGA,
                             intensity=intensity,
                             ovsat=ovSat,
                             ov=ov_1,
                             zl=zl,
                             cum_infilt=cum_inf,
                             cum_time=cum_time)

    cum_inf = green_output[0]
    inf = green_output[1]
    runoff = green_output[2]
    leach = green_output[3]

    if cum_time % 10 == 0:
        print("########### time step: ", cum_time)
        print(cum_inf)

"""
