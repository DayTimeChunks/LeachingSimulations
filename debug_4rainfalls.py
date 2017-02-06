from leach_hydrology import *
from pestmob import *

"""
This debug script will pass initial
concentrations to the pesticide mobilization function (pest_test2).
This function should take a list of initial concentrations.
"""

d = (14.93 * 2)  # Diameter of falcon tube (mm)
area = ((d / 2) ** 2) * 3.1416  # (mm2)
zl = soil_height = 30  # Mixing layer depth in mm

''' Hydrological controlling parameters - all in mm/h'''
# Alteck (Martine Trautmann, sampled pre-event)
kSat_crop = 2.24  # mm/min = 0.224 cm/min (13.45 cm/h - Crop Soil)
kSat_crop2 = kSat_crop/100
ovSat = 0.45  # Saturated water content (assumed)
ov_1 = 0.25   # Initial water content m3. m-3
ov_2 = 0.40   # Initial water content m3. m-3
ovSat_crop = 0.45  # Saturated water content (assumed)
psi_crop = 1100  # soil suction Alteck mm
#  (Lefrancq, 2014: 61.7 cm , p.160; 110 cm, p.189)

water_data2 = leachsim(dtGA=1,
                       ov=ov_2,
                       kSat=kSat_crop2,
                       psi=psi_crop,
                       soil_height=23)

# Convert volumes from mm3 to mL in hydroplot function
cum_time_30min = water_data2[:, 0]
cum_inf_135mmh = water_data2[:, 4]
cum_inf_55mmh = water_data2[:, 5]
cum_inf_30mmh = water_data2[:, 6]
cum_leach_135mmh = water_data2[:, 7]
cum_leach_55mmh = water_data2[:, 8]
cum_leach_30mmh = water_data2[:, 9]

roff_135mmh = water_data2[:, 10]
roff_55mmh = water_data2[:, 11]
roff_30mmh = water_data2[:, 12]

cum_roff_135mmh = water_data2[:, 13]
cum_roff_55mmh = water_data2[:, 14]
cum_roff_30mmh = water_data2[:, 15]

infil_135mmh = water_data2[:, 16]
infil_55mmh = water_data2[:, 17]
infil_30mmh = water_data2[:, 18]

# Hydro data
percol_data2 = stackdata3(cum_time_30min,
                          cum_leach_135mmh, cum_leach_55mmh, cum_leach_30mmh)

runoff_data2 = stackdata3(cum_time_30min,
                          cum_roff_135mmh, cum_roff_55mmh, cum_roff_30mmh)

infil_data2 = stackdata3(cum_time_30min,
                         infil_135mmh, infil_55mmh, infil_30mmh)

""" Observed Percolation Annual Crop """

# [sterile, untreat, sterile, untreat]
# all at 6 min, high inetnesity
# observed data in mL
leach_high_6min = np.array([14.192, 8.245, 2.410, 5.469])

# all at 12 min, med intensity
leach_med_12min = np.array([18.672, 19.0, 0.830, 11.407])

# all at 30min, med intensity
leach_med_30min = np.array([12.697, 2.473, 3.52, 20.291])

# all at 30min, low intensity
leach_low_30min = np.array([29.656, 9.375, 0.409, 3.385])




"""
hydroplot(runoff_data2,
          "Ponding at 135mmh", "Ponding at 55mmh", "Ponding at 30mmh",
          leach_high_6min,
          leach_med_12min, leach_med_30min,
          leach_low_30min)


"""


