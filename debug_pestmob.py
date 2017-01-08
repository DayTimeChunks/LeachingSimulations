from leach_hydrology import *
from pestmob import *


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

"""
water_data1 = leachsim(dtGA=1,
                       ov=0.30,
                       ovSat=ovSat_crop,
                       kSat=kSat_crop,
                       psi=psi_crop,
                       soil_height=30)

"""
water_data2 = leachsim(dtGA=1,
                       ov=ov_2,
                       ovSat=ovSat_crop,
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

# Hydro data
percol_data2 = stackdata3(cum_time_30min,
                          cum_leach_135mmh, cum_leach_55mmh, cum_leach_30mmh)

runoff_data2 = stackdata3(cum_time_30min,
                          cum_roff_135mmh, cum_roff_55mmh, cum_roff_30mmh)

infil_data2 = stackdata3(cum_time_30min,
                         infil_135mmh, infil_55mmh, infil_30mmh)


"""
hydroplot(runoff_data2,
          "Ponding at 135mmh", "Ponding at 55mmh", "Ponding at 30mmh",
          leach_high_6min,
          leach_med_12min, leach_med_30min,
          leach_low_30min)


"""
pb_crop = 0.99/10**3  # bulk density (g/cm^3) -> g/mm^3
porosity_crop = 0.61  # Crop soil
runoff_vel = 0.000001  # mm/min

# Soil characteristics (defined above)
fom_crop_sterile = 3.87/100.0
fom_crop_untreat = 5.51/100.0
foc_crop_sterile = 0.58*fom_crop_sterile
foc_crop_untreat = 0.58*fom_crop_untreat

# Soil characteristics (OC Black & Walkley)
foc_crop_untreat2 = 2.04/100
foc_crop_sterile2 = 0.70*foc_crop_untreat2

# Pesticide Koc
Koc_mexyl = [163.0, 50.0, 30]  # [(a) , (b), (c)] [ml/g]
Koc_mexyl = np.array(Koc_mexyl)/10**3  # [mm3/g]

# Kd (a) - NPIC @ http://npic.orst.edu/ingred/ppdmove.htm
Kd_mexylA_crop_sterile = Koc_mexyl[0]*foc_crop_sterile  # ml/g
Kd_mexylA_crop_untreat = Koc_mexyl[0]*foc_crop_untreat

# Kd (b) - PAN @ http://www.pesticideinfo.org/
Kd_mexylB_crop_sterile = Koc_mexyl[1]*foc_crop_sterile  # ml/g
Kd_mexylB_crop_untreat = Koc_mexyl[1]*foc_crop_untreat

# Kd (c) - https://toxnet.nlm.nih.gov/cgi-bin/sis/search/a?dbs+hsdb:@term+@DOCNO+7061
Kd_mexylC_crop_sterile = Koc_mexyl[2]*foc_crop_sterile  # ml/g
Kd_mexylC_crop_untreat = Koc_mexyl[2]*foc_crop_untreat

Kd_mexyl = [Kd_mexylA_crop_sterile, Kd_mexylA_crop_untreat,
            Kd_mexylB_crop_sterile, Kd_mexylB_crop_untreat,
            Kd_mexylC_crop_sterile, Kd_mexylC_crop_untreat]

# Initial mass
mx_ini_sterile_list = [1496.75, 1440.72, 1047.95, 1462.08,  # 0 days
                       1127.52, 1267.11, 994.09, 1050.48]  # 10 days

mx_ini_untreat_list = [1222.86, 1211.28, 751.13, 979.82,  # 0 days
                       1006.54, 1025.43, 938.23, 830.68]  # 10 days

mx_ini_sterile = sum(mx_ini_sterile_list) / len(mx_ini_sterile_list)
mx_ini_untreated = sum(mx_ini_untreat_list) / len(mx_ini_untreat_list)

mx_ini = [mx_ini_sterile, mx_ini_untreated]

# Observed Metalaxyl Output - Crop
# high, med-12, med-30, low
mx_obs_sterile_crop = np.array([(8.35 + 37.57) / 2.0,
                                (290.27) / 1.0, (93.29 + 82.15) / 2.0,
                                (285.29) / 1.0])

mx_obs_untreat_crop = np.array([(175.44 + 40.03) / 2.0,
                                (272.48 + 168.51) / 2.0, (35.12 + 146.10) / 2.0,
                                (86.01 + 76.49) / 2.0])

mx_sol = [mx_obs_sterile_crop, mx_obs_untreat_crop]


cum_mx_crop = pest_test(Kd_mexyl,
                        pb_crop,
                        ovSat_crop,
                        percol_data2, runoff_data2, infil_data2,
                        area, soil_height,
                        mx_ini,
                        mx_sol,
                        d, runoff_vel, dtGA=1)





# pass a data set with evolution of mass in overland flow for each intensity
cum_time_30min = cum_mx_crop[:, 0]

# Cumulative leachate sterilized
cum_mass_leach_st_135mmh = cum_mx_crop[:, 1]
cum_mass_leach_st_55mmh = cum_mx_crop[:, 2]
cum_mass_leach_st_30mmh = cum_mx_crop[:, 3]

# Cumulative leachate untreated
cum_mass_leach_un_135mmh = cum_mx_crop[:, 4]
cum_mass_leach_un_55mmh = cum_mx_crop[:, 5]
cum_mass_leach_un_30mmh = cum_mx_crop[:, 6]

# Leachate sterilized
mass_leach_st_135mmh = cum_mx_crop[:, 7]
mass_leach_st_55mmh = cum_mx_crop[:, 8]
mass_leach_st_30mmh = cum_mx_crop[:, 9]

# Leachate untreated
mass_leach_un_135mmh = cum_mx_crop[:, 10]
mass_leach_un_55mmh = cum_mx_crop[:, 11]
mass_leach_un_30mmh = cum_mx_crop[:, 12]


# Group each compartment for graphing
cum_leach_mx_crop = stackdata6(cum_time_30min,
                               cum_mass_leach_st_135mmh, cum_mass_leach_st_55mmh, cum_mass_leach_st_30mmh,
                               cum_mass_leach_un_135mmh, cum_mass_leach_un_55mmh, cum_mass_leach_un_30mmh)

# NO RUN OFF DATA
"""
mass_percol2 = stackdata6(cum_time_30min,
        cum_mass_leach_st_135mmh, cum_mass_leach_st_55mmh, cum_mass_leach_st_30mmh,
        cum_mass_leach_un_135mmh, cum_mass_leach_un_55mmh, cum_mass_leach_un_30mmh)

mass_runoff_st_135mmh = cum_mx_crop[:, 13]
mass_runoff_st_55mmh = cum_mx_crop[:, 14]
mass_runoff_st_30mmh = cum_mx_crop[:, 15]

mass_runoff_un_135mmh = cum_mx_crop[:, 16]
mass_runoff_un_55mmh = cum_mx_crop[:, 17]
mass_runoff_un_30mmh = cum_mx_crop[:, 18]

data = stackdata6(cum_time_30min,
                  mass_runoff_st_135mmh, mass_runoff_st_55mmh, mass_runoff_st_30mmh,
                  mass_runoff_un_135mmh, mass_runoff_un_55mmh, mass_runoff_un_30mmh)

# print(data[:][:, 0])  # time axis
print(data[:][:, 1])  # mass_runoff_st_135mmh
print(len(data[0]))
"""


#pestiplot(cum_leach_mx_crop, mx_obs_sterile_crop, mx_obs_untreat_crop,
#          'Cumulative Metalaxyl - Annual Crop Soil')

