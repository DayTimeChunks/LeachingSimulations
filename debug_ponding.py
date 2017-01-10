
from leach_hydrology import *
from pestmob import *
from mixinglayer import *

intensities = [2.25, 0.92, 0.5]  # mm/min = [0.225, 0.092, 0.05] cm/min
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

time_size_135mmh = water_data2[:, 19]
time_size_55mmh = water_data2[:, 20]
time_size_30mmh = water_data2[:, 21]

# Hydro data
percol_data2 = stackdata3(cum_time_30min,
                          cum_leach_135mmh, cum_leach_55mmh, cum_leach_30mmh)

runoff_data2 = stackdata3(cum_time_30min,
                          cum_roff_135mmh, cum_roff_55mmh, cum_roff_30mmh)

infil_data2 = stackdata3(cum_time_30min,
                         infil_135mmh, infil_55mmh, infil_30mmh)

time_sizes2 = [time_size_135mmh, time_size_55mmh, time_size_30mmh]

# print(len(time_size_135mmh), len(cum_leach_135mmh))
#####################################################################
# END OF HYDRO  -> START Contaminants
#####################################################################

# Soil characteristics
pb_crop = 0.99/10**3  # bulk density (g/cm^3) -> g/mm^3
porosity_crop = 0.61  # Crop soil
soil_height2 = 23

# Assumed (used to calculate Reynolds number)
runoff_vel = 100  # mm/min

# Fraction organic matter and carbon (Scenario 1)
fom_crop_sterile = 3.87/100.0
fom_crop_untreat = 5.51/100.0
foc_crop_sterile = 0.58*fom_crop_sterile
foc_crop_untreat = 0.58*fom_crop_untreat

# Soil characteristics (OC Black & Walkley - Scenario 2)
foc_crop_untreat2 = 2.04/100
foc_crop_sterile2 = 0.70*foc_crop_untreat2

# Pesticide Koc
Koc_mexyl = [163.0, 50.0, 30]  # [(a) , (b), (c)] [ml/g]
Koc_mexyl = np.array(Koc_mexyl)*10**3  # [mm3/g]

# Kd (a) - NPIC @ http://npic.orst.edu/ingred/ppdmove.htm
Kd_mexylA_crop_sterile = Koc_mexyl[0]*foc_crop_sterile
Kd_mexylA_crop_untreat = Koc_mexyl[0]*foc_crop_untreat

# Kd (b) - PAN @ http://www.pesticideinfo.org/
Kd_mexylB_crop_sterile = Koc_mexyl[1]*foc_crop_sterile
Kd_mexylB_crop_untreat = Koc_mexyl[1]*foc_crop_untreat

# Kd (c) - https://toxnet.nlm.nih.gov/cgi-bin/sis/search/a?dbs+hsdb:@term+@DOCNO+7061
Kd_mexylC_crop_sterile = Koc_mexyl[2]*foc_crop_sterile
Kd_mexylC_crop_untreat = Koc_mexyl[2]*foc_crop_untreat

Kd_mexyl = [Kd_mexylA_crop_sterile, Kd_mexylA_crop_untreat,
            Kd_mexylB_crop_sterile, Kd_mexylB_crop_untreat,
            Kd_mexylC_crop_sterile, Kd_mexylC_crop_untreat]


# Initial mass - 2nd pulse
mx_ini_sterile_list = [1496.75, 1440.72, 1047.95, 1462.08,  # 0 days
                       1127.52, 1267.11, 994.09, 1050.48]  # 10 days

mx_ini_untreat_list = [1222.86, 1211.28, 751.13, 979.82,  # 0 days
                       1006.54, 1025.43, 938.23, 830.68]  # 10 days

mx_ini_sterile = sum(mx_ini_sterile_list) / len(mx_ini_sterile_list)
mx_ini_untreated = sum(mx_ini_untreat_list) / len(mx_ini_untreat_list)

mx_ini = [mx_ini_sterile, mx_ini_untreated]

# Observed Metalaxyl Leachate - Crop
# high, med-12, med-30, low
mx_obs_sterile_crop = np.array([(8.35 + 37.57) / 2.0,
                                (290.27) / 1.0, (93.29 + 82.15) / 2.0,
                                (285.29) / 1.0])

mx_obs_untreat_crop = np.array([(175.44 + 40.03) / 2.0,
                                (272.48 + 168.51) / 2.0, (35.12 + 146.10) / 2.0,
                                (86.01 + 76.49) / 2.0])

mx_sol_leach = [mx_obs_sterile_crop, mx_obs_untreat_crop]

# Observed Metalaxyl in Ponding - Crop
# high, med-12, med-30, low
mx_obs_sterile_crop_roff = np.array([(5.74 + 4.07) / 2.0,
                                     (9.16) / 1.0, (4.34 + 14.03) / 2.0,
                                     (12.42) / 1.0])

mx_obs_untreat_crop_roff = np.array([(4.72 + 3.17) / 2.0,
                                     (1.84) / 1.0, (8.90 + 0.14) / 2.0,
                                     (5.82 + 9.65) / 1.0])

mx_sol_pond = [mx_obs_sterile_crop_roff, mx_obs_untreat_crop_roff]

# Any length unit input must be: "mm"
cum_mx_crop2 = pest_test(Kd_mexyl,
                         intensities,
                         pb_crop,
                         ovSat_crop,
                         percol_data2, runoff_data2, infil_data2, time_sizes2,
                         area, soil_height2,
                         mx_ini,
                         mx_sol_leach, mx_sol_pond,
                         d, runoffvelocity=100,
                         KFILM=True)

# Time axis
cum_time_30min = cum_mx_crop2[:, 0]

# Ponding sterilized
mass_runoff_st_135mmh = cum_mx_crop2[:, 13]
mass_runoff_st_55mmh = cum_mx_crop2[:, 14]
mass_runoff_st_30mmh = cum_mx_crop2[:, 15]

# Ponding untreated
mass_runoff_un_135mmh = cum_mx_crop2[:, 16]
mass_runoff_un_55mmh = cum_mx_crop2[:, 17]
mass_runoff_un_30mmh = cum_mx_crop2[:, 18]

mass_overflow2 = stackdata6(cum_time_30min,
        mass_runoff_st_135mmh, mass_runoff_st_55mmh, mass_runoff_st_30mmh,
        mass_runoff_un_135mmh, mass_runoff_un_55mmh, mass_runoff_un_30mmh)

# Organized individually:
# [sterile, untreat, sterile_aged, untreat_aged]
# all at 6 min, high inetnesity
roff2_mx_high_6min = np.array([5.7, 4.7, 4.1, 3.2])

# all at 12 min, med intensity
roff2_mx_med_12min = np.array([0., 1.8, 9.2, 0.])

# all at 30min, med intensity
roff2_mx_med_30min = np.array([4.3, 8.9, 14.0, 0.1])

# all at 30min, low intensity
roff2_mx_low_30min = np.array([0., 5.8, 12.4, 9.6])

pestiplot_all(mass_overflow2,
              roff2_mx_high_6min,
              roff2_mx_med_12min, roff2_mx_med_30min,
              roff2_mx_low_30min,
              'Mass in Ponding Metalaxyl - 2nd Pulse Crop Soil',
              'Metalaxyl')