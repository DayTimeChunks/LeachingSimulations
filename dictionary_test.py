from pestmob import *
from leach_hydrology import *


# Alteck (Martine Trautmann, sampled pre-event)
kSat_crop = 2.24  # mm/min = 0.224 cm/min (13.45 cm/h - Crop Soil)
kSat_crop2 = kSat_crop/100
ovSat = 0.45  # Saturated water content (assumed)
ov_1 = 0.25   # Initial water content m3. m-3
ov_2 = 0.40   # Initial water content m3. m-3
ovSat_crop = 0.45  # Saturated water content (assumed)
psi_crop = 1100  # soil suction Alteck mm
#  (Lefrancq, 2014: 61.7 cm , p.160; 110 cm, p.189)

d = (14.93 * 2)  # Diameter of falcon tube (mm)
area = ((d / 2) ** 2) * 3.1416  # (mm2)
soil_height2 = 23  # mm

# Assumed (used to calculate Reynolds number)
runoff_vel = 10.0  # mm/min

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

#  Dictionary contains:
#  Scenario: (initial_mas, leached_mass_observed, ponded_mass_obs)
pest_dict_S = {
    'highS0d': (1496.75, 8.35, 5.7),
    'highS1d': (1127.52, 37.57, 4.1),
    'med12S0d': (1440.72, 290.3, 0.),
    'med12S1d': (1267.11, 0., 9.2),
    'med30S0d': (1047.95, 93.3, 4.3),
    'med30S1d': (994.09, 82.2, 14.0),
    'lowS0d': (1462.08, 285.3, 0.),
    'lowS1d': (1050.48, 0., 12.4)
}

pest_dict_L = {
    'highL0d': (1222.86, 175.44, 4.7),
    'highL1d': (1006.54, 40.03, 3.2),
    'med12L0d': (1211.28, 272.5, 1.8),
    'med12L1d': (1025.43, 168.5, 0.),
    'med30L0d': (751.13, 35.1, 8.9),
    'med30L1d': (938.23, 146.1, 0.1),
    'lowL0d': (979.82, 86.0, 5.8),
    'lowL1d': (830.68, 76.5, 9.6)
}

pb_crop = 0.99/10**3  # bulk density (g/cm^3) -> g/mm^3
porosity_crop = 0.61  # Crop soil
runoff_vel = 1.0  # mm/min

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
Koc_mexyl = np.array(Koc_mexyl)*10**3 # [mm3/g]

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


sterile_sim = pest_test2(Kd_mexyl,
                         mx_dict_S,
                         pb_crop,
                         ovSat_crop,
                         percol_data2, runoff_data2, infil_data2, time_sizes2,
                         area, soil_height2,
                         d, runoff_vel)

print(sterile_sim[:, 1])


# Time axis
cum_time_30min = sterile_sim[:, 0]

# Cumulative leachate sterilized
a_high_0d_cum_mass_out_dt = sterile_sim[:, 1]
b_high_1d_cum_mass_out_dt = sterile_sim[:, 2]

c_med12_0d_cum_mass_out_dt = sterile_sim[:, 3]
d_med12_1d_cum_mass_out_dt = sterile_sim[:, 4]

e_med30_0d_cum_mass_out_dt = sterile_sim[:, 5]
f_med30_1d_cum_mass_out_dt = sterile_sim[:, 6]

g_low_0d_cum_mass_out_dt = sterile_sim[:, 7]
h_h_low_1d_cum_mass_out_dt = sterile_sim[:, 8]

# Ponded mass
a_high_0d_overmass_dt = sterile_sim[:, 9]
b_high_1d_overmass_dt = sterile_sim[:, 10]
c_med12_0d_overmass_dt = sterile_sim[:, 11]
d_med12_1d_overmass_dt = sterile_sim[:, 12]

e_med30_0d_overmass_dt = sterile_sim[:, 13]
f_med30_1d_overmass_dt = sterile_sim[:, 14]
g_low_0d_overmass_dt = sterile_sim[:, 15]
h_h_low_1d_overmass_dt = sterile_sim[:, 16]

mass_percol2 = stackdata8(cum_time_30min,
                          a_high_0d_cum_mass_out_dt, b_high_1d_cum_mass_out_dt,
                          c_med12_0d_cum_mass_out_dt, d_med12_1d_cum_mass_out_dt,
                          e_med30_0d_cum_mass_out_dt, f_med30_1d_cum_mass_out_dt,
                          g_low_0d_cum_mass_out_dt, h_h_low_1d_cum_mass_out_dt)






