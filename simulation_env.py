# -*- coding: cp1252 -*-

from greenampt import *
from pestievent import *
from pestmob import *
import numpy as np
from hydroplots import *
from leach_hydrology import *

""" Contaminant characteristics & conditions """
Koc_smeto = 185.0  # [ml/g]
Koc_mexyl = 163.0  # [ml/g]

fom_vine_sterile = 3.53/100.0
fom_vine_untreat = 2.93/100.0
fom_crop_sterile = 3.87/100.0
fom_crop_untreat = 5.51/100.0

foc_vine_sterile = 0.58*fom_vine_sterile
foc_vine_untreat = 0.58*fom_vine_untreat
foc_crop_sterile = 0.58*fom_crop_sterile
foc_crop_untreat = 0.58*fom_crop_untreat

Kd_smeto_crop_sterile = Koc_smeto*foc_crop_sterile  # ml/g
Kd_smeto_crop_untreat = Koc_smeto*foc_crop_untreat  # ml/g
Kd_smeto_vine_sterile = Koc_smeto*foc_vine_sterile
Kd_smeto_vine_untreat = Koc_smeto*foc_vine_untreat

Kd_mexyl_crop_sterile = Koc_mexyl*foc_crop_sterile  # ml/g
Kd_mexyl_crop_untreat = Koc_mexyl*foc_crop_untreat
Kd_mexyl_vine_sterile = Koc_mexyl*foc_vine_sterile  # ml/g
Kd_mexyl_vine_untreat = Koc_mexyl*foc_vine_untreat  # ml/g

# Zinc log(Kd) range: 1.5 - 6.9
Kd_zinc_sterile = 10**2.9  # [mL/g] = [L/Kg] (Allison and Allison, 2005 - EPA/600/R-05/074)
Kd_zinc_untreat = 10**3.9  # [mL/g] = [L/Kg] (Allison and Allison, 2005 - EPA/600/R-05/074)


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
psi_crop = 110  # soil suction Alteck [cm]
psi_vine = 110  # soil suction Rouffach (guess)

""" Microcosm """
d = (1.493 * 2)  # Diameter of falcon tube (cm)
area = ((d / 2) ** 2) * 3.1416  # (cmÂ²)
zl = soil_height = 3.0  # Mixing layer depth in cm

################################################################
# Annual Crop ##################################################
################################################################

""" Hydorlogy calculation - Annual Crop Soil """
water_data = leachsim(ovSat=ovSat_crop,
                      kSat=kSat_crop,
                      psi=psi_crop)


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

""" Observed Hydrlogy Annual Crop """
leach_high_6min = np.array([16.253, 12.958, 17.536, 14.29])  # all at 6 min
leach_med_12min = np.array([10.089, 5.902, 13.981, 10.602])  # all at 12 min
leach_med_30min = np.array([49.197, 40.402, 45.772, 47.201])  # all at 30min
leach_low_30min = np.array([20.037, 17.508, 22.376, 20.085])  # all at 30min

hydroplot(data, leach_high_6min, leach_med_12min, leach_med_30min, leach_low_30min)


""" Copper simulation w Retardation Factor - 1st Pulse, Annual Crop """
# Initial mass
mass_ini_sterile = (1627 + 1107) / float(2)  # all intensities:{0d, 10d)
mass_ini_untreated = (1184 + 1177) / float(2)  # all intesities:{0d, 10d)

# Observed Cupper Output at 6min, 12min, 30min, 30min @ 135, 55, 55, 30mm/h
cu_sol_sterile = np.array([11.29, 11.63, 306.80, 21.08])
cu_sol_untreat = np.array([0.5, 1.405, 37.0, 1])

# Copper log(Kd) range: 0.1 - 7.0, max-mean = 5.5
Kd_copper_1 = 10**0.65  # [mL/g] = [L/Kg] (Allison and Allison, 2005 - EPA/600/R-05/074)
Kd_copper_2 = 10**0.6  # [mL/g] = [L/Kg] (Allison and Allison, 2005 - EPA/600/R-05/074)
Kd_copper_3 = 10**0.55
Kd_copper_4 = 10**0.1
Kd_copper_5 = 10**0.05
Kd_copper_6 = 10**0.001

Kd_copper = [Kd_copper_1, Kd_copper_2, Kd_copper_3, Kd_copper_4, Kd_copper_5]
mass_ini = [mass_ini_sterile, mass_ini_untreated]

cu_sol = [cu_sol_sterile, cu_sol_untreat]

cum_copper_Alteck = pest_test(Kd_copper,
                              pb_crop,
                              ovSat_crop,
                              water_data,
                              area, soil_height,
                              mass_ini,
                              cu_sol)

pestiplot(cum_copper_Alteck, cu_sol_sterile, cu_sol_untreat, 'Cumulative Cu - Annual Crop Soil')


cum_copper_Alteck = freundlich(Kd_copper, 0.87,
                              pb_crop,
                              ovSat_crop,
                              water_data,
                              area, soil_height,
                              mass_ini,
                              cu_sol)

pestiplot(cum_copper_Alteck, cu_sol_sterile, cu_sol_untreat, 'Cumulative Cu - Annual Crop Soil, Freundlich')


"""
################
# Zinc simulation w Retardation Factor - 1st Pulse

mass_ini_sterile = (3106 + 2594) / float(2)  # all intensities:{0d, 10d)
mass_ini_untreated = (2636 + 2586) / float(2)

Kd_zinc = [Kd_zinc_sterile, Kd_zinc_untreat]
mass_ini = [mass_ini_sterile, mass_ini_untreated]
cum_zinc_Alteck = pesti_ret(Kd_zinc,
                            pb_crop,
                            ovSat_crop,
                            water_data,
                            area, soil_height,
                            mass_ini)

# Observed Zinc Output
zn_sol_sterile = np.array([17, 17, 404.25, 30.13])
zn_sol_untreat = np.array([2.19, 2, 36.85, 3.09])

pestiplot(cum_zinc_Alteck, zn_sol_sterile, zn_sol_untreat, 'Cumulative Zn - Annual Crop Soil')

################
# Metalaxyl simulation w Retardation Factor - 1st Pulse, Annual Crop
# Initial mass
mx_ini_sterile = (1818.1+1472.7) / float(2)
mx_ini_untreated = (1518.1+1413.3) / float(2)

Kd_mexyl_crop = [Kd_mexyl_crop_sterile, Kd_mexyl_crop_untreat]
mass_ini = [mx_ini_sterile, mx_ini_untreated]
cum_mx_crop = pesti_ret(Kd_mexyl_crop,
                        pb_crop,
                        ovSat_crop,
                        water_data,
                        area, soil_height,
                        mass_ini)

# Observed Metalaxyl Output - Crop
mx_obs_sterile_crop = np.array([(138.1+207.1)/2.0,
                                (201.0+50.4)/2.0, (641.8+356.8)/2.0,
                                (177.0+293.5)/2.0])  # high, med-12, med-30, low
mx_obs_untreat_crop = np.array([(145.4+283.5)/2.0,
                                (158.4+262.3)/2.0, (674.9+360.2)/2.0,
                                (418.2+480.9)/2.0])

pestiplot(cum_mx_crop, mx_obs_sterile_crop, mx_obs_untreat_crop, 'Cumulative Metalaxyl - Annual Crop Soil')

################
# S-metolachlor simulation w Retardation Factor - 1st Pulse, Annual Crop

# Initial mass
mr_ini_sterile = (5176.9+4213.1)/float(2)  # all intensities:{0d, 10d)
mr_ini_untreated = (3460.8+2832.7)/float(2)  # all intesities:{0d, 10d)

Kd_smeto_crop = [Kd_smeto_crop_sterile, Kd_smeto_crop_untreat]
mass_ini = [mr_ini_sterile, mr_ini_untreated]
cum_mr_crop = pesti_ret(Kd_smeto_crop,
                        pb_crop,
                        ovSat_crop,
                        water_data,
                        area, soil_height,
                        mass_ini)

# Observed Metolachlor Output - Crop
mr_obs_sterile_crop = np.array([(65.7+77.7)/2.0,
                                (79.5+16.9)/2.0, (327.4+153.4)/2.0,
                                (70.9+110.7)/2.0])  # high, med-12, med-30, low
mr_obs_untreat_crop = np.array([(53.5+89.2)/2.0,
                                (49.4+70.5)/2.0, (281.7+116.4)/2.0,
                                (142.7+127.2)/2.0])
pestiplot(cum_mr_crop, mr_obs_sterile_crop, mr_obs_untreat_crop, 'Cumulative S-metolachlor - Annual Crop Soil')

"""

"""
################################################################
# Vineyard #####################################################
################################################################
# Hydorlogy calculation - Vineyard Crop Soil
water_data = leachsim(ovSat=ovSat_vine,
                      kSat=kSat_vine,
                      psi=psi_crop)

cum_time_30min = water_data[:, 0]
cum_inf_135mmh = water_data[:, 4]
cum_inf_55mmh = water_data[:, 5]
cum_inf_30mmh = water_data[:, 6]
cum_leach_135mmh = water_data[:, 7]
cum_leach_55mmh = water_data[:, 8]
cum_leach_30mmh = water_data[:, 9]

dataVine = stackdata6(cum_time_30min,
                 cum_inf_135mmh, cum_inf_55mmh, cum_inf_30mmh,
                 cum_leach_135mmh, cum_leach_55mmh, cum_leach_30mmh)

# Observed Hydrlogy Vineyard
leach_high_6min = np.array([13.609, 13.610, 17.676, 17.705])  # all at 6 min
leach_med_12min = np.array([13.787, 11.112, 11.858, 11.294])  # all at 12 min
leach_med_30min = np.array([48.185, 46.402, 48.164, 47.032])  # all at 30min
leach_low_30min = np.array([22.595, 19.082, 21.285, 20.871])  # all at 30min

hydroplot(dataVine, leach_high_6min, leach_med_12min, leach_med_30min, leach_low_30min)

###########
# Copper simulation w Retardation Factor - 1st Pulse, Vineyard Crop
# Initial mass
cu_ini_sterile = (3735 + 4008) / float(2)  # all intensities:{0d, 10d)
cu_ini_untreated = (4011 + 3860) / float(2)  # all intesities:{0d, 10d)

cum_copper_vine = pesti_ret(Kd_copper_sterile, Kd_copper_untreat,
                            pb_vine,
                            ovSat_vine,
                            water_data,
                            area, soil_height,
                            cu_ini_sterile, cu_ini_untreated)

# Observed Cupper Output - Vine
cu_obs_sterile_vine = np.array([(23.4+38.5)/2.0,
                                (28.3+29.7)/2.0, (611.1+588.5)/2.0,
                                (48.4+51.8)/2.0])  # high, med, med, low
cu_obs_untreat_vine = np.array([(1.4+0.6)/2.0,
                                (1.0+0.4)/2.0, (47.3+33.3)/2.0,
                                (2.7+0.8)/2.0])

pestiplot(cum_copper_vine, cu_obs_sterile_vine, cu_obs_untreat_vine, 'Cumulative Cu - Vineyard Soil')

###########
# Zinc simulation w Retardation Factor - 1st Pulse, Vineyard Crop
# Initial mass
zn_ini_sterile = (2811 + 3122) / float(2)
zn_ini_untreated = (3216 + 3031) / float(2)


cum_zinc_vine = pesti_ret(Kd_copper_sterile, Kd_copper_untreat,
                          pb_vine,
                          ovSat_vine,
                          water_data,
                          area, soil_height,
                          zn_ini_sterile, zn_ini_untreated)


# Observed Zinc Output - Vine
zn_obs_sterile_vine = np.array([(0.6+2.4)/2.0,
                                (0.8+0)/2.0, (24.2+4.3)/2.0,
                                (1.8+0)/2.0])  # high, med-12, med-30, low
zn_obs_untreat_vine = np.array([0,
                                0, (0.5+10.0)/2.0,
                                0])

pestiplot(cum_zinc_vine, zn_obs_sterile_vine, zn_obs_untreat_vine, 'Cumulative Zn - Vineyard Soil')

###########
# Metalaxyl simulation w Retardation Factor - 1st Pulse, Vineyard Crop
# Initial mass
mx_ini_sterile = (2105.8+1504.9) / float(2)  # all intensities:{0d, 10d)
mx_ini_untreated = (1741.7+1696.7) / float(2)  # all intesities:{0d, 10d)

cum_mx_vine = pesti_ret(Kd_mexyl_vine_sterile, Kd_mexyl_vine_untreat,
                        pb_vine,
                        ovSat_vine,
                        water_data,
                        area, soil_height,
                        mx_ini_sterile, mx_ini_untreated)

# Observed Metalaxyl Output - Vine
mx_obs_sterile_vine = np.array([(307.0+269.5)/2.0,
                                (279.2+294.9)/2.0, (313.8+529.8)/2.0,
                                (419.4+346.1)/2.0])  # high, med-12, med-30, low
mx_obs_untreat_vine = np.array([(271.3+116.9)/2.0,
                                (304.3+276.9)/2.0, (616.4+409.1)/2.0,
                                (365.0+429.5)/2.0])

pestiplot(cum_mx_vine, mx_obs_sterile_vine, mx_obs_untreat_vine, 'Cumulative Metalaxyl - Vineyard Soil')

###########
# S-metolachlor simulation w Retardation Factor - 1st Pulse, Vineyard Crop
# Initial mass
mr_ini_sterile = (7343.7+3052.6)/float(2)  # all intensities:{0d, 10d)
mr_ini_untreated = (4038.5+4904.5)/float(2)  # all intesities:{0d, 10d)

cum_mr_vine = pesti_ret(Kd_smeto_vine_sterile, Kd_smeto_vine_untreat,
                        pb_vine,
                        ovSat_vine,
                        water_data,
                        area, soil_height,
                        mr_ini_sterile, mr_ini_untreated)

# Observed Metolachlor Output - Vine
mr_obs_sterile_vine = np.array([(148.5+111.5)/2.0,
                                (160.1+110.7)/2.0, (188.3+235.8)/2.0,
                                (237.4+142.8)/2.0])  # high, med-12, med-30, low
mr_obs_untreat_vine = np.array([(111.8+37.2)/2.0,
                                (111.3+79.3)/2.0, (278.4+156.1)/2.0,
                                (142.0+134.7)/2.0])
pestiplot(cum_mr_vine, mr_obs_sterile_vine, mr_obs_untreat_vine, 'Cumulative S-metolachlor - Vineyard Soil')


"""
