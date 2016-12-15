from math import exp
from hydroplots import *
import numpy as np
from mixinglayer import *


def copper_ret(Kd_min, Kd_max, pb, ov_sat, water_data, area, soil_height):
    """
    Water data:
    """
    vol_h2o_sat = area * soil_height * ov_sat
    cum_time_30min = water_data[:, 0]
    leach_6min = water_data[:, 1]  # leached volume per timestep @ high intensity
    leach_12min = water_data[:, 2]
    leach_30min = water_data[:, 3]
    leached_vol = [leach_6min, leach_12min, leach_30min]

    """ Copper Lab Data"""
    #  Mass initial (untreated)
    mass_ini_sterile = (1627 + 1107) / float(2)  # all intensities:{0d, 10d)
    mass_ini_untreated = (1184 + 1177) / float(2)  # all intesities:{0d, 10d)

    """
    Simulations for min Kd (three rainfall intesities)
    Kd_min is matched to iniital mass in the sterile soils (mass_ini_sterile: {higher output})
    """
    conc_old = mass_ini_sterile/vol_h2o_sat
    r_factor_min = 1 + (pb * Kd_min) / ov_sat
    #  conc_old[i] : initial concentration in soil for each rainfall event replica
    for i in range(3):
        mass_left_old = conc_old * vol_h2o_sat
        cum_mass_cu_out = 0
        cum_mass_cu_out_dt = []
        for t in range(len(cum_time_30min)):
            #  water_data[:, t] = leached volume at time-step t
            conc_new = conc_old*exp(-(leached_vol[i][t]) / (r_factor_min*vol_h2o_sat))
            mass_left_new = conc_new*vol_h2o_sat
            mass_out = mass_left_old - mass_left_new
            cum_mass_cu_out += mass_out
            cum_mass_cu_out_dt.append(cum_mass_cu_out)
            conc_old = conc_new
        if i == 0:
            highint_cum_mass_st_cu_out_dt = cum_mass_cu_out_dt
            #print("cum mass out @ 6min, I=135 mm/h: ", cum_mass_cu_out_dt[23])
        elif i == 1:
            medint_cum_mass_st_cu_out_dt = cum_mass_cu_out_dt
            #print("cum mass out @ 12min, I=55 mm/h: ", cum_mass_cu_out_dt[47])
            #print("cum mass out @ 30min, I=55 mm/h: ", cum_mass_cu_out_dt[-1])
        elif i == 2:
            lowint_cum_mass_st_cu_out_dt = cum_mass_cu_out_dt
            #print("cum mass out @ 30min, I=30 mm/h: ", cum_mass_cu_out_dt[-1])
        else:
            print("Do you have more than three Cu leaching scenarios?")

    """
    Simulations for max Kd (three rainfall intesities)
    Kd_max is matched to iniital mass in the untreated soils (mass_ini_sterile: {lower output})
    """
    conc_old = mass_ini_untreated / vol_h2o_sat
    r_factor = 1 + (pb * Kd_max) / ov_sat
    for i in range(3):
        mass_left_old = conc_old * vol_h2o_sat
        cum_mass_cu_out = 0
        cum_mass_cu_out_dt = []
        for t in range(len(cum_time_30min)):
            #  water_data[:, t] = leached volume at time-step t
            conc_new = conc_old * exp(-(leached_vol[i][t]) / (r_factor * vol_h2o_sat))
            mass_left_new = conc_new * vol_h2o_sat
            mass_out = mass_left_old - mass_left_new
            cum_mass_cu_out += mass_out
            cum_mass_cu_out_dt.append(cum_mass_cu_out)
            conc_old = conc_new
        if i == 0:
            highint_cum_mass_un_cu_out_dt = cum_mass_cu_out_dt
            #print("cum mass out @ 6min, I=135 mm/h: ", cum_mass_cu_out_dt[23])
        elif i == 1:
            medint_cum_mass_un_cu_out_dt = cum_mass_cu_out_dt
            #print("cum mass out @ 12min, I=55 mm/h: ", cum_mass_cu_out_dt[47])
            #print("cum mass out @ 30min, I=55 mm/h: ", cum_mass_cu_out_dt[-1])
        elif i == 2:
            lowint_cum_mass_un_cu_out_dt = cum_mass_cu_out_dt
            #print("cum mass out @ 30min, I=30 mm/h: ", cum_mass_cu_out_dt[-1])
        else:
            print("Do you have more than three Cu leaching scenarios?")

    return stackdata6(cum_time_30min,
                     highint_cum_mass_st_cu_out_dt, medint_cum_mass_st_cu_out_dt, lowint_cum_mass_st_cu_out_dt,
                     highint_cum_mass_un_cu_out_dt, medint_cum_mass_un_cu_out_dt, lowint_cum_mass_un_cu_out_dt)


