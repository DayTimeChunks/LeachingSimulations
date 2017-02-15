from math import exp, log10
from hydroplots import *
from mixinglayer import *
import matplotlib.pyplot as plt


def pesti_ret(Kd,  # Kd_sterile, Kd_untreat,
              pb,
              ov_sat,
              water_data,
              area, soil_height,
              mass_ini):

    """
    Water data: leaching, 1st pulse
    """
    vol_h2o_sat = area * soil_height * ov_sat
    cum_time_30min = water_data[:, 0]
    leach_135mmh = water_data[:, 1]  # leached volume per timestep @ high intensity
    leach_55mmh = water_data[:, 2]
    leach_30mmh = water_data[:, 3]
    leached_vol = [leach_135mmh, leach_55mmh, leach_30mmh]
    soil_m = 54.047  # [g]

    """ Cummualtive Output series """
    highint_cum_mass_st_out_dt = []
    medint_cum_mass_st_out_dt = []
    lowint_cum_mass_st_out_dt = []

    highint_cum_mass_un_out_dt = []
    medint_cum_mass_un_out_dt = []
    lowint_cum_mass_un_out_dt = []

    """ Output series """
    highint_mass_st_out_dt = []
    medint_mass_st_out_dt = []
    lowint_mass_st_out_dt = []

    highint_mass_un_out_dt = []
    medint_mass_un_out_dt = []
    lowint_mass_un_out_dt = []

    for k in range(len(Kd)):  # for each Kd tested (sterile & untreat)
        r_factor = 1 + (pb * Kd[k]) / ov_sat
        for i in range(3):
            conc_soil_old = mass_ini[k] / soil_m
            conc_liq_old = conc_soil_old / Kd[k]
            mass_old = conc_soil_old * soil_m
            cum_mass_out = 0
            mass_out_dt = []
            cum_mass_out_dt = []
            for t in range(len(cum_time_30min)):
                # McGrath leaching:
                conc_liq_new = conc_liq_old * exp(-(leached_vol[i][t]) / (r_factor * vol_h2o_sat))
                mass_out = (conc_liq_old - conc_liq_new) * vol_h2o_sat
                # re-equilibrate
                mass_new = mass_old - mass_out
                conc_soil_old = mass_new / soil_m
                conc_liq_old = conc_soil_old / Kd[k]
                cum_mass_out += mass_out
                mass_out_dt.append(mass_out)
                cum_mass_out_dt.append(cum_mass_out)
                mass_old = mass_new

            if i == 0 and k == 0:
                highint_cum_mass_st_out_dt = cum_mass_out_dt
                highint_mass_st_out_dt = mass_out_dt
                print("Mass balance for intensity ", i + 1,
                      "and Kd: ", Kd[k], "is ",
                      abs(mass_ini[k] - (cum_mass_out_dt[-1] + mass_old)) < 1 * 10 ** (-6))
            elif i == 1 and k == 0:
                medint_cum_mass_st_out_dt = cum_mass_out_dt
                medint_mass_st_out_dt = mass_out_dt
                print("Mass balance for intensity ", i + 1,
                      "and Kd: ", Kd[k], "is ",
                      abs(mass_ini[k] - (cum_mass_out_dt[-1] + mass_old)) < 1 * 10 ** (-6))
            elif i == 2 and k == 0:
                lowint_cum_mass_st_out_dt = cum_mass_out_dt
                lowint_mass_st_out_dt = mass_out_dt
                print("Mass balance for intensity ", i + 1,
                      "and Kd: ", Kd[k], "is ",
                      abs(mass_ini[k] - (cum_mass_out_dt[-1] + mass_old)) < 1 * 10 ** (-6))
            elif i == 0 and k == 1:
                highint_cum_mass_un_out_dt = cum_mass_out_dt
                highint_mass_un_out_dt = mass_out_dt
                print("Mass balance for intensity ", i + 1,
                      "and Kd: ", Kd[k], "is ",
                      abs(mass_ini[k] - (cum_mass_out_dt[-1] + mass_old)) < 1 * 10 ** (-6))
            elif i == 1 and k == 1:
                medint_cum_mass_un_out_dt = cum_mass_out_dt
                medint_mass_un_out_dt = mass_out_dt
                print("Mass balance for intensity ", i + 1,
                      "and Kd: ", Kd[k], "is ",
                      abs(mass_ini[k] - (cum_mass_out_dt[-1] + mass_old)) < 1 * 10 ** (-6))
            elif i == 2 and k == 1:
                lowint_cum_mass_un_out_dt = cum_mass_out_dt
                lowint_mass_un_out_dt = mass_out_dt
                print("Mass balance for intensity ", i + 1,
                      "and Kd: ", Kd[k], "is ",
                      abs(mass_ini[k] - (cum_mass_out_dt[-1] + mass_old)) < 1 * 10 ** (-6))
            else:
                print("Do you have more than 3 leaching intensities or 2 Kd's?")

    return [stackdata6(cum_time_30min,
                       highint_cum_mass_st_out_dt, medint_cum_mass_st_out_dt, lowint_cum_mass_st_out_dt,
                       highint_cum_mass_un_out_dt, medint_cum_mass_un_out_dt, lowint_cum_mass_un_out_dt),
            stackdata6(cum_time_30min,
                       highint_mass_st_out_dt, medint_mass_st_out_dt, lowint_mass_st_out_dt,
                       highint_mass_un_out_dt, medint_mass_un_out_dt, lowint_mass_un_out_dt)]


def pest_test(Kd,  # Kd list range to test
              intensities,
              pb,
              ov_sat,
              percol_data, pond_data, infil_data, time_sizes,
              area, soil_height,
              mass_ini,
              pest_sol_leach, pest_sol_pond,
              d, runoffvelocity,
              KFILM = True):
    """
    :param Kd: ensure units are mm3/g
    :param pb: ensure units are g/mm3
    :param ov_sat:
    :param percol_data:
    :param pond_data:
    :param infil_data:
    :param area:
    :param soil_height:
    :param mass_ini:
    :param pest_sol_leach:
    :param d:
    :param runoffvelocity:
    :return:
    """

    # Leaching, ponding and infiltration in vol (cm3)
    # are stored as lists to iterate through

    soil_m = 54.047  # [g] Average soil mass used in the experiment.

    vol_h2o_sat = area * soil_height * ov_sat  # mm3
    cum_time_30min = pond_data[:, 0]

    leach_135mmh = percol_data[:, 1]  # leached volume per timestep @ high intensity
    leach_55mmh = percol_data[:, 2]
    leach_30mmh = percol_data[:, 3]
    leached_vol = [leach_135mmh, leach_55mmh, leach_30mmh]

    pond_135mmh = pond_data[:, 1]  # ponded volume per timestep @ high intensity
    pond_55mmh = pond_data[:, 2]
    pond_30mmh = pond_data[:, 3]
    ponded_vol = [pond_135mmh, pond_55mmh, pond_30mmh]

    infil_135mmh = infil_data[:, 1]  # infil. volume per timestep @ high intensity
    infil_55mmh = infil_data[:, 2]
    infil_30mmh = infil_data[:, 3]
    infil_vol = [infil_135mmh, infil_55mmh, infil_30mmh]  # cm3

    """ Cummualtive Output series """
    highint_cum_mass_st_out_dt = [0.0] * len(cum_time_30min)
    medint_cum_mass_st_out_dt = [0.0] * len(cum_time_30min)
    lowint_cum_mass_st_out_dt = [0.0] * len(cum_time_30min)

    highint_cum_mass_un_out_dt = [0.0] * len(cum_time_30min)
    medint_cum_mass_un_out_dt = [0.0] * len(cum_time_30min)
    lowint_cum_mass_un_out_dt = [0.0] * len(cum_time_30min)

    """ Output series """
    highint_mass_st_out_dt = [0.0] * len(cum_time_30min)
    medint_mass_st_out_dt = [0.0] * len(cum_time_30min)
    lowint_mass_st_out_dt = [0.0] * len(cum_time_30min)

    highint_mass_un_out_dt = [0.0] * len(cum_time_30min)
    medint_mass_un_out_dt = [0.0] * len(cum_time_30min)
    lowint_mass_un_out_dt = [0.0] * len(cum_time_30min)

    for m in range(len(mass_ini)):  # for each initial mass (sterile & untreat)
        error = 10**9
        #  Check best Kd fit from selected range
        for k in range(len(Kd)):
            r_factor = 1 + (pb * Kd[k]) / ov_sat

            K_L = kfilm(d, runoffvelocity)  # mm/min
            #  Start a new intensity scenario
            for i in range(len(intensities)):
                mass_tot = mass_ini[m]
                # Whelan et al. 1987 (p4.9)
                conc_soil_old = \
                    ((mass_ini[m] / soil_m)*pb*Kd[k])/(ov_sat + pb*Kd[k])  # ug/g
                conc_liq_old = conc_soil_old / Kd[k]  # ug/mm3

                # Leaching
                cum_mass_out = 0
                mass_out_dt = []
                cum_mass_out_dt = []

                # Runoff
                mass_in_overflow = 0
                conc_in_overflow = 0
                mass_overflow_dt = []  # non cumulative as final value is total value
                #  Compute leaching for each time step under scenario "i"
                for t in range(len(cum_time_30min)):
                    # McGrath leaching loss:
                    conc_liq_new = conc_liq_old * exp(-(leached_vol[i][t]) / (r_factor * vol_h2o_sat))
                    # ug/mm3 * exp(- mm3/[-]*mm3) = ug/mm3
                    mass_liq_new = conc_liq_new*vol_h2o_sat  # ug
                    mass_out = (conc_liq_old - conc_liq_new) * vol_h2o_sat
                    # ug
                    cum_mass_out += mass_out

                    # Store temporal and cumulative mass leached
                    mass_out_dt.append(mass_out)
                    cum_mass_out_dt.append(cum_mass_out)

                    # Overland flow exchange/transfer
                    if ponded_vol[i][t] > 0:
                        if KFILM:
                            # Mixing layer model from Havis (1986) and Havis et al. (1992)
                            # Kfilm equation from Bennett and Myers, 1982; in Shi et al., 2011, p.1220)
                            # mass_transfered_overflow = K_L * (conc_liq_new-conc_in_overflow) * area * time_sizes[i][t]
                            mass_transfered_overflow = K_L * conc_liq_new * area * time_sizes[i][t]
                            # mm/min * ug/mm3 * mm2 * dt = [ug]
                            if mass_transfered_overflow < 0:
                                if conc_in_overflow*ponded_vol[i][t] < mass_transfered_overflow:
                                    mass_transfered_overflow = conc_in_overflow*ponded_vol[i][t]
                                else:
                                    print ("mass in of: ", conc_in_overflow*ponded_vol[i][t])
                                    print ("mass -> pond: ", mass_transfered_overflow, "at: ", t)
                                    print ("conc liq (t-1): ", conc_liq_old)
                                    print ("conc liq (t): ", conc_liq_new)
                                    print("conc in of (t-1): ", mass_overflow_dt[t-1]/vol_h2o_sat)
                                    print ("conc in of (t): ", conc_in_overflow)

                            mass_in_overflow += mass_transfered_overflow
                            conc_in_overflow = mass_in_overflow/ponded_vol[i][t]  # ug/mm3
                            mass_liq_new -= mass_transfered_overflow
                            # conc_liq_new = mass_liq_new/vol_h2o_sat
                            # conc_liq_new += conc_in_overflow*infil_vol[i][t]/vol_h2o_sat
                            # mass_in_overflow -= conc_in_overflow*infil_vol[i][t]

#                        if mass_in_overflow < 0:
#                            print("Error in mixing layer approach, mass is negative")
                        else:
                            print("Error, over land flow transfer not defined")

                    # Store the temporal mass in overflow
                    mass_overflow_dt.append(mass_in_overflow)

                    # re-equilibrate
                    mass_tot_new = mass_tot - mass_out
                    conc_soil_old = (((mass_tot_new - mass_in_overflow) / soil_m)*pb*Kd[k])/(ov_sat + pb*Kd[k])
                    conc_liq_old = conc_soil_old / Kd[k]
                    mass_tot = mass_tot_new

                # Calculate error of intensity i, @ (kd, m)
                if i == 0 and m == 0:
                    highint_error6_st = (cum_mass_out_dt[5] - pest_sol_leach[m][i]) ** 2 + \
                                        (mass_overflow_dt[5] - pest_sol_pond[m][i])**2
                    temp_out_high_st = mass_out_dt
                    temp_cum_out_high_st = cum_mass_out_dt
                    temp_overflow_high_st = mass_overflow_dt
                    # print("Mass balance for intensity ", i + 1, "and Kd: ", Kd[k], "is ",
                    #      abs(mass_ini[m] - mass_tot - cum_mass_out_dt[-1]) < 1*10**-6)
                elif i == 1 and m == 0:
                    medint_error12_st = (cum_mass_out_dt[11] - pest_sol_leach[m][i]) ** 2 + \
                                        (mass_overflow_dt[11] - pest_sol_pond[m][i])**2
                    medint_error30_st = (cum_mass_out_dt[-1] - pest_sol_leach[m][i + 1]) ** 2 + \
                                        (mass_overflow_dt[-1] - pest_sol_pond[m][i + 1])**2
                    temp_cum_out_med_st = cum_mass_out_dt
                    temp_out_med_st = mass_out_dt
                    temp_overflow_med_st = mass_overflow_dt
                    # print("Mass balance for intensity ", i + 1, "and Kd: ", Kd[k], "is ",
                    #      abs(mass_ini[m] - mass_tot - cum_mass_out_dt[-1]) < 1 * 10 ** -6)
                elif i == 2 and m == 0:
                    lowint_error30_st = (cum_mass_out_dt[-1] - pest_sol_leach[m][i + 1]) ** 2 + \
                                        (mass_overflow_dt[-1] - pest_sol_pond[m][i + 1]) ** 2
                    temp_cum_out_low_st = cum_mass_out_dt
                    temp_out_low_st = mass_out_dt
                    temp_overflow_low_st = mass_overflow_dt
                elif i == 0 and m == 1:
                    highint_error6_un = (cum_mass_out_dt[5] - pest_sol_leach[m][i]) ** 2 + \
                                        (mass_overflow_dt[5] - pest_sol_pond[m][i]) ** 2
                    temp_cum_out_high_un = cum_mass_out_dt
                    temp_out_high_un = mass_out_dt
                    temp_overflow_high_un = mass_overflow_dt
                elif i == 1 and m == 1:
                    medint_error12_un = (cum_mass_out_dt[11] - pest_sol_leach[m][i]) ** 2 + \
                                        (mass_overflow_dt[11] - pest_sol_pond[m][i])**2
                    medint_error30_un = (cum_mass_out_dt[-1] - pest_sol_leach[m][i + 1]) ** 2 + \
                                        (mass_overflow_dt[-1] - pest_sol_pond[m][i+1])**2
                    temp_cum_out_med_un = cum_mass_out_dt
                    temp_out_med_un = mass_out_dt
                    temp_overflow_med_un = mass_overflow_dt
                elif i == 2 and m == 1:
                    lowint_error30_un = (cum_mass_out_dt[-1] - pest_sol_leach[m][i + 1]) ** 2 + \
                                        (mass_overflow_dt[-1] - pest_sol_pond[m][i+1])**2
                    temp_cum_out_low_un = cum_mass_out_dt
                    temp_out_low_un = mass_out_dt
                    temp_overflow_low_un = mass_overflow_dt

            if m == 0:
                error_test = (((highint_error6_st +
                              medint_error12_st + medint_error30_st +
                              lowint_error30_st)/4) ** 0.5)/(max(pest_sol_leach[m])-min(pest_sol_leach[m]))

                if error_test < error:
                    error = error_test
                    error_st = error
                    kd_sterile = Kd[k]
                    num_st = k
                    # print("error = ", error,
                    #      "log Kd (sterile) = ", log10(kd_sterile), "Num: ", num_st+1)
                    highint_cum_mass_st_out_dt = temp_cum_out_high_st
                    highint_mass_st_out_dt = temp_out_high_st
                    highint_overmass_st_dt = temp_overflow_high_st
                    medint_cum_mass_st_out_dt = temp_cum_out_med_st
                    medint_mass_st_out_dt = temp_out_med_st
                    medint_overmass_st_dt = temp_overflow_med_st
                    lowint_cum_mass_st_out_dt = temp_cum_out_low_st
                    lowint_mass_st_out_dt = temp_out_low_st
                    lowint_overmass_st_dt = temp_overflow_low_st
                else:
                    continue

            elif m == 1:
                error_test = (((highint_error6_un +
                              medint_error12_un + medint_error30_un +
                              lowint_error30_un)/4) ** 0.5)/(max(pest_sol_leach[m])-min(pest_sol_leach[m]))
                if error_test < error:
                    error = error_test
                    error_un = error
                    kd_untreat = Kd[k]
                    num_un = k
                    # print("error = ", error,
                    #      "log Kd (untreat) = ", log10(kd_untreat), "Num: ", num_un+1)
                    highint_cum_mass_un_out_dt = temp_cum_out_high_un
                    highint_mass_un_out_dt = temp_out_high_un
                    highint_overmass_un_dt = temp_overflow_high_un
                    medint_cum_mass_un_out_dt = temp_cum_out_med_un
                    medint_mass_un_out_dt = temp_out_med_un
                    medint_overmass_un_dt = temp_overflow_med_un
                    lowint_cum_mass_un_out_dt = temp_cum_out_low_un
                    lowint_mass_un_out_dt = temp_out_low_un
                    lowint_overmass_un_dt = temp_overflow_low_un
                else:
                    continue

    print("Best log Kd (sterile): ", log10(kd_sterile), "( Num: ", num_st+1, ")", "\n",
          "Error: ", error_st)
    print("Best log Kd (untreat): ", log10(kd_untreat), "( Num: ", num_un+1, ")", "\n",
          "Error: ", error_un)

    return stackdata18(cum_time_30min,
                       highint_cum_mass_st_out_dt, medint_cum_mass_st_out_dt, lowint_cum_mass_st_out_dt,
                       highint_cum_mass_un_out_dt, medint_cum_mass_un_out_dt, lowint_cum_mass_un_out_dt,
                       highint_mass_st_out_dt, medint_mass_st_out_dt, lowint_mass_st_out_dt,
                       highint_mass_un_out_dt, medint_mass_un_out_dt, lowint_mass_un_out_dt,
                       highint_overmass_st_dt, medint_overmass_st_dt, lowint_overmass_st_dt,
                       highint_overmass_un_dt, medint_overmass_un_dt, lowint_overmass_un_dt)


def pest_test2(Kd,  # Kd list range to test
               pest_dict,  # Initial mass and observed data
               pb,
               ov_sat,
               percol_data, pond_data, time_sizes,
               area, soil_height,
               d, runoffvelocity,
               KFILM = True):
    """
    This function takes only four hydrological time series,
    based on a single Ksat.

    :param Kd: ensure units are mm3/g
    :param pb: ensure units are g/mm3
    :param pest_dict:
    :param ov_sat:
    :param percol_data:
    :param pond_data:
    :param time_sizes:
    :param area:
    :param soil_height:
    :param d:
    :param runoffvelocity:
    :param KFILM:
    :return: 16 time series of leached (0d and 10d) and ponded mass (0d and 10d).
    """

    # Leaching, ponding and infiltration in vol (cm3)
    # are stored as lists to iterate through

    soil_m = 54.047  # [g] Average soil mass used in the experiment.

    vol_h2o_sat = area * soil_height * ov_sat  # mm3
    cum_time_30min = pond_data[:, 0]

    leach_135mmh = percol_data[:, 1]  # leached volume per timestep @ high intensity
    leach_55mmhA = percol_data[:, 2]
    leach_55mmhB = percol_data[:, 3]
    leach_30mmh = percol_data[:, 4]
    leached_vol = [leach_135mmh, leach_55mmhA, leach_55mmhB, leach_30mmh]

    pond_135mmh = pond_data[:, 1]  # ponded volume per timestep @ high intensity
    pond_55mmhA = pond_data[:, 2]
    pond_55mmhB = pond_data[:, 3]
    pond_30mmh = pond_data[:, 4]
    ponded_vol = [pond_135mmh, pond_55mmhA, pond_55mmhB, pond_30mmh]

    # infil_135mmh = infil_data[:, 1]  # infil. volume per timestep @ high intensity
    # infil_55mmh = infil_data[:, 2]
    # infil_30mmh = infil_data[:, 3]
    # infil_vol = [infil_135mmh, infil_55mmh, infil_30mmh]  # cm3

    """ Cummualtive Output series """
    highint_cum_mass_st_out_dt = [0.0] * len(cum_time_30min)
    medint_cum_mass_st_out_dt = [0.0] * len(cum_time_30min)
    lowint_cum_mass_st_out_dt = [0.0] * len(cum_time_30min)

    highint_cum_mass_un_out_dt = [0.0] * len(cum_time_30min)
    medint_cum_mass_un_out_dt = [0.0] * len(cum_time_30min)
    lowint_cum_mass_un_out_dt = [0.0] * len(cum_time_30min)

    """ Output series """
    highint_mass_st_out_dt = [0.0] * len(cum_time_30min)
    medint_mass_st_out_dt = [0.0] * len(cum_time_30min)
    lowint_mass_st_out_dt = [0.0] * len(cum_time_30min)

    highint_mass_un_out_dt = [0.0] * len(cum_time_30min)
    medint_mass_un_out_dt = [0.0] * len(cum_time_30min)
    lowint_mass_un_out_dt = [0.0] * len(cum_time_30min)

    # Check best Kd fit from selected range
    error = 10 ** 9  # define a max error
    for k in range(len(Kd)):
        r_factor = 1 + (pb * Kd[k]) / ov_sat
        K_L = kfilm(d, runoffvelocity)  # mm/min

        #  for index, (name, (mass_i, leach_obs, pond_obs)) in enumerate(pest_dict.iteritems()): # python 2.7
        for index, (name, (mass_i, leach_obs, pond_obs)) in enumerate(pest_dict.items()):  # python 3.

            # Assign appropriate rainfall intensity reference
            if index < 3:
                i = 0
            elif index < 6:
                i = 1
            else:
                i = 2

            mass_tot = mass_i
            # Whelan et al. 1987 (p4.9)
            conc_soil_old = \
                ((mass_tot/ soil_m)*pb*Kd[k])/(ov_sat + pb*Kd[k])  # ug/g
            conc_liq_old = conc_soil_old / Kd[k]  # ug/mm3

            # Leaching
            cum_mass_out = 0
            mass_out_dt = []
            cum_mass_out_dt = []

            # Runoff
            mass_in_overflow = 0
            conc_in_overflow = 0
            mass_overflow_dt = []  # non cumulative as final value is total value

            #  Compute leaching for each time step under scenario of "index"
            for t in range(len(cum_time_30min)):
                # McGrath leaching loss:
                conc_liq_new = conc_liq_old * exp(-(leached_vol[i][t]) / (r_factor * vol_h2o_sat))
                # ug/mm3 * exp(- mm3/[-]*mm3) = ug/mm3
                mass_liq_new = conc_liq_new*vol_h2o_sat  # ug
                mass_out = (conc_liq_old - conc_liq_new) * vol_h2o_sat
                # ug
                cum_mass_out += mass_out

                # Store temporal and cumulative mass leached
                mass_out_dt.append(mass_out)
                cum_mass_out_dt.append(cum_mass_out)

                # Overland flow exchange/transfer
                if ponded_vol[i][t] > 0:
                    if KFILM:
                        # Mixing layer model from Havis (1986) and Havis et al. (1992)
                        # Kfilm equation from Bennett and Myers, 1982; in Shi et al., 2011, p.1220)
                        # mass_transfered_overflow = K_L * (conc_liq_new-conc_in_overflow) * area * time_sizes[i][t]
                        mass_transfered_overflow = K_L * conc_liq_new * area * time_sizes[i][t]
                        # mm/min * ug/mm3 * mm2 * dt = [ug]
                        if mass_transfered_overflow < 0:
                            if conc_in_overflow*ponded_vol[i][t] < mass_transfered_overflow:
                                mass_transfered_overflow = conc_in_overflow*ponded_vol[i][t]
                            else:
                                print ("mass in of: ", conc_in_overflow*ponded_vol[i][t])
                                print ("mass -> pond: ", mass_transfered_overflow, "at: ", t)
                                print ("conc liq (t-1): ", conc_liq_old)
                                print ("conc liq (t): ", conc_liq_new)
                                print("conc in of (t-1): ", mass_overflow_dt[t-1]/vol_h2o_sat)
                                print ("conc in of (t): ", conc_in_overflow)

                        mass_in_overflow += mass_transfered_overflow
                        conc_in_overflow = mass_in_overflow/ponded_vol[i][t]  # ug/mm3
                        mass_liq_new -= mass_transfered_overflow
                        # conc_liq_new = mass_liq_new/vol_h2o_sat
                        # conc_liq_new += conc_in_overflow*infil_vol[i][t]/vol_h2o_sat
                        # mass_in_overflow -= conc_in_overflow*infil_vol[i][t]

#                        if mass_in_overflow < 0:
#                            print("Error in mixing layer approach, mass is negative")
                    else:
                        print("Error, over land flow transfer not defined")

                # Store the temporal mass in overflow
                mass_overflow_dt.append(mass_in_overflow)

                # re-equilibrate
                mass_tot_new = mass_tot - mass_out
                conc_soil_old = (((mass_tot_new - mass_in_overflow) / soil_m)*pb*Kd[k])/(ov_sat + pb*Kd[k])
                conc_liq_old = conc_soil_old / Kd[k]
                mass_tot = mass_tot_new

            if index == 0:
                high_0d_error = (cum_mass_out_dt[5] - leach_obs) ** 2 + \
                                    (mass_overflow_dt[5] - pond_obs)**2
                temp_out_high_0d = mass_out_dt
                temp_cum_out_high_0d = cum_mass_out_dt
                temp_overflow_high_0d = mass_overflow_dt

            # high int, 6 min; 0 days, living. Scenario A2 - h6L0d
            elif index == 1:
                high_1d_error = (cum_mass_out_dt[5] - leach_obs) ** 2 + \
                                    (mass_overflow_dt[5] - pond_obs)**2
                temp_out_high_1d = mass_out_dt
                temp_cum_out_high_1d = cum_mass_out_dt
                temp_overflow_high_1d = mass_overflow_dt

            # i(0), m(8) = high int, 6 min; 10 days, sterile. Scenario A3 - h6S1d
            elif index == 2:
                med12_0d_error = (cum_mass_out_dt[11] - leach_obs) ** 2 + \
                                    (mass_overflow_dt[11] - pond_obs)**2
                temp_out_med12_0d = mass_out_dt
                temp_cum_out_med12_0d = cum_mass_out_dt
                temp_overflow_med12_0d = mass_overflow_dt

            # i(0), m(12) = high int, 6 min; 10 days, living. Scenario A4 - h61d
            elif index == 3:
                med12_1d__error = (cum_mass_out_dt[11] - leach_obs) ** 2 + \
                                (mass_overflow_dt[11] - pond_obs) ** 2
                temp_out_med12_1d = mass_out_dt
                temp_cum_out_med12_1d = cum_mass_out_dt
                temp_overflow_med12_1d = mass_overflow_dt

            elif index == 4:
                med30_0d__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                  (mass_overflow_dt[-1] - pond_obs) ** 2
                temp_out_med30_0d = mass_out_dt
                temp_cum_out_med30_0d = cum_mass_out_dt
                temp_overflow_med30_0d = mass_overflow_dt
            elif index == 5:
                med30_1d__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                  (mass_overflow_dt[-1] - pond_obs) ** 2
                temp_out_med30_1d = mass_out_dt
                temp_cum_out_med30_1d = cum_mass_out_dt
                temp_overflow_med30_1d = mass_overflow_dt
            elif index == 6:
                low_0d__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                  (mass_overflow_dt[-1] - pond_obs) ** 2
                temp_out_low_0d = mass_out_dt
                temp_cum_out_low_0d = cum_mass_out_dt
                temp_overflow_low_0d = mass_overflow_dt
            elif index == 7:
                low_1d__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                  (mass_overflow_dt[-1] - pond_obs) ** 2
                temp_out_low_1d = mass_out_dt
                temp_cum_out_low_1d = cum_mass_out_dt
                temp_overflow_low_1d = mass_overflow_dt
            else:
                print("Index number error")

        error_test = (((high_0d_error + high_1d_error +
                        med12_0d_error + med12_1d__error +
                        med30_0d__error + med30_1d__error +
                        low_0d__error + low_1d__error)/8.0) ** 0.5)

        if error_test < error:
            error = error_test
            kd_chosen = Kd[k]
            num_kd = k
            # print("error = ", error,
            #      "log Kd = ", log10(kd_chosen), "Num: ", num_kd+1)
            # Store final variables

            # a) Cumulative mass leached
            high_0d_cum_mass_out_dt = temp_cum_out_high_0d
            high_1d_cum_mass_out_dt = temp_cum_out_high_1d
            med12_0d_cum_mass_out_dt = temp_cum_out_med12_0d
            med12_1d_cum_mass_out_dt = temp_cum_out_med12_1d
            med30_0d_cum_mass_out_dt = temp_cum_out_med30_0d
            med30_1d_cum_mass_out_dt = temp_cum_out_med30_1d
            low_0d_cum_mas_out_dt = temp_cum_out_low_0d
            low_1d_cum_mas_out_dt = temp_cum_out_low_1d

            # b) Mass in ponded water
            high_0d_overmass_dt = temp_overflow_high_0d
            high_1d_overmass_dt = temp_overflow_high_1d
            med12_0d_overmass_dt = temp_overflow_med12_0d
            med12_1d_overmass_dt = temp_overflow_med12_1d
            med30_0d_overmass_dt = temp_overflow_med30_0d
            med30_1d_overmass_dt = temp_overflow_med30_1d
            low_0d_overmass_dt = temp_overflow_low_0d
            low_1d_overmass_dt = temp_overflow_low_1d

        else:
            continue

    print("Best log Kd: ", log10(kd_chosen), "( Num: ", num_kd+1, ")", "\n",
          "Error: ", error)

    return stackdata16(cum_time_30min,
                       high_0d_cum_mass_out_dt, high_1d_cum_mass_out_dt,
                       med12_0d_cum_mass_out_dt, med12_1d_cum_mass_out_dt,
                       med30_0d_cum_mass_out_dt, med30_1d_cum_mass_out_dt,
                       low_0d_cum_mas_out_dt, low_1d_cum_mas_out_dt,
                       high_0d_overmass_dt, high_1d_overmass_dt,
                       med12_0d_overmass_dt, med12_1d_overmass_dt,
                       med30_0d_overmass_dt, med30_1d_overmass_dt,
                       low_0d_overmass_dt, low_1d_overmass_dt)


def pest_linear_x(
        Kd,  # Kd list range to test
        x,
        pest_dict,  # Initial mass and observed data
        pb_i, pb_f,
        percol_data, pond_data,
        time_sizes,
        area, soil_height,
        d, runoffvelocity,
        soil,
        isFirstCycle,
        isLiving,
        KFILM=True):
    """
    This function takes 8 hydrological time series,
    based optimal Ksat values.

    The dictionary is either Sterile or Living

    :param Kd: ensure units are mm3/g
    :param pb: ensure units are g/mm3
    :param pest_dict:
    :param ov_sat:
    :param percol_data: 8 time series
    :param pond_data: 8 time series
    :param time_sizes: 8 time series
    :param area:
    :param soil_height:
    :param d:
    :param runoffvelocity:
    :param KFILM:
    :return:
        16 time series of leached (0d and 10d) and ponded mass (0d and 10d).
        Sees no distinction between sterile or live
    """

    # Leaching, ponding and infiltration in vol (cm3)
    # are stored as lists to iterate through

    if isFirstCycle:
        if soil == 'Alteck':
            ovSat = 0.61
            ov = 0.20
        elif soil == 'Rouff':
            ovSat = 0.55
            ov = 0.20
    else:
        # Saturated water content (0.61 - Alteck; 0.55 - Rouff)
        if soil == 'Alteck':
            ovSat = 0.61
            ov = 0.61 - 0.039
        elif soil == 'Rouff':
            ovSat = 0.55
            ov = 0.55 - 0.038

    soil_m = 54.047  # [g] Average soil mass used in the experiment.

    vol_h2o_sat = area * soil_height * ovSat  # mm3
    cum_time_30min = pond_data[:, 0]

    min_mass = 200
    max_mass = 0
    sum_mass = 0
    masses = 0
    for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) \
            in enumerate(sorted(pest_dict.items())):
        # print (index, (name, (mass_i, leach_obs, pond_obs)))

        if float(leach_obs) > 0.:
            masses += 1
            sum_mass += float(leach_obs)
        if float(leach_obs) > max_mass:
            max_mass = float(leach_obs)
        if float(leach_obs) < min_mass:
            min_mass = float(leach_obs)

    mean_mass = sum_mass/masses

    SStot_fresh = 0
    SStot_aged = 0
    for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) in enumerate(pest_dict.items()):
        if index % 2 == 0:
            if float(leach_obs) > 0:
                SStot_fresh += (float(leach_obs) - mean_mass)**2
            else:
                SStot_fresh += mean_mass**2
        else:
            if float(leach_obs) > 0:
                SStot_aged += (float(leach_obs) - mean_mass)**2
            else:
                SStot_aged += mean_mass**2

    if isFirstCycle:
        leach_135mmh_fresh = percol_data[:, 1]  # leached volume per timestep @ high intensity
        leach_55mmhA_fresh = percol_data[:, 2]
        leach_55mmhB_fresh = percol_data[:, 2]
        leach_30mmh_fresh = percol_data[:, 3]

        leach_135mmh_aged = percol_data[:, 1]  # leached volume per timestep @ high intensity
        leach_55mmhA_aged = percol_data[:, 2]
        leach_55mmhB_aged = percol_data[:, 2]
        leach_30mmh_aged = percol_data[:, 3]

        pond_135mmh_fresh = pond_data[:, 1]  # ponded volume per timestep @ high intensity
        pond_55mmhA_fresh = pond_data[:, 2]
        pond_55mmhB_fresh = pond_data[:, 2]
        pond_30mmh_fresh = pond_data[:, 3]

        pond_135mmh_aged = pond_data[:, 1]  # ponded volume per timestep @ high intensity
        pond_55mmhA_aged = pond_data[:, 2]
        pond_55mmhB_aged = pond_data[:, 2]
        pond_30mmh_aged = pond_data[:, 3]

    else:
        print("Error? X-factor should not be used for 2nd cycle...")
        if isLiving:
            leach_135mmh_fresh = percol_data[:, 9]  # leached volume per timestep @ high intensity
            leach_55mmhA_fresh = percol_data[:, 10]
            leach_55mmhB_fresh = percol_data[:, 11]
            leach_30mmh_fresh = percol_data[:, 12]

            leach_135mmh_aged = percol_data[:, 13]  # leached volume per timestep @ high intensity
            leach_55mmhA_aged = percol_data[:, 14]
            leach_55mmhB_aged = percol_data[:, 15]
            leach_30mmh_aged = percol_data[:, 16]

            pond_135mmh_fresh = pond_data[:, 9]  # ponded volume per timestep @ high intensity
            pond_55mmhA_fresh = pond_data[:, 10]
            pond_55mmhB_fresh = pond_data[:, 11]
            pond_30mmh_fresh = pond_data[:, 12]

            pond_135mmh_aged = pond_data[:, 13]  # ponded volume per timestep @ high intensity
            pond_55mmhA_aged = pond_data[:, 14]
            pond_55mmhB_aged = pond_data[:, 15]
            pond_30mmh_aged = pond_data[:, 16]
        else:
            leach_135mmh_fresh = percol_data[:, 1]  # leached volume per timestep @ high intensity
            leach_55mmhA_fresh = percol_data[:, 2]
            leach_55mmhB_fresh = percol_data[:, 3]
            leach_30mmh_fresh = percol_data[:, 4]

            leach_135mmh_aged = percol_data[:, 5]  # leached volume per timestep @ high intensity
            leach_55mmhA_aged = percol_data[:, 6]
            leach_55mmhB_aged = percol_data[:, 7]
            leach_30mmh_aged = percol_data[:, 8]

            pond_135mmh_fresh = pond_data[:, 1]  # ponded volume per timestep @ high intensity
            pond_55mmhA_fresh = pond_data[:, 2]
            pond_55mmhB_fresh = pond_data[:, 3]
            pond_30mmh_fresh = pond_data[:, 4]

            pond_135mmh_aged = pond_data[:, 5]  # ponded volume per timestep @ high intensity
            pond_55mmhA_aged = pond_data[:, 6]
            pond_55mmhB_aged = pond_data[:, 7]
            pond_30mmh_aged = pond_data[:, 8]

    leached_vol = [leach_135mmh_fresh, leach_135mmh_aged,
                   leach_55mmhA_fresh, leach_55mmhA_aged,
                   leach_55mmhB_fresh, leach_55mmhB_aged,
                   leach_30mmh_fresh, leach_30mmh_aged]

    ponded_vol = [pond_135mmh_fresh, pond_135mmh_aged,
                  pond_55mmhA_fresh, pond_55mmhA_aged,
                  pond_55mmhB_fresh, pond_55mmhB_aged,
                  pond_30mmh_fresh, pond_30mmh_aged]

    # infil_135mmh = infil_data[:, 1]  # infil. volume per timestep @ high intensity
    # infil_55mmh = infil_data[:, 2]
    # infil_30mmh = infil_data[:, 3]
    # infil_vol = [infil_135mmh, infil_55mmh, infil_30mmh]  # cm3

    # define a max error
    error_fresh = 10 ** 9
    error_aged = 10 ** 9
    # Check best Kd fit from selected range
    for k in range(len(Kd)):
        K_L = kfilm(d, runoffvelocity)  # mm/min
        #  for index, (name, (mass_i, leach_obs, pond_obs)) in enumerate(pest_dict.iteritems()): # python 2.7
        for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) \
                in enumerate(sorted(pest_dict.items())):  # python 3.

            # Leaching
            cum_mass_out = 0
            mass_out_dt = []
            cum_mass_out_dt = []

            # Runoff
            mass_in_overflow = 0
            conc_in_overflow = 0
            mass_overflow_dt = []

            mass_tot = mass_i

            if index % 2 == 0:
                # Whelan et al. 1987 (p4.9)
                conc_soil_old = \
                    ((mass_tot / soil_m) * pb_i * Kd[k]/x / (ov + pb_i * (Kd[k]/x)))  # ug/g
                conc_liq_old = \
                    ((mass_tot / soil_m) * pb_i) / (ov + pb_i * Kd[k]/x)  # ug/g

                r_factor = 1 + (pb_i * Kd[k]/x) / ovSat
                conc_liq_old = conc_soil_old / (Kd[k]/x)  # ug/mm3

                #  Compute leaching for each time step under scenario of "index"
                hasleached = 0
                for t in range(len(cum_time_30min)):
                    # Update r_factor, based on increasing pb.
                    if leached_vol[index][t] > 0:
                        hasleached += 1
                    if hasleached > 12:
                        pb_i = pb_f
                        r_factor = 1 + (pb_i * Kd[k]/x) / ovSat

                    # McGrath leaching loss:
                    conc_liq_new = conc_liq_old * exp(-(leached_vol[index][t]) / (r_factor * vol_h2o_sat))
                    # ug/mm3 * exp(- mm3/[-]*mm3) = ug/mm3
                    mass_liq_new = conc_liq_new * vol_h2o_sat  # ug
                    mass_out = (conc_liq_old - conc_liq_new) * vol_h2o_sat
                    # ug
                    cum_mass_out += mass_out

                    # Store temporal and cumulative mass leached
                    mass_out_dt.append(mass_out)
                    cum_mass_out_dt.append(cum_mass_out)

                    # Overland flow exchange/transfer
                    if ponded_vol[index][t] > 0:
                        if KFILM:
                            # Mixing layer model from Havis (1986) and Havis et al. (1992)
                            # Kfilm equation from Bennett and Myers, 1982; in Shi et al., 2011, p.1220)
                            # mass_transfered_overflow = K_L * (conc_liq_new-conc_in_overflow) * area * time_sizes[index][t]
                            mass_transfered_overflow = K_L * conc_liq_new * area * time_sizes[index][t]
                            # mm/min * ug/mm3 * mm2 * dt = [ug]
                            if mass_transfered_overflow < 0:
                                if conc_in_overflow * ponded_vol[index][t] < mass_transfered_overflow:
                                    mass_transfered_overflow = conc_in_overflow * ponded_vol[index][t]
                                else:
                                    print ("mass in of: ", conc_in_overflow * ponded_vol[index][t])
                                    print ("mass -> pond: ", mass_transfered_overflow, "at: ", t)
                                    print ("conc liq (t-1): ", conc_liq_old)
                                    print ("conc liq (t): ", conc_liq_new)
                                    print("conc in of (t-1): ", mass_overflow_dt[t - 1] / vol_h2o_sat)
                                    print ("conc in of (t): ", conc_in_overflow)

                            mass_in_overflow += mass_transfered_overflow
                            conc_in_overflow = mass_in_overflow / ponded_vol[index][t]  # ug/mm3
                            mass_liq_new -= mass_transfered_overflow
                            conc_liq_new = mass_liq_new/vol_h2o_sat
                            # conc_liq_new += conc_in_overflow*infil_vol[index][t]/vol_h2o_sat
                            # mass_in_overflow -= conc_in_overflow*infil_vol[index][t]

                        #                        if mass_in_overflow < 0:
                        #                            print("Error in mixing layer approach, mass is negative")
                        else:
                            print("Error, over land flow transfer not defined")

                    # Store the temporal mass in overflow
                    mass_overflow_dt.append(mass_in_overflow)

                    # re-equilibrate
                    mass_tot_new = mass_tot - mass_out
                    # conc_soil_new = conc_liq_new * Kd[k]
                    # conc_soil_old = Kd[k]*conc_liq_new

                    conc_liq_old = conc_liq_new
                    mass_tot = mass_tot_new

            else:
                # Whelan et al. 1987 (p4.9)
                conc_soil_old = \
                    ((mass_tot / soil_m) * pb_i * Kd[k] / (ov + pb_i * (Kd[k])))  # ug/g
                conc_liq_old = \
                    ((mass_tot / soil_m) * pb_i) / (ov + pb_i * Kd[k])  # ug/g

                r_factor = 1 + (pb_i * Kd[k]) / ovSat
                conc_liq_old = conc_soil_old / (Kd[k])  # ug/mm3

                #  Compute leaching for each time step under scenario of "index"
                hasleached = 0
                for t in range(len(cum_time_30min)):
                    # Update r_factor, based on increasing pb.
                    if leached_vol[index][t] > 0:
                        hasleached += 1
                    if hasleached > 12:
                        pb_i = pb_f
                        r_factor = 1 + (pb_i * Kd[k]) / ovSat

                    # McGrath leaching loss:
                    conc_liq_new = conc_liq_old * exp(-(leached_vol[index][t]) / (r_factor * vol_h2o_sat))
                    # ug/mm3 * exp(- mm3/[-]*mm3) = ug/mm3
                    mass_liq_new = conc_liq_new * vol_h2o_sat  # ug
                    mass_out = (conc_liq_old - conc_liq_new) * vol_h2o_sat
                    # ug
                    cum_mass_out += mass_out

                    # Store temporal and cumulative mass leached
                    mass_out_dt.append(mass_out)
                    cum_mass_out_dt.append(cum_mass_out)

                    # Overland flow exchange/transfer
                    if ponded_vol[index][t] > 0:
                        if KFILM:
                            # Mixing layer model from Havis (1986) and Havis et al. (1992)
                            # Kfilm equation from Bennett and Myers, 1982; in Shi et al., 2011, p.1220)
                            mass_transfered_overflow = K_L * conc_liq_new * area * time_sizes[index][t]
                            # mm/min * ug/mm3 * mm2 * dt = [ug]
                            if mass_transfered_overflow < 0:
                                if conc_in_overflow * ponded_vol[index][t] < mass_transfered_overflow:
                                    mass_transfered_overflow = conc_in_overflow * ponded_vol[index][t]
                                else:
                                    print ("mass in of: ", conc_in_overflow * ponded_vol[index][t])
                                    print ("mass -> pond: ", mass_transfered_overflow, "at: ", t)
                                    print ("conc liq (t-1): ", conc_liq_old)
                                    print ("conc liq (t): ", conc_liq_new)
                                    print("conc in of (t-1): ", mass_overflow_dt[t - 1] / vol_h2o_sat)
                                    print ("conc in of (t): ", conc_in_overflow)

                            mass_in_overflow += mass_transfered_overflow
                            conc_in_overflow = mass_in_overflow / ponded_vol[index][t]  # ug/mm3
                            mass_liq_new -= mass_transfered_overflow
                            conc_liq_new = mass_liq_new / vol_h2o_sat
                        else:
                            print("Error, over land flow transfer not defined")

                    # Store the temporal mass in overflow
                    mass_overflow_dt.append(mass_in_overflow)
                    mass_tot_new = mass_tot - mass_out
                    conc_liq_old = conc_liq_new
                    mass_tot = mass_tot_new

            if index == 0:
                try:
                    high_fresh_error = (cum_mass_out_dt[5] - leach_obs) ** 2 + \
                                       (mass_overflow_dt[5] - pond_obs) ** 2
                    SS0 = (leach_obs - cum_mass_out_dt[5]) ** 2
                except TypeError:
                    try:
                        high_fresh_error = (cum_mass_out_dt[5] - leach_obs) ** 2
                        SS0 = (leach_obs - cum_mass_out_dt[5]) ** 2
                        # high_fresh_error_prc = ((cum_mass_out_dt[5] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            high_fresh_error = (cum_mass_out_dt[5] - pond_obs) ** 2
                            SS0 = 0.

                        except TypeError:
                            high_fresh_error = 0.
                            SS0 = 0.

                temp_out_high_0d = mass_out_dt
                temp_cum_out_high_0d = cum_mass_out_dt
                temp_overflow_high_0d = mass_overflow_dt

            elif index == 1:
                try:
                    high_aged_error = (cum_mass_out_dt[5] - leach_obs) ** 2 + \
                                    (mass_overflow_dt[5] - pond_obs) ** 2
                    SS1 = (leach_obs - cum_mass_out_dt[5]) ** 2

                except TypeError:
                    try:
                        high_aged_error = (cum_mass_out_dt[5] - leach_obs) ** 2
                        SS1 = (leach_obs - cum_mass_out_dt[5]) ** 2
                        # high_aged_error_prc = ((cum_mass_out_dt[5] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            high_aged_error = (mass_overflow_dt[5] - pond_obs) ** 2
                            SS1 = 0.
                            # high_aged_error_prc = 0.
                        except TypeError:
                            high_aged_error = 0.
                            SS1 = 0.
                            # high_aged_error_prc = 0.

                temp_out_high_1d = mass_out_dt
                temp_cum_out_high_1d = cum_mass_out_dt
                temp_overflow_high_1d = mass_overflow_dt

            elif index == 2:
                try:
                    med12_fresh_error = (cum_mass_out_dt[11] - leach_obs) ** 2 + \
                                        (mass_overflow_dt[11] - pond_obs) ** 2
                    SS2 = (leach_obs - cum_mass_out_dt[11]) ** 2

                except TypeError:
                    try:
                        med12_fresh_error = (cum_mass_out_dt[11] - leach_obs) ** 2
                        SS2 = (leach_obs - cum_mass_out_dt[11]) ** 2
                        # med12_fresh_error_prc = ((cum_mass_out_dt[11] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            med12_fresh_error = (mass_overflow_dt[11] - pond_obs) ** 2
                            SS2 = 0.
                            # med12_fresh_error_prc = 0
                        except TypeError:
                            med12_fresh_error = 0.
                            SS2 = 0.
                            # med12_fresh_error_prc = 0

                temp_out_med12_0d = mass_out_dt
                temp_cum_out_med12_0d = cum_mass_out_dt
                temp_overflow_med12_0d = mass_overflow_dt

            elif index == 3:
                try:
                    med12_aged__error = (cum_mass_out_dt[11] - leach_obs) ** 2 + \
                                        (mass_overflow_dt[11] - pond_obs) ** 2
                    SS3 = (leach_obs - cum_mass_out_dt[11]) ** 2

                except TypeError:
                    try:
                        med12_aged__error = (cum_mass_out_dt[11] - leach_obs) ** 2
                        SS3 = (leach_obs - cum_mass_out_dt[11]) ** 2
                        # med12_aged__error_prc = ((cum_mass_out_dt[11] - leach_obs) / leach_obs) * 100
                        # print("Index:", index)
                        # print("Pest dict: ", pest_dict['d_med12_1d'][1])
                        # print("SS3: ", cum_mass_out_dt[11], leach_obs)
                    except TypeError:
                        try:
                            med12_aged__error = (mass_overflow_dt[11] - pond_obs) ** 2
                            SS3 = 0.
                            # med12_aged__error_prc = 0.
                        except TypeError:
                            med12_aged__error = 0
                            SS3 = 0.
                            # med12_aged__error_prc = 0.

                temp_out_med12_1d = mass_out_dt
                temp_cum_out_med12_1d = cum_mass_out_dt
                temp_overflow_med12_1d = mass_overflow_dt

            elif index == 4:
                try:
                    med30_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                  (mass_overflow_dt[-1] - pond_obs) ** 2
                    SS4 = (leach_obs - cum_mass_out_dt[-1])**2
                except TypeError:
                    try:
                        med30_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        SS4 = (leach_obs - cum_mass_out_dt[-1])**2
                        # med30_fresh__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            med30_fresh__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                            SS4 = 0.
                            # med30_fresh__error_prc = 0.
                        except TypeError:
                            med30_fresh__error = 0
                            SS4 = 0.
                            # med30_fresh__error_prc = 0.

                temp_out_med30_0d = mass_out_dt
                temp_cum_out_med30_0d = cum_mass_out_dt
                temp_overflow_med30_0d = mass_overflow_dt

            elif index == 5:
                try:
                    med30_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                  (mass_overflow_dt[-1] - pond_obs) ** 2
                    SS5 = (leach_obs - cum_mass_out_dt[-1])**2

                except TypeError:
                    try:
                        med30_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        SS5 = (cum_mass_out_dt[-1] - leach_obs)**2
                        # med30_aged__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            med30_aged__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                            SS5 = 0.
                            # med30_aged__error_prc = 0.
                        except TypeError:
                            med30_aged__error = 0.
                            SS5 = 0.
                            # med30_aged__error_prc = 0.

                temp_out_med30_1d = mass_out_dt
                temp_cum_out_med30_1d = cum_mass_out_dt
                temp_overflow_med30_1d = mass_overflow_dt

            elif index == 6:
                try:
                    low_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                (mass_overflow_dt[-1] - pond_obs) ** 2
                    SS6 = (cum_mass_out_dt[-1] - leach_obs) ** 2

                except TypeError:
                    try:
                        low_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        SS6 = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        # low_fresh__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            low_fresh__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                            SS6 = 0.
                            # low_fresh__error_prc = 0.
                        except TypeError:
                            low_fresh__error = 0.
                            SS6 = 0.
                            # low_fresh__error_prc = 0.

                temp_out_low_0d = mass_out_dt
                temp_cum_out_low_0d = cum_mass_out_dt
                temp_overflow_low_0d = mass_overflow_dt

            elif index == 7:
                try:
                    low_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                (mass_overflow_dt[-1] - pond_obs) ** 2
                    SS7 = (cum_mass_out_dt[-1] - leach_obs) ** 2
                except TypeError:
                    try:
                        low_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        SS7 = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        # low_aged__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            low_aged__error = (mass_overflow_dt[-1] - pond_obs)
                            SS7 = 0.
                            # low_aged__error_prc = 0.
                        except TypeError:
                            low_aged__error = 0.
                            SS7 = 0.
                            # low_aged__error_prc = 0.

                temp_out_low_1d = mass_out_dt
                temp_cum_out_low_1d = cum_mass_out_dt
                temp_overflow_low_1d = mass_overflow_dt
            else:
                print("Index number error")

        SS_test = SS0 + SS1 + SS2 + SS3 + SS4 + SS5 + SS6 + SS7
        SS_fresh = SS0 + SS2 + SS4 + SS6
        SS_aged = SS1 + SS3 + SS5 + SS7

        if SS_fresh < error_fresh:
            error_fresh = SS_fresh
            SSres_fresh = SS0 + SS2 + SS4 + SS6
            r_squared_fresh = 1 - (SSres_fresh/SStot_fresh)
            kd_chosen_fresh = Kd[k]
            x_factor = x
            num_kd_fresh = k

            # a) Cumulative mass leached
            high_fresh_cum_mass_out_dt = temp_cum_out_high_0d
            try:
                high_fresh_error_prc = ((high_fresh_cum_mass_out_dt[5] - pest_dict['a_high_0d'][1]) / pest_dict['a_high_0d'][1]) * 100
            except TypeError:
                high_fresh_error_prc = "No obs | < LD "

            med12_fresh_cum_mass_out_dt = temp_cum_out_med12_0d
            try:
                med12_fresh_error_prc = ((med12_fresh_cum_mass_out_dt[11] - pest_dict['c_med12_0d'][1]) / pest_dict['c_med12_0d'][1]) * 100
            except TypeError:
                med12_fresh_error_prc = "No obs | < LD "

            med30_fresh_cum_mass_out_dt = temp_cum_out_med30_0d
            try:
                med30_fresh__error_prc = ((med30_fresh_cum_mass_out_dt[-1] - pest_dict['e_med30_0d'][1]) / pest_dict['e_med30_0d'][1]) * 100
            except TypeError:
                med30_fresh__error_prc = "No obs | < LD "

            low_fresh_cum_mas_out_dt = temp_cum_out_low_0d
            try:
                low_fresh__error_prc = ((low_fresh_cum_mas_out_dt[-1] - pest_dict['g_low_0d'][1]) / pest_dict['g_low_0d'][1]) * 100
            except TypeError:
                low_fresh__error_prc = "No obs | < LD "

            # b) Mass in ponded water
            high_fresh_overmass_dt = temp_overflow_high_0d
            med12_fresh_overmass_dt = temp_overflow_med12_0d
            med30_fresh_overmass_dt = temp_overflow_med30_0d
            low_fresh_overmass_dt = temp_overflow_low_0d

        if SS_aged < error_aged:
            error_aged = SS_aged
            SSres_aged = SS1 + SS3 + SS5 + SS7
            r_squared_aged = 1 - (SSres_aged/SStot_aged)
            kd_chosen_aged = Kd[k]
            num_kd_aged = k

            # a) Cumulative mass leached
            high_aged_cum_mass_out_dt = temp_cum_out_high_1d
            try:
                high_aged_error_prc = ((high_aged_cum_mass_out_dt[5] - pest_dict['b_high_1d'][1]) /
                                       pest_dict['b_high_1d'][1]) * 100
            except TypeError:
                high_aged_error_prc = "No obs | < LD "

            med12_aged_cum_mass_out_dt = temp_cum_out_med12_1d
            try:
                med12_aged__error_prc = ((med12_aged_cum_mass_out_dt[11] - pest_dict['d_med12_1d'][1]) /
                                         pest_dict['d_med12_1d'][1]) * 100
            except TypeError:
                med12_aged__error_prc = "No obs | < LD "

            med30_aged_cum_mass_out_dt = temp_cum_out_med30_1d
            try:
                med30_aged__error_prc = ((med30_aged_cum_mass_out_dt[-1] - pest_dict['f_med30_1d'][1]) /
                                         pest_dict['f_med30_1d'][1]) * 100
            except TypeError:
                med30_aged__error_prc = "No obs | < LD "

            low_aged_cum_mas_out_dt = temp_cum_out_low_1d
            try:
                low_aged__error_prc = (
                                      (low_aged_cum_mas_out_dt[-1] - pest_dict['h_low_1d'][1]) / pest_dict['h_low_1d'][1]) * 100
            except TypeError:
                low_aged__error_prc = "No obs | < LD "

            high_aged_overmass_dt = temp_overflow_high_1d
            med12_aged_overmass_dt = temp_overflow_med12_1d
            med30_aged_overmass_dt = temp_overflow_med30_1d
            low_aged_overmass_dt = temp_overflow_low_1d

        else:
            continue

    if isFirstCycle:
        print("1st Pulse")
    else:
        print("2nd Pulse")

    if isLiving:
        print("Living")
        fom_crop_untreat = 5.51 / 100.0
        foc = 0.58 * fom_crop_untreat
    else:
        print("Sterile")
        fom_crop_sterile = 3.87 / 100.0
        foc = 0.58 * fom_crop_sterile

    print("--------------------------------------------")
    print("Koc tested: ", np.array(Kd)/10**3/foc)  # Convert back to cm3/g
    print("Best log Kd (Fresh): ", log10(kd_chosen_fresh/10**3), "cm3/g", "( Num: ", num_kd_fresh + 1, ")", "\n",
          "x factor: ", x, "\n",
          "R2: ", r_squared_fresh)
    print("Best log Kd (Aged): ", log10(kd_chosen_aged / 10 ** 3), "cm3/g", "( Num: ", num_kd_aged + 1, ")", "\n",
          "x factor: No factor considered. ", "\n",
          "R2: ", r_squared_aged)
    print("Scenario - modality - Predicted error prcnt (%) | Predicted | Observed |")
    print("--------------------------------------------")
    print("(A) 135 mm/h - Fresh ", high_fresh_error_prc, "|", high_fresh_cum_mass_out_dt[5], "|", pest_dict['a_high_0d'][1])
    print("(A) 135 mm/h - Aged ", high_aged_error_prc, "|", high_aged_cum_mass_out_dt[5], "|", pest_dict['b_high_1d'][1])
    print("(B) 55 mm/h - Fresh ", med12_fresh_error_prc, "|", med12_fresh_cum_mass_out_dt[11], "|", pest_dict['c_med12_0d'][1])
    print("(B) 55 mm/h - Aged ", med12_aged__error_prc, "|", med12_aged_cum_mass_out_dt[11], "|", pest_dict['d_med12_1d'][1])
    print("(C) 55 mm/h - Fresh ", med30_fresh__error_prc, "|", med30_fresh_cum_mass_out_dt[-1], "|", pest_dict['e_med30_0d'][1])
    print("(C) 55 mm/h - Aged ", med30_aged__error_prc, "|", med12_aged_cum_mass_out_dt[-1], "|", pest_dict['f_med30_1d'][1])
    print("(D) 30 mm/h - Fresh ", low_fresh__error_prc, "|", high_fresh_cum_mass_out_dt[-1], "|", pest_dict['g_low_0d'][1])
    print("(D) 30 mm/h - Aged ", low_aged__error_prc, "|", high_aged_cum_mass_out_dt[-1], "|", pest_dict['h_low_1d'][1])

    # print ("SS1 to SS7: ", SS1, SS2, SS3, SS4, SS5, SS6, SS7)
    # print ("SSres: ", SSres)
    # print ("SStot: ", SStot)
    # print("mean mass: ", mean_mass, "sum_mass:", sum_mass)

    return stackdata16(cum_time_30min,
                       high_fresh_cum_mass_out_dt, high_aged_cum_mass_out_dt, # leached
                       med12_fresh_cum_mass_out_dt, med12_aged_cum_mass_out_dt,
                       med30_fresh_cum_mass_out_dt, med30_aged_cum_mass_out_dt,
                       low_fresh_cum_mas_out_dt, low_aged_cum_mas_out_dt,
                       high_fresh_overmass_dt, high_aged_overmass_dt,  # ponded
                       med12_fresh_overmass_dt, med12_aged_overmass_dt,
                       med30_fresh_overmass_dt, med30_aged_overmass_dt,
                       low_fresh_overmass_dt, low_aged_overmass_dt)


def pest_linear(
        Kd,  # Kd list range to test
        pest_dict,  # Initial mass and observed data
        pb_i, pb_f,
        percol_data, pond_data,
        time_sizes,
        area, soil_height,
        d, runoffvelocity,
        soil,
        isFirstCycle,
        isLiving,
        KFILM=True):
    """
    This function takes 8 hydrological time series,
    based optimal Ksat values.

    The dictionary is either Sterile or Living

    :param Kd: ensure units are mm3/g
    :param pb: ensure units are g/mm3
    :param pest_dict:
    :param ov_sat:
    :param percol_data: 8 time series
    :param pond_data: 8 time series
    :param time_sizes: 8 time series
    :param area:
    :param soil_height:
    :param d:
    :param runoffvelocity:
    :param KFILM:
    :return:
        16 time series of leached (0d and 10d) and ponded mass (0d and 10d).
        Sees no distinction between sterile or live
    """

    # Leaching, ponding and infiltration in vol (cm3)
    # are stored as lists to iterate through

    if isFirstCycle:
        if soil == 'Alteck':
            ovSat = 0.61
            ov = 0.20
        elif soil == 'Rouff':
            ovSat = 0.55
            ov = 0.20
    else:
        # Saturated water content (0.61 - Alteck; 0.55 - Rouff)
        if soil == 'Alteck':
            ovSat = 0.61
            ov = 0.61 - 0.039
        elif soil == 'Rouff':
            ovSat = 0.55
            ov = 0.55 - 0.038

    soil_m = 54.047  # [g] Average soil mass used in the experiment.

    vol_h2o_sat = area * soil_height * ovSat  # mm3
    cum_time_30min = pond_data[:, 0]

    min_mass = 200
    max_mass = 0
    sum_mass = 0
    masses = 0
    for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) \
            in enumerate(sorted(pest_dict.items())):
        # print (index, (name, (mass_i, leach_obs, pond_obs)))

        if float(leach_obs) > 0.:
            masses += 1
            sum_mass += float(leach_obs)
        if float(leach_obs) > max_mass:
            max_mass = float(leach_obs)
        if float(leach_obs) < min_mass:
            min_mass = float(leach_obs)

    mean_mass = sum_mass/masses

    SStot = 0
    for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) in enumerate(pest_dict.items()):
        if float(leach_obs) > 0:
            SStot += (float(leach_obs) - mean_mass)**2
        else:
            SStot += mean_mass**2

    if isFirstCycle:
        leach_135mmh_fresh = percol_data[:, 1]  # leached volume per timestep @ high intensity
        leach_55mmhA_fresh = percol_data[:, 2]
        leach_55mmhB_fresh = percol_data[:, 2]
        leach_30mmh_fresh = percol_data[:, 3]

        leach_135mmh_aged = percol_data[:, 1]  # leached volume per timestep @ high intensity
        leach_55mmhA_aged = percol_data[:, 2]
        leach_55mmhB_aged = percol_data[:, 2]
        leach_30mmh_aged = percol_data[:, 3]

        pond_135mmh_fresh = pond_data[:, 1]  # ponded volume per timestep @ high intensity
        pond_55mmhA_fresh = pond_data[:, 2]
        pond_55mmhB_fresh = pond_data[:, 2]
        pond_30mmh_fresh = pond_data[:, 3]

        pond_135mmh_aged = pond_data[:, 1]  # ponded volume per timestep @ high intensity
        pond_55mmhA_aged = pond_data[:, 2]
        pond_55mmhB_aged = pond_data[:, 2]
        pond_30mmh_aged = pond_data[:, 3]

    else:
        if isLiving:
            leach_135mmh_fresh = percol_data[:, 9]  # leached volume per timestep @ high intensity
            leach_55mmhA_fresh = percol_data[:, 10]
            leach_55mmhB_fresh = percol_data[:, 11]
            leach_30mmh_fresh = percol_data[:, 12]

            leach_135mmh_aged = percol_data[:, 13]  # leached volume per timestep @ high intensity
            leach_55mmhA_aged = percol_data[:, 14]
            leach_55mmhB_aged = percol_data[:, 15]
            leach_30mmh_aged = percol_data[:, 16]

            pond_135mmh_fresh = pond_data[:, 9]  # ponded volume per timestep @ high intensity
            pond_55mmhA_fresh = pond_data[:, 10]
            pond_55mmhB_fresh = pond_data[:, 11]
            pond_30mmh_fresh = pond_data[:, 12]

            pond_135mmh_aged = pond_data[:, 13]  # ponded volume per timestep @ high intensity
            pond_55mmhA_aged = pond_data[:, 14]
            pond_55mmhB_aged = pond_data[:, 15]
            pond_30mmh_aged = pond_data[:, 16]
        else:
            leach_135mmh_fresh = percol_data[:, 1]  # leached volume per timestep @ high intensity
            leach_55mmhA_fresh = percol_data[:, 2]
            leach_55mmhB_fresh = percol_data[:, 3]
            leach_30mmh_fresh = percol_data[:, 4]

            leach_135mmh_aged = percol_data[:, 5]  # leached volume per timestep @ high intensity
            leach_55mmhA_aged = percol_data[:, 6]
            leach_55mmhB_aged = percol_data[:, 7]
            leach_30mmh_aged = percol_data[:, 8]

            pond_135mmh_fresh = pond_data[:, 1]  # ponded volume per timestep @ high intensity
            pond_55mmhA_fresh = pond_data[:, 2]
            pond_55mmhB_fresh = pond_data[:, 3]
            pond_30mmh_fresh = pond_data[:, 4]

            pond_135mmh_aged = pond_data[:, 5]  # ponded volume per timestep @ high intensity
            pond_55mmhA_aged = pond_data[:, 6]
            pond_55mmhB_aged = pond_data[:, 7]
            pond_30mmh_aged = pond_data[:, 8]

    leached_vol = [leach_135mmh_fresh, leach_135mmh_aged,
                   leach_55mmhA_fresh, leach_55mmhA_aged,
                   leach_55mmhB_fresh, leach_55mmhB_aged,
                   leach_30mmh_fresh, leach_30mmh_aged]

    ponded_vol = [pond_135mmh_fresh, pond_135mmh_aged,
                  pond_55mmhA_fresh, pond_55mmhA_aged,
                  pond_55mmhB_fresh, pond_55mmhB_aged,
                  pond_30mmh_fresh, pond_30mmh_aged]

    # infil_135mmh = infil_data[:, 1]  # infil. volume per timestep @ high intensity
    # infil_55mmh = infil_data[:, 2]
    # infil_30mmh = infil_data[:, 3]
    # infil_vol = [infil_135mmh, infil_55mmh, infil_30mmh]  # cm3

    # Check best Kd fit from selected range
    error = 10 ** 9  # define a max error
    for k in range(len(Kd)):
        K_L = kfilm(d, runoffvelocity)  # mm/min
        #  for index, (name, (mass_i, leach_obs, pond_obs)) in enumerate(pest_dict.iteritems()): # python 2.7
        for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) \
                in enumerate(sorted(pest_dict.items())):  # python 3.

            # Leaching
            cum_mass_out = 0
            mass_out_dt = []
            cum_mass_out_dt = []

            # Runoff
            mass_in_overflow = 0
            conc_in_overflow = 0
            mass_overflow_dt = []

            mass_tot = mass_i

            # Whelan et al. 1987 (p4.9)
            conc_soil_old = \
                ((mass_tot / soil_m) * pb_i * Kd[k]) / (ov + pb_i * Kd[k])  # ug/g
            conc_liq_old = \
                ((mass_tot / soil_m) * pb_i) / (ov + pb_i * Kd[k])  # ug/g

            r_factor = 1 + (pb_i * Kd[k]) / ovSat
            conc_liq_old = conc_soil_old / Kd[k]  # ug/mm3

            #  Compute leaching for each time step under scenario of "index"
            hasleached = 0
            for t in range(len(cum_time_30min)):
                # Update r_factor, based on increasing pb.
                if leached_vol[index][t] > 0:
                    hasleached += 1
                if hasleached > 12:
                    pb_i = pb_f
                    r_factor = 1 + (pb_i * Kd[k]) / ovSat

                # McGrath leaching loss:
                conc_liq_new = conc_liq_old * exp(-(leached_vol[index][t]) / (r_factor * vol_h2o_sat))
                # ug/mm3 * exp(- mm3/[-]*mm3) = ug/mm3
                mass_liq_new = conc_liq_new * vol_h2o_sat  # ug
                mass_out = (conc_liq_old - conc_liq_new) * vol_h2o_sat
                # ug
                cum_mass_out += mass_out

                # Store temporal and cumulative mass leached
                mass_out_dt.append(mass_out)
                cum_mass_out_dt.append(cum_mass_out)

                # Overland flow exchange/transfer
                if ponded_vol[index][t] > 0:
                    if KFILM:
                        # Mixing layer model from Havis (1986) and Havis et al. (1992)
                        # Kfilm equation from Bennett and Myers, 1982; in Shi et al., 2011, p.1220)
                        # mass_transfered_overflow = K_L * (conc_liq_new-conc_in_overflow) * area * time_sizes[index][t]
                        mass_transfered_overflow = K_L * conc_liq_new * area * time_sizes[index][t]
                        # mm/min * ug/mm3 * mm2 * dt = [ug]
                        if mass_transfered_overflow < 0:
                            if conc_in_overflow * ponded_vol[index][t] < mass_transfered_overflow:
                                mass_transfered_overflow = conc_in_overflow * ponded_vol[index][t]
                            else:
                                print ("mass in of: ", conc_in_overflow * ponded_vol[index][t])
                                print ("mass -> pond: ", mass_transfered_overflow, "at: ", t)
                                print ("conc liq (t-1): ", conc_liq_old)
                                print ("conc liq (t): ", conc_liq_new)
                                print("conc in of (t-1): ", mass_overflow_dt[t - 1] / vol_h2o_sat)
                                print ("conc in of (t): ", conc_in_overflow)

                        mass_in_overflow += mass_transfered_overflow
                        conc_in_overflow = mass_in_overflow / ponded_vol[index][t]  # ug/mm3
                        mass_liq_new -= mass_transfered_overflow
                        conc_liq_new = mass_liq_new/vol_h2o_sat
                        # conc_liq_new += conc_in_overflow*infil_vol[index][t]/vol_h2o_sat
                        # mass_in_overflow -= conc_in_overflow*infil_vol[index][t]

                    #                        if mass_in_overflow < 0:
                    #                            print("Error in mixing layer approach, mass is negative")
                    else:
                        print("Error, over land flow transfer not defined")

                # Store the temporal mass in overflow
                mass_overflow_dt.append(mass_in_overflow)

                # re-equilibrate
                mass_tot_new = mass_tot - mass_out
                # conc_soil_new = conc_liq_new * Kd[k]
                # conc_soil_old = Kd[k]*conc_liq_new

                conc_liq_old = conc_liq_new
                mass_tot = mass_tot_new

            if index == 0:
                try:
                    high_fresh_error = (cum_mass_out_dt[5] - leach_obs) ** 2 + \
                                       (mass_overflow_dt[5] - pond_obs) ** 2
                    SS0 = (leach_obs - cum_mass_out_dt[5]) ** 2
                except TypeError:
                    try:
                        high_fresh_error = (cum_mass_out_dt[5] - leach_obs) ** 2
                        SS0 = (leach_obs - cum_mass_out_dt[5]) ** 2
                        # high_fresh_error_prc = ((cum_mass_out_dt[5] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            high_fresh_error = (cum_mass_out_dt[5] - pond_obs) ** 2
                            SS0 = 0.

                        except TypeError:
                            high_fresh_error = 0.
                            SS0 = 0.


                temp_out_high_0d = mass_out_dt
                temp_cum_out_high_0d = cum_mass_out_dt
                temp_overflow_high_0d = mass_overflow_dt

            elif index == 1:
                try:
                    high_aged_error = (cum_mass_out_dt[5] - leach_obs) ** 2 + \
                                    (mass_overflow_dt[5] - pond_obs) ** 2
                    SS1 = (leach_obs - cum_mass_out_dt[5]) ** 2

                except TypeError:
                    try:
                        high_aged_error = (cum_mass_out_dt[5] - leach_obs) ** 2
                        SS1 = (leach_obs - cum_mass_out_dt[5]) ** 2
                        # high_aged_error_prc = ((cum_mass_out_dt[5] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            high_aged_error = (mass_overflow_dt[5] - pond_obs) ** 2
                            SS1 = 0.
                            # high_aged_error_prc = 0.
                        except TypeError:
                            high_aged_error = 0.
                            SS1 = 0.
                            # high_aged_error_prc = 0.

                temp_out_high_1d = mass_out_dt
                temp_cum_out_high_1d = cum_mass_out_dt
                temp_overflow_high_1d = mass_overflow_dt

            elif index == 2:
                try:
                    med12_fresh_error = (cum_mass_out_dt[11] - leach_obs) ** 2 + \
                                        (mass_overflow_dt[11] - pond_obs) ** 2
                    SS2 = (leach_obs - cum_mass_out_dt[11]) ** 2

                except TypeError:
                    try:
                        med12_fresh_error = (cum_mass_out_dt[11] - leach_obs) ** 2
                        SS2 = (leach_obs - cum_mass_out_dt[11]) ** 2
                        # med12_fresh_error_prc = ((cum_mass_out_dt[11] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            med12_fresh_error = (mass_overflow_dt[11] - pond_obs) ** 2
                            SS2 = 0.
                            # med12_fresh_error_prc = 0
                        except TypeError:
                            med12_fresh_error = 0.
                            SS2 = 0.
                            # med12_fresh_error_prc = 0

                temp_out_med12_0d = mass_out_dt
                temp_cum_out_med12_0d = cum_mass_out_dt
                temp_overflow_med12_0d = mass_overflow_dt

            elif index == 3:
                try:
                    med12_aged__error = (cum_mass_out_dt[11] - leach_obs) ** 2 + \
                                        (mass_overflow_dt[11] - pond_obs) ** 2
                    SS3 = (leach_obs - cum_mass_out_dt[11]) ** 2

                except TypeError:
                    try:
                        med12_aged__error = (cum_mass_out_dt[11] - leach_obs) ** 2
                        SS3 = (leach_obs - cum_mass_out_dt[11]) ** 2
                        # med12_aged__error_prc = ((cum_mass_out_dt[11] - leach_obs) / leach_obs) * 100
                        # print("Index:", index)
                        # print("Pest dict: ", pest_dict['d_med12_1d'][1])
                        # print("SS3: ", cum_mass_out_dt[11], leach_obs)
                    except TypeError:
                        try:
                            med12_aged__error = (mass_overflow_dt[11] - pond_obs) ** 2
                            SS3 = 0.
                            # med12_aged__error_prc = 0.
                        except TypeError:
                            med12_aged__error = 0
                            SS3 = 0.
                            # med12_aged__error_prc = 0.

                temp_out_med12_1d = mass_out_dt
                temp_cum_out_med12_1d = cum_mass_out_dt
                temp_overflow_med12_1d = mass_overflow_dt

            elif index == 4:
                try:
                    med30_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                  (mass_overflow_dt[-1] - pond_obs) ** 2
                    SS4 = (leach_obs - cum_mass_out_dt[-1])**2
                except TypeError:
                    try:
                        med30_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        SS4 = (leach_obs - cum_mass_out_dt[-1])**2
                        # med30_fresh__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            med30_fresh__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                            SS4 = 0.
                            # med30_fresh__error_prc = 0.
                        except TypeError:
                            med30_fresh__error = 0
                            SS4 = 0.
                            # med30_fresh__error_prc = 0.

                temp_out_med30_0d = mass_out_dt
                temp_cum_out_med30_0d = cum_mass_out_dt
                temp_overflow_med30_0d = mass_overflow_dt

            elif index == 5:
                try:
                    med30_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                  (mass_overflow_dt[-1] - pond_obs) ** 2
                    SS5 = (leach_obs - cum_mass_out_dt[-1])**2

                except TypeError:
                    try:
                        med30_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        SS5 = (cum_mass_out_dt[-1] - leach_obs)**2
                        # med30_aged__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            med30_aged__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                            SS5 = 0.
                            # med30_aged__error_prc = 0.
                        except TypeError:
                            med30_aged__error = 0.
                            SS5 = 0.
                            # med30_aged__error_prc = 0.

                temp_out_med30_1d = mass_out_dt
                temp_cum_out_med30_1d = cum_mass_out_dt
                temp_overflow_med30_1d = mass_overflow_dt

            elif index == 6:
                try:
                    low_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                (mass_overflow_dt[-1] - pond_obs) ** 2
                    SS6 = (cum_mass_out_dt[-1] - leach_obs) ** 2

                except TypeError:
                    try:
                        low_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        SS6 = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        # low_fresh__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            low_fresh__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                            SS6 = 0.
                            # low_fresh__error_prc = 0.
                        except TypeError:
                            low_fresh__error = 0.
                            SS6 = 0.
                            # low_fresh__error_prc = 0.

                temp_out_low_0d = mass_out_dt
                temp_cum_out_low_0d = cum_mass_out_dt
                temp_overflow_low_0d = mass_overflow_dt
            elif index == 7:
                try:
                    low_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                (mass_overflow_dt[-1] - pond_obs) ** 2
                    SS7 = (cum_mass_out_dt[-1] - leach_obs) ** 2
                except TypeError:
                    try:
                        low_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        SS7 = (cum_mass_out_dt[-1] - leach_obs) ** 2
                        # low_aged__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            low_aged__error = (mass_overflow_dt[-1] - pond_obs)
                            SS7 = 0.
                            # low_aged__error_prc = 0.
                        except TypeError:
                            low_aged__error = 0.
                            SS7 = 0.
                            # low_aged__error_prc = 0.

                temp_out_low_1d = mass_out_dt
                temp_cum_out_low_1d = cum_mass_out_dt
                temp_overflow_low_1d = mass_overflow_dt
            else:
                print("Index number error")

        error_test = (((high_fresh_error + high_aged_error +
                        med12_fresh_error + med12_aged__error +
                        med30_fresh__error + med30_aged__error +
                        low_fresh__error + low_aged__error) / 8.0) ** 0.5)/(max_mass - min_mass)

        SS_test = SS0 + SS1 + SS2 + SS3 + SS4 + SS5 + SS6 + SS7

        # if error_test < error:
        if SS_test < error:
            # error = error_test
            error = SS_test
            SSres = SS0 + SS1 + SS2 + SS3 + SS4 + SS5 + SS6 + SS7
            r_squared = 1 - (SSres/SStot)
            kd_chosen = Kd[k]
            num_kd = k
            # print("error = ", error,
            #      "log Kd = ", log10(kd_chosen), "Num: ", num_kd+1)
            # Store final variables

            # a) Cumulative mass leached
            high_fresh_cum_mass_out_dt = temp_cum_out_high_0d
            try:
                high_fresh_error_prc = ((high_fresh_cum_mass_out_dt[5] - pest_dict['a_high_0d'][1]) / pest_dict['a_high_0d'][1]) * 100
            except TypeError:
                high_fresh_error_prc = "No obs | < LD "

            high_aged_cum_mass_out_dt = temp_cum_out_high_1d
            try:
                high_aged_error_prc = ((high_aged_cum_mass_out_dt[5] - pest_dict['b_high_1d'][1]) / pest_dict['b_high_1d'][1]) * 100
            except TypeError:
                high_aged_error_prc = "No obs | < LD "

            med12_fresh_cum_mass_out_dt = temp_cum_out_med12_0d
            try:
                med12_fresh_error_prc = ((med12_fresh_cum_mass_out_dt[11] - pest_dict['c_med12_0d'][1]) / pest_dict['c_med12_0d'][1]) * 100
            except TypeError:
                med12_fresh_error_prc = "No obs | < LD "

            med12_aged_cum_mass_out_dt = temp_cum_out_med12_1d
            try:
                med12_aged__error_prc = ((med12_aged_cum_mass_out_dt[11] - pest_dict['d_med12_1d'][1]) / pest_dict['d_med12_1d'][1]) * 100
            except TypeError:
                med12_aged__error_prc = "No obs | < LD "

            med30_fresh_cum_mass_out_dt = temp_cum_out_med30_0d
            try:
                med30_fresh__error_prc = ((med30_fresh_cum_mass_out_dt[-1] - pest_dict['e_med30_0d'][1]) / pest_dict['e_med30_0d'][1]) * 100
            except TypeError:
                med30_fresh__error_prc = "No obs | < LD "

            med30_aged_cum_mass_out_dt = temp_cum_out_med30_1d
            try:
                med30_aged__error_prc = ((med30_aged_cum_mass_out_dt[-1] - pest_dict['f_med30_1d'][1]) / pest_dict['f_med30_1d'][1]) * 100
            except TypeError:
                med30_aged__error_prc = "No obs | < LD "

            low_fresh_cum_mas_out_dt = temp_cum_out_low_0d
            try:
                low_fresh__error_prc = ((low_fresh_cum_mas_out_dt[-1] - pest_dict['g_low_0d'][1]) / pest_dict['g_low_0d'][1]) * 100
            except TypeError:
                low_fresh__error_prc = "No obs | < LD "

            low_aged_cum_mas_out_dt = temp_cum_out_low_1d
            try:
                low_aged__error_prc = ((low_aged_cum_mas_out_dt[-1] - pest_dict['h_low_1d'][1]) / pest_dict['h_low_1d'][1]) * 100
            except TypeError:
                low_aged__error_prc = "No obs | < LD "

            # b) Mass in ponded water
            high_fresh_overmass_dt = temp_overflow_high_0d
            high_aged_overmass_dt = temp_overflow_high_1d
            med12_fresh_overmass_dt = temp_overflow_med12_0d
            med12_aged_overmass_dt = temp_overflow_med12_1d
            med30_fresh_overmass_dt = temp_overflow_med30_0d
            med30_aged_overmass_dt = temp_overflow_med30_1d
            low_fresh_overmass_dt = temp_overflow_low_0d
            low_aged_overmass_dt = temp_overflow_low_1d

        else:
            continue

    if isFirstCycle:
        print("1st Pulse")
    else:
        print("2nd Pulse")

    if isLiving:
        print("Living")
        fom_crop_untreat = 5.51 / 100.0
        foc = 0.58 * fom_crop_untreat
    else:
        print("Sterile")
        fom_crop_sterile = 3.87 / 100.0
        foc = 0.58 * fom_crop_sterile

    print("--------------------------------------------")
    print("Koc tested: ", np.array(Kd)/10**3/foc)  # Convert back to cm3/g
    print("Best log Kd: ", log10(kd_chosen/10**3), "cm3/g", "( Num: ", num_kd + 1, ")", "\n",
          "R2: ", r_squared)
    print("Scenario - modality - Predicted error prcnt (%) | Predicted | Observed |")
    print("--------------------------------------------")
    print("(A) 135 mm/h - Fresh ", high_fresh_error_prc, "|", high_fresh_cum_mass_out_dt[5], "|", pest_dict['a_high_0d'][1])
    print("(A) 135 mm/h - Aged ", high_aged_error_prc, "|", high_aged_cum_mass_out_dt[5], "|", pest_dict['b_high_1d'][1])
    print("(B) 55 mm/h - Fresh ", med12_fresh_error_prc, "|", med12_fresh_cum_mass_out_dt[11], "|", pest_dict['c_med12_0d'][1])
    print("(B) 55 mm/h - Aged ", med12_aged__error_prc, "|", med12_aged_cum_mass_out_dt[11], "|", pest_dict['d_med12_1d'][1])
    print("(C) 55 mm/h - Fresh ", med30_fresh__error_prc, "|", med30_fresh_cum_mass_out_dt[-1], "|", pest_dict['e_med30_0d'][1])
    print("(C) 55 mm/h - Aged ", med30_aged__error_prc, "|", med12_aged_cum_mass_out_dt[-1], "|", pest_dict['f_med30_1d'][1])
    print("(D) 30 mm/h - Fresh ", low_fresh__error_prc, "|", high_fresh_cum_mass_out_dt[-1], "|", pest_dict['g_low_0d'][1])
    print("(D) 30 mm/h - Aged ", low_aged__error_prc, "|", high_aged_cum_mass_out_dt[-1], "|", pest_dict['h_low_1d'][1])

    # print ("SS1 to SS7: ", SS1, SS2, SS3, SS4, SS5, SS6, SS7)
    # print ("SSres: ", SSres)
    # print ("SStot: ", SStot)
    # print("mean mass: ", mean_mass, "sum_mass:", sum_mass)

    return stackdata16(cum_time_30min,
                       high_fresh_cum_mass_out_dt, high_aged_cum_mass_out_dt, # leached
                       med12_fresh_cum_mass_out_dt, med12_aged_cum_mass_out_dt,
                       med30_fresh_cum_mass_out_dt, med30_aged_cum_mass_out_dt,
                       low_fresh_cum_mas_out_dt, low_aged_cum_mas_out_dt,
                       high_fresh_overmass_dt, high_aged_overmass_dt,  # ponded
                       med12_fresh_overmass_dt, med12_aged_overmass_dt,
                       med30_fresh_overmass_dt, med30_aged_overmass_dt,
                       low_fresh_overmass_dt, low_aged_overmass_dt)


def pest_langmuir(
        Kd,  # Kd list range to test
        Qrange,  # a list of fractions to obtain Q (i.e. Cmax)
        pest_dict,  # Initial mass and observed data
        pb_i, pb_f,
        percol_data, pond_data,
        time_sizes,
        area, soil_height,
        d, runoffvelocity,
        soil,
        isFirstCycle,
        isLiving,
        KFILM=True):
    """
    This function takes 8 hydrological time series,
    based optimal Ksat values.

    The dictionary is either Sterile or Living

    :param Kd: ensure units are mm3/g
    :param pb: ensure units are g/mm3
    :param pest_dict:
    :param ov_sat:
    :param percol_data: 8 time series
    :param pond_data: 8 time series
    :param time_sizes: 8 time series
    :param area:
    :param soil_height:
    :param d:
    :param runoffvelocity:
    :param KFILM:
    :return:
        16 time series of leached (0d and 10d) and ponded mass (0d and 10d).
        Sees no distinction between sterile or live
    """

    # Leaching, ponding and infiltration in vol (cm3)
    # are stored as lists to iterate through

    if isFirstCycle:
        if soil == 'Alteck':
            ovSat = 0.61
            ov = 0.20
        elif soil == 'Rouff':
            ovSat = 0.55
            ov = 0.20
    else:
        # Saturated water content (0.61 - Alteck; 0.55 - Rouff)
        if soil == 'Alteck':
            ovSat = 0.61
            ov = 0.61 - 0.039
        elif soil == 'Rouff':
            ovSat = 0.55
            ov = 0.55 - 0.038

    soil_m = 54.047  # [g] Average soil mass used in the experiment.

    vol_h2o_sat = area * soil_height * ovSat  # mm3
    cum_time_30min = pond_data[:, 0]

    min_mass = 200
    max_mass = 0
    sum_mass = 0
    masses = 0
    for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) \
            in enumerate(sorted(pest_dict.items())):
        # print (index, (name, (mass_i, leach_obs, pond_obs)))

        if float(leach_obs) > 0.:
            masses += 1
            sum_mass += float(leach_obs)
        if float(leach_obs) > max_mass:
            max_mass = float(leach_obs)
        if float(leach_obs) < min_mass:
            min_mass = float(leach_obs)

    mean_mass = sum_mass/masses

    SStot = 0
    for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) \
            in enumerate(pest_dict.items()):
        if float(leach_obs) > 0:
            SStot += (float(leach_obs) - mean_mass)**2
        else:
            SStot += mean_mass**2

    if isFirstCycle:
        leach_135mmh_fresh = percol_data[:, 1]  # leached volume per timestep @ high intensity
        leach_55mmhA_fresh = percol_data[:, 2]
        leach_55mmhB_fresh = percol_data[:, 2]
        leach_30mmh_fresh = percol_data[:, 3]

        leach_135mmh_aged = percol_data[:, 1]  # leached volume per timestep @ high intensity
        leach_55mmhA_aged = percol_data[:, 2]
        leach_55mmhB_aged = percol_data[:, 2]
        leach_30mmh_aged = percol_data[:, 3]

        pond_135mmh_fresh = pond_data[:, 1]  # ponded volume per timestep @ high intensity
        pond_55mmhA_fresh = pond_data[:, 2]
        pond_55mmhB_fresh = pond_data[:, 2]
        pond_30mmh_fresh = pond_data[:, 3]

        pond_135mmh_aged = pond_data[:, 1]  # ponded volume per timestep @ high intensity
        pond_55mmhA_aged = pond_data[:, 2]
        pond_55mmhB_aged = pond_data[:, 2]
        pond_30mmh_aged = pond_data[:, 3]

    else:
        if isLiving:
            leach_135mmh_fresh = percol_data[:, 9]  # leached volume per timestep @ high intensity
            leach_55mmhA_fresh = percol_data[:, 10]
            leach_55mmhB_fresh = percol_data[:, 11]
            leach_30mmh_fresh = percol_data[:, 12]

            leach_135mmh_aged = percol_data[:, 13]  # leached volume per timestep @ high intensity
            leach_55mmhA_aged = percol_data[:, 14]
            leach_55mmhB_aged = percol_data[:, 15]
            leach_30mmh_aged = percol_data[:, 16]

            pond_135mmh_fresh = pond_data[:, 9]  # ponded volume per timestep @ high intensity
            pond_55mmhA_fresh = pond_data[:, 10]
            pond_55mmhB_fresh = pond_data[:, 11]
            pond_30mmh_fresh = pond_data[:, 12]

            pond_135mmh_aged = pond_data[:, 13]  # ponded volume per timestep @ high intensity
            pond_55mmhA_aged = pond_data[:, 14]
            pond_55mmhB_aged = pond_data[:, 15]
            pond_30mmh_aged = pond_data[:, 16]
        else:
            leach_135mmh_fresh = percol_data[:, 1]  # leached volume per timestep @ high intensity
            leach_55mmhA_fresh = percol_data[:, 2]
            leach_55mmhB_fresh = percol_data[:, 3]
            leach_30mmh_fresh = percol_data[:, 4]

            leach_135mmh_aged = percol_data[:, 5]  # leached volume per timestep @ high intensity
            leach_55mmhA_aged = percol_data[:, 6]
            leach_55mmhB_aged = percol_data[:, 7]
            leach_30mmh_aged = percol_data[:, 8]

            pond_135mmh_fresh = pond_data[:, 1]  # ponded volume per timestep @ high intensity
            pond_55mmhA_fresh = pond_data[:, 2]
            pond_55mmhB_fresh = pond_data[:, 3]
            pond_30mmh_fresh = pond_data[:, 4]

            pond_135mmh_aged = pond_data[:, 5]  # ponded volume per timestep @ high intensity
            pond_55mmhA_aged = pond_data[:, 6]
            pond_55mmhB_aged = pond_data[:, 7]
            pond_30mmh_aged = pond_data[:, 8]

    leached_vol = [leach_135mmh_fresh, leach_135mmh_aged,
                   leach_55mmhA_fresh, leach_55mmhA_aged,
                   leach_55mmhB_fresh, leach_55mmhB_aged,
                   leach_30mmh_fresh, leach_30mmh_aged]

    ponded_vol = [pond_135mmh_fresh, pond_135mmh_aged,
                  pond_55mmhA_fresh, pond_55mmhA_aged,
                  pond_55mmhB_fresh, pond_55mmhB_aged,
                  pond_30mmh_fresh, pond_30mmh_aged]

    # infil_135mmh = infil_data[:, 1]  # infil. volume per timestep @ high intensity
    # infil_55mmh = infil_data[:, 2]
    # infil_30mmh = infil_data[:, 3]
    # infil_vol = [infil_135mmh, infil_55mmh, infil_30mmh]  # cm3

    error = 10 ** 9  # define a max error
    # Check best Qmax fit from selected range
    for Q in Qrange:
        # Check best Kd fit from selected range
        for k in range(len(Kd)):
            K_L = kfilm(d, runoffvelocity)  # mm/min
            #  for index, (name, (mass_i, leach_obs, pond_obs)) in enumerate(pest_dict.iteritems()): # python 2.7
            for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) \
                    in enumerate(sorted(pest_dict.items())):  # python 3.

                # Leaching
                cum_mass_out = 0
                mass_out_dt = []
                cum_mass_out_dt = []

                # Runoff
                mass_in_overflow = 0
                conc_in_overflow = 0
                mass_overflow_dt = []

                mass_tot = mass_i

                # Whelan et al. 1987 (p4.9)
                conc_soil_old = \
                    ((mass_tot / soil_m) * pb_i * Kd[k]) / (ov + pb_i * Kd[k])  # ug/g
                """
                conc_liq_old = \
                    ((mass_tot / soil_m) * pb_i) / (ov + pb_i * Kd[k])  # ug/g
                """

                conc_liq_old = conc_soil_old / (Kd[k] * (Q - conc_soil_old))  # ug/mm3
                r_factor = 1 + (pb_i / ovSat) * (Kd[k] * Q * conc_liq_old) / (1 + Kd[k] * conc_liq_old) ** 2

                #  Compute leaching for each time step under scenario of "index"
                hasleached = 0
                for t in range(len(cum_time_30min)):
                    # Update r_factor, based on increasing pb.
                    if leached_vol[index][t] > 0:
                        hasleached += 1
                    if hasleached > 12:
                        pb_i = pb_f
                        r_factor = 1 + (pb_i / ovSat) * (Kd[k] * Q * conc_liq_old) / (1 + Kd[k] * conc_liq_old) ** 2

                    # McGrath leaching loss:
                    conc_liq_new = conc_liq_old * exp(-(leached_vol[index][t]) / (r_factor * vol_h2o_sat))
                    # ug/mm3 * exp(- mm3/[-]*mm3) = ug/mm3
                    mass_liq_new = conc_liq_new * vol_h2o_sat  # ug
                    mass_out = (conc_liq_old - conc_liq_new) * vol_h2o_sat
                    # ug
                    cum_mass_out += mass_out

                    # Store temporal and cumulative mass leached
                    mass_out_dt.append(mass_out)
                    cum_mass_out_dt.append(cum_mass_out)

                    # Overland flow exchange/transfer
                    if ponded_vol[index][t] > 0:
                        if KFILM:
                            # Mixing layer model from Havis (1986) and Havis et al. (1992)
                            # Kfilm equation from Bennett and Myers, 1982; in Shi et al., 2011, p.1220)
                            # mass_transfered_overflow = K_L * (conc_liq_new-conc_in_overflow) * area * time_sizes[index][t]
                            mass_transfered_overflow = K_L * conc_liq_new * area * time_sizes[index][t]
                            # mm/min * ug/mm3 * mm2 * dt = [ug]
                            if mass_transfered_overflow < 0:
                                if conc_in_overflow * ponded_vol[index][t] < mass_transfered_overflow:
                                    mass_transfered_overflow = conc_in_overflow * ponded_vol[index][t]
                                else:
                                    print ("mass in of: ", conc_in_overflow * ponded_vol[index][t])
                                    print ("mass -> pond: ", mass_transfered_overflow, "at: ", t)
                                    print ("conc liq (t-1): ", conc_liq_old)
                                    print ("conc liq (t): ", conc_liq_new)
                                    print("conc in of (t-1): ", mass_overflow_dt[t - 1] / vol_h2o_sat)
                                    print ("conc in of (t): ", conc_in_overflow)

                            mass_in_overflow += mass_transfered_overflow
                            conc_in_overflow = mass_in_overflow / ponded_vol[index][t]  # ug/mm3
                            mass_liq_new -= mass_transfered_overflow
                            conc_liq_new = mass_liq_new/vol_h2o_sat
                            # conc_liq_new += conc_in_overflow*infil_vol[index][t]/vol_h2o_sat
                            # mass_in_overflow -= conc_in_overflow*infil_vol[index][t]

                        #                        if mass_in_overflow < 0:
                        #                            print("Error in mixing layer approach, mass is negative")
                        else:
                            print("Error, over land flow transfer not defined")

                    # Store the temporal mass in overflow
                    mass_overflow_dt.append(mass_in_overflow)

                    # re-equilibrate
                    mass_tot_new = mass_tot - mass_out
                    conc_liq_old = conc_liq_new
                    mass_tot = mass_tot_new

                if index == 0:
                    try:
                        high_fresh_error = (cum_mass_out_dt[5] - leach_obs) ** 2 + \
                                           (mass_overflow_dt[5] - pond_obs) ** 2
                        SS0 = (leach_obs - cum_mass_out_dt[5]) ** 2
                    except TypeError:
                        try:
                            high_fresh_error = (cum_mass_out_dt[5] - leach_obs) ** 2
                            SS0 = (leach_obs - cum_mass_out_dt[5]) ** 2
                            # high_fresh_error_prc = ((cum_mass_out_dt[5] - leach_obs) / leach_obs) * 100
                        except TypeError:
                            try:
                                high_fresh_error = (cum_mass_out_dt[5] - pond_obs) ** 2
                                SS0 = 0.

                            except TypeError:
                                high_fresh_error = 0.
                                SS0 = 0.


                    temp_out_high_0d = mass_out_dt
                    temp_cum_out_high_0d = cum_mass_out_dt
                    temp_overflow_high_0d = mass_overflow_dt

                elif index == 1:
                    try:
                        high_aged_error = (cum_mass_out_dt[5] - leach_obs) ** 2 + \
                                        (mass_overflow_dt[5] - pond_obs) ** 2
                        SS1 = (leach_obs - cum_mass_out_dt[5]) ** 2

                    except TypeError:
                        try:
                            high_aged_error = (cum_mass_out_dt[5] - leach_obs) ** 2
                            SS1 = (leach_obs - cum_mass_out_dt[5]) ** 2
                            # high_aged_error_prc = ((cum_mass_out_dt[5] - leach_obs) / leach_obs) * 100
                        except TypeError:
                            try:
                                high_aged_error = (mass_overflow_dt[5] - pond_obs) ** 2
                                SS1 = 0.
                                # high_aged_error_prc = 0.
                            except TypeError:
                                high_aged_error = 0.
                                SS1 = 0.
                                # high_aged_error_prc = 0.

                    temp_out_high_1d = mass_out_dt
                    temp_cum_out_high_1d = cum_mass_out_dt
                    temp_overflow_high_1d = mass_overflow_dt

                elif index == 2:
                    try:
                        med12_fresh_error = (cum_mass_out_dt[11] - leach_obs) ** 2 + \
                                            (mass_overflow_dt[11] - pond_obs) ** 2
                        SS2 = (leach_obs - cum_mass_out_dt[11]) ** 2

                    except TypeError:
                        try:
                            med12_fresh_error = (cum_mass_out_dt[11] - leach_obs) ** 2
                            SS2 = (leach_obs - cum_mass_out_dt[11]) ** 2
                            # med12_fresh_error_prc = ((cum_mass_out_dt[11] - leach_obs) / leach_obs) * 100
                        except TypeError:
                            try:
                                med12_fresh_error = (mass_overflow_dt[11] - pond_obs) ** 2
                                SS2 = 0.
                                # med12_fresh_error_prc = 0
                            except TypeError:
                                med12_fresh_error = 0.
                                SS2 = 0.
                                # med12_fresh_error_prc = 0

                    temp_out_med12_0d = mass_out_dt
                    temp_cum_out_med12_0d = cum_mass_out_dt
                    temp_overflow_med12_0d = mass_overflow_dt

                elif index == 3:
                    try:
                        med12_aged__error = (cum_mass_out_dt[11] - leach_obs) ** 2 + \
                                            (mass_overflow_dt[11] - pond_obs) ** 2
                        SS3 = (leach_obs - cum_mass_out_dt[11]) ** 2

                    except TypeError:
                        try:
                            med12_aged__error = (cum_mass_out_dt[11] - leach_obs) ** 2
                            SS3 = (leach_obs - cum_mass_out_dt[11]) ** 2
                            # med12_aged__error_prc = ((cum_mass_out_dt[11] - leach_obs) / leach_obs) * 100
                            # print("Index:", index)
                            # print("Pest dict: ", pest_dict['d_med12_1d'][1])
                            # print("SS3: ", cum_mass_out_dt[11], leach_obs)
                        except TypeError:
                            try:
                                med12_aged__error = (mass_overflow_dt[11] - pond_obs) ** 2
                                SS3 = 0.
                                # med12_aged__error_prc = 0.
                            except TypeError:
                                med12_aged__error = 0
                                SS3 = 0.
                                # med12_aged__error_prc = 0.

                    temp_out_med12_1d = mass_out_dt
                    temp_cum_out_med12_1d = cum_mass_out_dt
                    temp_overflow_med12_1d = mass_overflow_dt

                elif index == 4:
                    try:
                        med30_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                      (mass_overflow_dt[-1] - pond_obs) ** 2
                        SS4 = (leach_obs - cum_mass_out_dt[-1])**2
                    except TypeError:
                        try:
                            med30_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                            SS4 = (leach_obs - cum_mass_out_dt[-1])**2
                            # med30_fresh__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                        except TypeError:
                            try:
                                med30_fresh__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                                SS4 = 0.
                                # med30_fresh__error_prc = 0.
                            except TypeError:
                                med30_fresh__error = 0
                                SS4 = 0.
                                # med30_fresh__error_prc = 0.

                    temp_out_med30_0d = mass_out_dt
                    temp_cum_out_med30_0d = cum_mass_out_dt
                    temp_overflow_med30_0d = mass_overflow_dt

                elif index == 5:
                    try:
                        med30_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                      (mass_overflow_dt[-1] - pond_obs) ** 2
                        SS5 = (leach_obs - cum_mass_out_dt[-1])**2

                    except TypeError:
                        try:
                            med30_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                            SS5 = (cum_mass_out_dt[-1] - leach_obs)**2
                            # med30_aged__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                        except TypeError:
                            try:
                                med30_aged__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                                SS5 = 0.
                                # med30_aged__error_prc = 0.
                            except TypeError:
                                med30_aged__error = 0.
                                SS5 = 0.
                                # med30_aged__error_prc = 0.

                    temp_out_med30_1d = mass_out_dt
                    temp_cum_out_med30_1d = cum_mass_out_dt
                    temp_overflow_med30_1d = mass_overflow_dt

                elif index == 6:
                    try:
                        low_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                    (mass_overflow_dt[-1] - pond_obs) ** 2
                        SS6 = (cum_mass_out_dt[-1] - leach_obs) ** 2

                    except TypeError:
                        try:
                            low_fresh__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                            SS6 = (cum_mass_out_dt[-1] - leach_obs) ** 2
                            # low_fresh__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                        except TypeError:
                            try:
                                low_fresh__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                                SS6 = 0.
                                # low_fresh__error_prc = 0.
                            except TypeError:
                                low_fresh__error = 0.
                                SS6 = 0.
                                # low_fresh__error_prc = 0.

                    temp_out_low_0d = mass_out_dt
                    temp_cum_out_low_0d = cum_mass_out_dt
                    temp_overflow_low_0d = mass_overflow_dt

                elif index == 7:
                    try:
                        low_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2 + \
                                    (mass_overflow_dt[-1] - pond_obs) ** 2
                        SS7 = (cum_mass_out_dt[-1] - leach_obs) ** 2
                    except TypeError:
                        try:
                            low_aged__error = (cum_mass_out_dt[-1] - leach_obs) ** 2
                            SS7 = (cum_mass_out_dt[-1] - leach_obs) ** 2
                            # low_aged__error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                        except TypeError:
                            try:
                                low_aged__error = (mass_overflow_dt[-1] - pond_obs)
                                SS7 = 0.
                                # low_aged__error_prc = 0.
                            except TypeError:
                                low_aged__error = 0.
                                SS7 = 0.
                                # low_aged__error_prc = 0.

                    temp_out_low_1d = mass_out_dt
                    temp_cum_out_low_1d = cum_mass_out_dt
                    temp_overflow_low_1d = mass_overflow_dt

                else:
                    print("Index number error")

            SS_test = SS0 + SS1 + SS2 + SS3 + SS4 + SS5 + SS6 + SS7

            # if error_test < error:
            if SS_test < error:
                # error = error_test
                error = SS_test
                SSres = SS0 + SS1 + SS2 + SS3 + SS4 + SS5 + SS6 + SS7
                r_squared = 1 - (SSres/SStot)
                kd_chosen = Kd[k]
                Qmax = Q
                num_kd = k
                # print("error = ", error,
                #      "log Kd = ", log10(kd_chosen), "Num: ", num_kd+1)
                # Store final variables

                # a) Cumulative mass leached
                high_fresh_cum_mass_out_dt = temp_cum_out_high_0d
                try:
                    high_fresh_error_prc = ((high_fresh_cum_mass_out_dt[5] - pest_dict['a_high_0d'][1]) / pest_dict['a_high_0d'][1]) * 100
                except TypeError:
                    high_fresh_error_prc = "No obs | < LD "

                high_aged_cum_mass_out_dt = temp_cum_out_high_1d
                try:
                    high_aged_error_prc = ((high_aged_cum_mass_out_dt[5] - pest_dict['b_high_1d'][1]) / pest_dict['b_high_1d'][1]) * 100
                except TypeError:
                    high_aged_error_prc = "No obs | < LD "

                med12_fresh_cum_mass_out_dt = temp_cum_out_med12_0d
                try:
                    med12_fresh_error_prc = ((med12_fresh_cum_mass_out_dt[11] - pest_dict['c_med12_0d'][1]) / pest_dict['c_med12_0d'][1]) * 100
                except TypeError:
                    med12_fresh_error_prc = "No obs | < LD "

                med12_aged_cum_mass_out_dt = temp_cum_out_med12_1d
                try:
                    med12_aged__error_prc = ((med12_aged_cum_mass_out_dt[11] - pest_dict['d_med12_1d'][1]) / pest_dict['d_med12_1d'][1]) * 100
                except TypeError:
                    med12_aged__error_prc = "No obs | < LD "

                med30_fresh_cum_mass_out_dt = temp_cum_out_med30_0d
                try:
                    med30_fresh__error_prc = ((med30_fresh_cum_mass_out_dt[-1] - pest_dict['e_med30_0d'][1]) / pest_dict['e_med30_0d'][1]) * 100
                except TypeError:
                    med30_fresh__error_prc = "No obs | < LD "

                med30_aged_cum_mass_out_dt = temp_cum_out_med30_1d
                try:
                    med30_aged__error_prc = ((med30_aged_cum_mass_out_dt[-1] - pest_dict['f_med30_1d'][1]) / pest_dict['f_med30_1d'][1]) * 100
                except TypeError:
                    med30_aged__error_prc = "No obs | < LD "

                low_fresh_cum_mas_out_dt = temp_cum_out_low_0d
                try:
                    low_fresh__error_prc = ((low_fresh_cum_mas_out_dt[-1] - pest_dict['g_low_0d'][1]) / pest_dict['g_low_0d'][1]) * 100
                except TypeError:
                    low_fresh__error_prc = "No obs | < LD "

                low_aged_cum_mas_out_dt = temp_cum_out_low_1d
                try:
                    low_aged__error_prc = ((low_aged_cum_mas_out_dt[-1] - pest_dict['h_low_1d'][1]) / pest_dict['h_low_1d'][1]) * 100
                except TypeError:
                    low_aged__error_prc = "No obs | < LD "

                # b) Mass in ponded water
                high_fresh_overmass_dt = temp_overflow_high_0d
                high_aged_overmass_dt = temp_overflow_high_1d
                med12_fresh_overmass_dt = temp_overflow_med12_0d
                med12_aged_overmass_dt = temp_overflow_med12_1d
                med30_fresh_overmass_dt = temp_overflow_med30_0d
                med30_aged_overmass_dt = temp_overflow_med30_1d
                low_fresh_overmass_dt = temp_overflow_low_0d
                low_aged_overmass_dt = temp_overflow_low_1d

            else:
                continue

    if isFirstCycle:
        print("1st Pulse")
    else:
        print("2nd Pulse")

    if isLiving:
        print("Living")
        fom_crop_untreat = 5.51 / 100.0
        foc = 0.58 * fom_crop_untreat
    else:
        print("Sterile")
        fom_crop_sterile = 3.87 / 100.0
        foc = 0.58 * fom_crop_sterile

    print("--------------------------------------------")
    print("Koc tested: ", np.array(Kd)/10**3/foc)  # Convert back to cm3/g
    print("Best log Kd: ", log10(kd_chosen/10**3), "cm3/g", "( Num: ", num_kd + 1, ")", "\n",
          "R2: ", r_squared)
    print("Best Cmax: ", Qmax)
    print("Scenario - modality - Predicted error prcnt (%) | Predicted | Observed |")
    print("--------------------------------------------")
    print("(A) 135 mm/h - Fresh ", high_fresh_error_prc, "|", high_fresh_cum_mass_out_dt[5], "|", pest_dict['a_high_0d'][1])
    print("(A) 135 mm/h - Aged ", high_aged_error_prc, "|", high_aged_cum_mass_out_dt[5], "|", pest_dict['b_high_1d'][1])
    print("(B) 55 mm/h - Fresh ", med12_fresh_error_prc, "|", med12_fresh_cum_mass_out_dt[11], "|", pest_dict['c_med12_0d'][1])
    print("(B) 55 mm/h - Aged ", med12_aged__error_prc, "|", med12_aged_cum_mass_out_dt[11], "|", pest_dict['d_med12_1d'][1])
    print("(C) 55 mm/h - Fresh ", med30_fresh__error_prc, "|", med30_fresh_cum_mass_out_dt[-1], "|", pest_dict['e_med30_0d'][1])
    print("(C) 55 mm/h - Aged ", med30_aged__error_prc, "|", med12_aged_cum_mass_out_dt[-1], "|", pest_dict['f_med30_1d'][1])
    print("(D) 30 mm/h - Fresh ", low_fresh__error_prc, "|", high_fresh_cum_mass_out_dt[-1], "|", pest_dict['g_low_0d'][1])
    print("(D) 30 mm/h - Aged ", low_aged__error_prc, "|", high_aged_cum_mass_out_dt[-1], "|", pest_dict['h_low_1d'][1])

    # print ("SS1 to SS7: ", SS1, SS2, SS3, SS4, SS5, SS6, SS7)
    # print ("SSres: ", SSres)
    # print ("SStot: ", SStot)
    # print("mean mass: ", mean_mass, "sum_mass:", sum_mass)

    return stackdata16(cum_time_30min,
                       high_fresh_cum_mass_out_dt, high_aged_cum_mass_out_dt, # leached
                       med12_fresh_cum_mass_out_dt, med12_aged_cum_mass_out_dt,
                       med30_fresh_cum_mass_out_dt, med30_aged_cum_mass_out_dt,
                       low_fresh_cum_mas_out_dt, low_aged_cum_mas_out_dt,
                       high_fresh_overmass_dt, high_aged_overmass_dt,  # ponded
                       med12_fresh_overmass_dt, med12_aged_overmass_dt,
                       med30_fresh_overmass_dt, med30_aged_overmass_dt,
                       low_fresh_overmass_dt, low_aged_overmass_dt)



def freundlich(Kd, a, # Kd list range to test
               pb,
               ov_sat,
               water_data,
               area, soil_height,
               mass_ini,
               pest_sol):

    """
    Water data: leaching, 1st pulse
    """
    vol_h2o_sat = area * soil_height * ov_sat
    cum_time_30min = water_data[:, 0]
    leach_135mmh = water_data[:, 1]  # leached volume per timestep @ high intensity
    leach_55mmh = water_data[:, 2]
    leach_30mmh = water_data[:, 3]
    leached_vol = [leach_135mmh, leach_55mmh, leach_30mmh]
    soil_m = 54.047  # [g]

    """ Cummualtive Output series """
    highint_cum_mass_st_out_dt = [0.0] * len(cum_time_30min)
    medint_cum_mass_st_out_dt = [0.0] * len(cum_time_30min)
    lowint_cum_mass_st_out_dt = [0.0] * len(cum_time_30min)

    highint_cum_mass_un_out_dt = [0.0] * len(cum_time_30min)
    medint_cum_mass_un_out_dt = [0.0] * len(cum_time_30min)
    lowint_cum_mass_un_out_dt = [0.0] * len(cum_time_30min)

    """ Output series """
    highint_mass_st_out_dt = [0.0] * len(cum_time_30min)
    medint_mass_st_out_dt = [0.0] * len(cum_time_30min)
    lowint_mass_st_out_dt = [0.0] * len(cum_time_30min)

    highint_mass_un_out_dt = [0.0] * len(cum_time_30min)
    medint_mass_un_out_dt = [0.0] * len(cum_time_30min)
    lowint_mass_un_out_dt = [0.0] * len(cum_time_30min)

    for m in range(len(mass_ini)):  # for each initial mass (sterile & untreat)
        error = 10**9
        #  Check best Kd fit from selected range
        for k in range(len(Kd)):
            #  Start a new intensity scenario
            for i in range(3):
                mass_tot = mass_ini[m]
                # Assumed at equilibrium at t=0
                # based on initial Kd (# Whelan et al. 1987 (p4.9))
                conc_soil_old = \
                    ((mass_ini[m] / soil_m)*pb*Kd[k])/(ov_sat + pb*Kd[k])
                conc_liq_old = conc_soil_old / Kd[k]
                cum_mass_out = 0
                mass_out_dt = []
                cum_mass_out_dt = []
                #  Compute leaching for each time step under scenario "i"
                for t in range(len(cum_time_30min)):
                    # McGrath leaching, plus Freundlich:
                    r_factor = 1 + (pb * a * Kd[k] * conc_liq_old ** (a - 1)) / ov_sat
                    conc_liq_new = conc_liq_old * exp(-(leached_vol[i][t]) / (r_factor * vol_h2o_sat))
                    mass_out = (conc_liq_old - conc_liq_new) * vol_h2o_sat
                    cum_mass_out += mass_out
                    mass_out_dt.append(mass_out)
                    cum_mass_out_dt.append(cum_mass_out)
                    # re-equilibrate
                    mass_tot_new = mass_tot - mass_out
                    conc_soil_old = Kd[k]*conc_liq_new**a  ## Attention here.. not being updated.
                    conc_liq_old = conc_liq_new
                    mass_tot = mass_tot_new

                # Calculate error of intensity i, @ (kd, m)
                if i == 0 and m == 0:
                    highint_error6_st = (cum_mass_out_dt[5] - pest_sol[m][i]) ** 2
                    temp_cum_out_high_st = cum_mass_out_dt
                    temp_out_high_st = mass_out_dt
                    # print("Mass balance for intensity ", i + 1, "and Kd: ", Kd[k], "is ",
                    #      abs(mass_ini[m] - mass_tot - cum_mass_out_dt[-1]) < 1*10**-6)
                elif i == 1 and m == 0:
                    medint_error12_st = (cum_mass_out_dt[11] - pest_sol[m][i]) ** 2
                    medint_error30_st = (cum_mass_out_dt[-1] - pest_sol[m][i + 1]) ** 2
                    temp_cum_out_med_st = cum_mass_out_dt
                    temp_out_med_st = mass_out_dt
                    # print("Mass balance for intensity ", i + 1, "and Kd: ", Kd[k], "is ",
                    #      abs(mass_ini[m] - mass_tot - cum_mass_out_dt[-1]) < 1 * 10 ** -6)
                elif i == 2 and m == 0:
                    lowint_error30_st = (cum_mass_out_dt[-1] - pest_sol[m][i + 1]) ** 2
                    temp_cum_out_low_st = cum_mass_out_dt
                    temp_out_low_st = mass_out_dt
                elif i == 0 and m == 1:
                    highint_error6_un = (cum_mass_out_dt[5] - pest_sol[m][i]) ** 2
                    temp_cum_out_high_un = cum_mass_out_dt
                    temp_out_high_un = mass_out_dt
                elif i == 1 and m == 1:
                    medint_error12_un = (cum_mass_out_dt[11] - pest_sol[m][i]) ** 2
                    medint_error30_un = (cum_mass_out_dt[-1] - pest_sol[m][i + 1]) ** 2
                    temp_cum_out_med_un = cum_mass_out_dt
                    temp_out_med_un = mass_out_dt
                elif i == 2 and m == 1:
                    lowint_error30_un = (cum_mass_out_dt[-1] - pest_sol[m][i + 1]) ** 2
                    temp_cum_out_low_un = cum_mass_out_dt
                    temp_out_low_un = mass_out_dt

            if m == 0:
                error_test = (((highint_error6_st +
                                medint_error12_st + medint_error30_st +
                                lowint_error30_st) / 4) ** 0.5) / (max(pest_sol[m]) - min(pest_sol[m]))

                if error_test < error:
                    error = error_test
                    error_st = error
                    kd_sterile = Kd[k]
                    num_st = k
                    # print("error = ", error,
                    #      "log Kd (sterile) = ", log10(kd_sterile), "Num: ", num_st+1)
                    highint_cum_mass_st_out_dt = temp_cum_out_high_st
                    highint_mass_st_out_dt = temp_out_high_st
                    medint_cum_mass_st_out_dt = temp_cum_out_med_st
                    medint_mass_st_out_dt = temp_out_med_st
                    lowint_cum_mass_st_out_dt = temp_cum_out_low_st
                    lowint_mass_st_out_dt = temp_out_low_st
                else:
                    continue

            elif m == 1:
                error_test = (((highint_error6_un +
                                medint_error12_un + medint_error30_un +
                                lowint_error30_un) / 4) ** 0.5) / (max(pest_sol[m]) - min(pest_sol[m]))
                if error_test < error:
                    error = error_test
                    error_un = error
                    kd_untreat = Kd[k]
                    num_un = k
                    # print("error = ", error,
                    #      "log Kd (untreat) = ", log10(kd_untreat), "Num: ", num_un+1)
                    highint_cum_mass_un_out_dt = temp_cum_out_high_un
                    highint_mass_un_out_dt = temp_out_high_un
                    medint_cum_mass_un_out_dt = temp_cum_out_med_un
                    medint_mass_un_out_dt = temp_out_med_un
                    lowint_cum_mass_un_out_dt = temp_cum_out_low_un
                    lowint_mass_un_out_dt = temp_out_low_un
                else:
                    continue
            print("Best log Kd (sterile): ", log10(kd_sterile), "( Num: ", num_st + 1, ")", "\n",
                  "Error: ", error_st)
            print("Best log Kd (untreat): ", log10(kd_untreat), "( Num: ", num_un + 1, ")", "\n",
                  "Error: ", error_un)
            return [stackdata6(cum_time_30min,
                               highint_cum_mass_st_out_dt, medint_cum_mass_st_out_dt, lowint_cum_mass_st_out_dt,
                               highint_cum_mass_un_out_dt, medint_cum_mass_un_out_dt, lowint_cum_mass_un_out_dt),
                    stackdata6(cum_time_30min,
                               highint_mass_st_out_dt, medint_mass_st_out_dt, lowint_mass_st_out_dt,
                               highint_mass_un_out_dt, medint_mass_un_out_dt, lowint_mass_un_out_dt)]


