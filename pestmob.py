from math import exp, log10
from hydroplots import *
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
    leach_6min = water_data[:, 1]  # leached volume per timestep @ high intensity
    leach_12min = water_data[:, 2]
    leach_30min = water_data[:, 3]
    leached_vol = [leach_6min, leach_12min, leach_30min]
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
    leach_6min = water_data[:, 1]  # leached volume per timestep @ high intensity
    leach_12min = water_data[:, 2]
    leach_30min = water_data[:, 3]
    leached_vol = [leach_6min, leach_12min, leach_30min]
    soil_m = 54.047  # [g]

    """ Cummualtive Output series """
    highint_cum_mass_st_out_dt = [0.0] * 119
    medint_cum_mass_st_out_dt = [0.0] * 119
    lowint_cum_mass_st_out_dt = [0.0] * 119

    highint_cum_mass_un_out_dt = [0.0] * 119
    medint_cum_mass_un_out_dt = [0.0] * 119
    lowint_cum_mass_un_out_dt = [0.0] * 119

    """ Output series """
    highint_mass_st_out_dt = [0.0] * 119
    medint_mass_st_out_dt = [0.0] * 119
    lowint_mass_st_out_dt = [0.0] * 119

    highint_mass_un_out_dt = [0.0] * 119
    medint_mass_un_out_dt = [0.0] * 119
    lowint_mass_un_out_dt = [0.0] * 119

    for m in range(len(mass_ini)):  # for each initial mass (sterile & untreat)
        error = 10**9
        #  Check best Kd fit from selected range
        for k in range(len(Kd)):
            r_factor = 1 + (pb * Kd[k]) / ov_sat
            #  Start a new intensity scenario
            for i in range(3):
                mass_tot = mass_ini[m]
                conc_soil_old = \
                    ((mass_ini[m] / soil_m)*pb*Kd[k])/(ov_sat + pb*Kd[k])  # Whelan et al. 1987 (p4.9)
                conc_liq_old = conc_soil_old / Kd[k]
                cum_mass_out = 0
                mass_out_dt = []
                cum_mass_out_dt = []
                #  Compute leaching for each time step under scenario "i"
                for t in range(len(cum_time_30min)):
                    # McGrath leaching:
                    conc_liq_new = conc_liq_old * exp(-(leached_vol[i][t]) / (r_factor * vol_h2o_sat))
                    mass_out = (conc_liq_old - conc_liq_new) * vol_h2o_sat
                    cum_mass_out += mass_out
                    mass_out_dt.append(mass_out)
                    cum_mass_out_dt.append(cum_mass_out)
                    # re-equilibrate
                    mass_tot_new = mass_tot - mass_out
                    conc_soil_old = ((mass_tot_new / soil_m)*pb*Kd[k])/(ov_sat + pb*Kd[k])
                    conc_liq_old = conc_soil_old / Kd[k]
                    mass_tot = mass_tot_new

                # Calculate error of intensity i, @ (kd, m)
                if i == 0 and m == 0:
                    highint_error6_st = (cum_mass_out_dt[23] - pest_sol[m][i]) ** 2
                    temp_cum_out_high_st = cum_mass_out_dt
                    temp_out_high_st = mass_out_dt
                    # print("Mass balance for intensity ", i + 1, "and Kd: ", Kd[k], "is ",
                    #      abs(mass_ini[m] - mass_tot - cum_mass_out_dt[-1]) < 1*10**-6)
                elif i == 1 and m == 0:
                    medint_error12_st = (cum_mass_out_dt[47] - pest_sol[m][i]) ** 2
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
                    highint_error6_un = (cum_mass_out_dt[23] - pest_sol[m][i]) ** 2
                    temp_cum_out_high_un = cum_mass_out_dt
                    temp_out_high_un = mass_out_dt
                elif i == 1 and m == 1:
                    medint_error12_un = (cum_mass_out_dt[47] - pest_sol[m][i]) ** 2
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
                              lowint_error30_st)/4) ** 0.5)/(max(pest_sol[m])-min(pest_sol[m]))

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
                              lowint_error30_un)/4) ** 0.5)/(max(pest_sol[m])-min(pest_sol[m]))
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
    print("Best log Kd (sterile): ", log10(kd_sterile), "( Num: ", num_st+1, ")", "\n",
          "Error: ", error_st)
    print("Best log Kd (untreat): ", log10(kd_untreat), "( Num: ", num_un+1, ")", "\n",
          "Error: ", error_un)
    return [stackdata6(cum_time_30min,
                       highint_cum_mass_st_out_dt, medint_cum_mass_st_out_dt, lowint_cum_mass_st_out_dt,
                       highint_cum_mass_un_out_dt, medint_cum_mass_un_out_dt, lowint_cum_mass_un_out_dt),
            stackdata6(cum_time_30min,
                       highint_mass_st_out_dt, medint_mass_st_out_dt, lowint_mass_st_out_dt,
                       highint_mass_un_out_dt, medint_mass_un_out_dt, lowint_mass_un_out_dt)]


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
    leach_6min = water_data[:, 1]  # leached volume per timestep @ high intensity
    leach_12min = water_data[:, 2]
    leach_30min = water_data[:, 3]
    leached_vol = [leach_6min, leach_12min, leach_30min]
    soil_m = 54.047  # [g]

    """ Cummualtive Output series """
    highint_cum_mass_st_out_dt = [0.0] * 119
    medint_cum_mass_st_out_dt = [0.0] * 119
    lowint_cum_mass_st_out_dt = [0.0] * 119

    highint_cum_mass_un_out_dt = [0.0] * 119
    medint_cum_mass_un_out_dt = [0.0] * 119
    lowint_cum_mass_un_out_dt = [0.0] * 119

    """ Output series """
    highint_mass_st_out_dt = [0.0] * 119
    medint_mass_st_out_dt = [0.0] * 119
    lowint_mass_st_out_dt = [0.0] * 119

    highint_mass_un_out_dt = [0.0] * 119
    medint_mass_un_out_dt = [0.0] * 119
    lowint_mass_un_out_dt = [0.0] * 119

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
                    highint_error6_st = (cum_mass_out_dt[23] - pest_sol[m][i]) ** 2
                    temp_cum_out_high_st = cum_mass_out_dt
                    temp_out_high_st = mass_out_dt
                    # print("Mass balance for intensity ", i + 1, "and Kd: ", Kd[k], "is ",
                    #      abs(mass_ini[m] - mass_tot - cum_mass_out_dt[-1]) < 1*10**-6)
                elif i == 1 and m == 0:
                    medint_error12_st = (cum_mass_out_dt[47] - pest_sol[m][i]) ** 2
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
                    highint_error6_un = (cum_mass_out_dt[23] - pest_sol[m][i]) ** 2
                    temp_cum_out_high_un = cum_mass_out_dt
                    temp_out_high_un = mass_out_dt
                elif i == 1 and m == 1:
                    medint_error12_un = (cum_mass_out_dt[47] - pest_sol[m][i]) ** 2
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



