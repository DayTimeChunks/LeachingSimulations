from math import exp, log10
from hydroplots import *
from mixinglayer import *
import matplotlib.pyplot as plt


def pest_linear_x(
        Kd,  # Kd list range to test
        x,
        initial_condition,
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

    try:
        mean_mass = sum_mass/masses
    except ZeroDivisionError:
        mean_mass = 0

    SStot_fresh = 0
    SStot_aged = 0
    for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) \
            in enumerate(pest_dict.items()):
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

    # Assign water time series
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

            if initial_condition == 'mean':
                mass_tot = mass_i
            elif initial_condition == 'min':
                mass_tot = mass_i - initial_mass_error
            elif initial_condition == 'max':
                mass_tot = mass_i + initial_mass_error

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
                        # med12_aged_error_prc = ((cum_mass_out_dt[11] - leach_obs) / leach_obs) * 100
                        # print("Index:", index)
                        # print("Pest dict: ", pest_dict['d_med12_1d'][1])
                        # print("SS3: ", cum_mass_out_dt[11], leach_obs)
                    except TypeError:
                        try:
                            med12_aged__error = (mass_overflow_dt[11] - pond_obs) ** 2
                            SS3 = 0.
                            # med12_aged_error_prc = 0.
                        except TypeError:
                            med12_aged__error = 0
                            SS3 = 0.
                            # med12_aged_error_prc = 0.

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
                        # med30_aged_error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            med30_aged__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                            SS5 = 0.
                            # med30_aged_error_prc = 0.
                        except TypeError:
                            med30_aged__error = 0.
                            SS5 = 0.
                            # med30_aged_error_prc = 0.

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
                        # low_fresh_error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            low_fresh__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                            SS6 = 0.
                            # low_fresh_error_prc = 0.
                        except TypeError:
                            low_fresh__error = 0.
                            SS6 = 0.
                            # low_fresh_error_prc = 0.

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
                        # low_aged_error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                    except TypeError:
                        try:
                            low_aged__error = (mass_overflow_dt[-1] - pond_obs)
                            SS7 = 0.
                            # low_aged_error_prc = 0.
                        except TypeError:
                            low_aged__error = 0.
                            SS7 = 0.
                            # low_aged_error_prc = 0.

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
            try:
                r_squared_fresh = 1 - (SSres_fresh/SStot_fresh)
            except ZeroDivisionError:
                r_squared_fresh = 1
            kd_chosen_fresh = Kd[k]
            x_factor = x
            num_kd_fresh = k

            # a) Cumulative mass leached
            high_fresh_cum_mass_out_dt = temp_cum_out_high_0d
            try:
                high_fresh_error_prc = ((high_fresh_cum_mass_out_dt[5] - pest_dict['a_high_0d'][1]) / pest_dict['a_high_0d'][1]) * 100
            except TypeError:
                high_fresh_error_prc = "No obs"

            med12_fresh_cum_mass_out_dt = temp_cum_out_med12_0d
            try:
                med12_fresh_error_prc = ((med12_fresh_cum_mass_out_dt[11] - pest_dict['c_med12_0d'][1]) / pest_dict['c_med12_0d'][1]) * 100
            except TypeError:
                med12_fresh_error_prc = "No obs"

            med30_fresh_cum_mass_out_dt = temp_cum_out_med30_0d
            try:
                med30_fresh_error_prc = ((med30_fresh_cum_mass_out_dt[-1] - pest_dict['e_med30_0d'][1]) / pest_dict['e_med30_0d'][1]) * 100
            except TypeError:
                med30_fresh_error_prc = "No obs"

            low_fresh_cum_mas_out_dt = temp_cum_out_low_0d
            try:
                low_fresh_error_prc = ((low_fresh_cum_mas_out_dt[-1] - pest_dict['g_low_0d'][1]) / pest_dict['g_low_0d'][1]) * 100
            except TypeError:
                low_fresh_error_prc = "No obs"

            # b) Mass in ponded water
            high_fresh_overmass_dt = temp_overflow_high_0d
            med12_fresh_overmass_dt = temp_overflow_med12_0d
            med30_fresh_overmass_dt = temp_overflow_med30_0d
            low_fresh_overmass_dt = temp_overflow_low_0d

        if SS_aged < error_aged:
            error_aged = SS_aged
            SSres_aged = SS1 + SS3 + SS5 + SS7
            try:
                r_squared_aged = 1 - (SSres_aged/SStot_aged)
            except ZeroDivisionError:
                r_squared_aged = 1
            kd_chosen_aged = Kd[k]
            num_kd_aged = k

            # a) Cumulative mass leached
            high_aged_cum_mass_out_dt = temp_cum_out_high_1d
            try:
                high_aged_error_prc = ((high_aged_cum_mass_out_dt[5] - pest_dict['b_high_1d'][1]) /
                                       pest_dict['b_high_1d'][1]) * 100
            except TypeError:
                high_aged_error_prc = "No obs"

            med12_aged_cum_mass_out_dt = temp_cum_out_med12_1d
            try:
                med12_aged_error_prc = ((med12_aged_cum_mass_out_dt[11] - pest_dict['d_med12_1d'][1]) /
                                         pest_dict['d_med12_1d'][1]) * 100
            except TypeError:
                med12_aged_error_prc = "No obs"

            med30_aged_cum_mass_out_dt = temp_cum_out_med30_1d
            try:
                med30_aged_error_prc = ((med30_aged_cum_mass_out_dt[-1] - pest_dict['f_med30_1d'][1]) /
                                         pest_dict['f_med30_1d'][1]) * 100
            except TypeError:
                med30_aged_error_prc = "No obs"

            low_aged_cum_mas_out_dt = temp_cum_out_low_1d
            try:
                low_aged_error_prc = (
                                      (low_aged_cum_mas_out_dt[-1] - pest_dict['h_low_1d'][1]) / pest_dict['h_low_1d'][1]) * 100
            except TypeError:
                low_aged_error_prc = "No obs"

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

    if isLiving and soil == 'Alteck':
        print("Living Crop")
        fom_untreat = 5.51 / 100.0
        foc = 0.58 * fom_untreat
    elif not isLiving and soil == 'Alteck':
        print("Sterile Crop")
        fom_sterile = 3.87 / 100.0
        foc = 0.58 * fom_sterile

    elif isLiving and soil == 'Rouff':
        print("Living Vine")
        fom_vine_untreat = 2.93 / 100.0
        foc = 0.58 * fom_vine_untreat
    elif not isLiving and soil == 'Rouff':
        print("Sterile Vine")
        fom_vine_sterile = 3.53 / 100.0
        foc = 0.58 * fom_vine_sterile

    koc_fresh_chosen = (kd_chosen_fresh/x/10**3/foc)
    koc_aged_chosen = kd_chosen_aged/10**3/foc

    print("--------------------------------------------")
    print("Koc tested: ", np.array(Kd)/10**3/foc)  # Convert back to cm3/g
    print("Best Kd (Fresh): ", kd_chosen_fresh/10**3, "cm3/g", "( Num: ", num_kd_fresh + 1, ")", "\n",
          "x factor: ", x, "\n",
          "R2: ", r_squared_fresh, "\n",
          "Effective Koc (fresh) - [cm3/g]:", koc_fresh_chosen)
    print("--------------------------------------------")
    print("Best Kd (Aged): ", kd_chosen_aged / 10 ** 3, "cm3/g", "( Num: ", num_kd_aged + 1, ")", "\n",
          "x factor: No factor considered. ", "\n",
          "R2: ", r_squared_aged, "\n",
          "Effective Koc (aged) [cm3/g]:", koc_aged_chosen)
    print("--------------------------------------------")
    print("--------------------------------------------")
    print("Scenario - modality - Predicted error prcnt (%) | Predicted | Observed |")
    print("--------------------------------------------")
    print("(A) 135 mm/h - Fresh ", high_fresh_error_prc, "|", high_fresh_cum_mass_out_dt[5], "|", pest_dict['a_high_0d'][1])
    print("(A) 135 mm/h - Aged ", high_aged_error_prc, "|", high_aged_cum_mass_out_dt[5], "|", pest_dict['b_high_1d'][1])
    print("(B) 55 mm/h - Fresh ", med12_fresh_error_prc, "|", med12_fresh_cum_mass_out_dt[11], "|", pest_dict['c_med12_0d'][1])
    print("(B) 55 mm/h - Aged ", med12_aged_error_prc, "|", med12_aged_cum_mass_out_dt[11], "|", pest_dict['d_med12_1d'][1])
    print("(C) 55 mm/h - Fresh ", med30_fresh_error_prc, "|", med30_fresh_cum_mass_out_dt[-1], "|", pest_dict['e_med30_0d'][1])
    print("(C) 55 mm/h - Aged ", med30_aged_error_prc, "|", med12_aged_cum_mass_out_dt[-1], "|", pest_dict['f_med30_1d'][1])
    print("(D) 30 mm/h - Fresh ", low_fresh_error_prc, "|", high_fresh_cum_mass_out_dt[-1], "|", pest_dict['g_low_0d'][1])
    print("(D) 30 mm/h - Aged ", low_aged_error_prc, "|", high_aged_cum_mass_out_dt[-1], "|", pest_dict['h_low_1d'][1])

    # print ("SS1 to SS7: ", SS1, SS2, SS3, SS4, SS5, SS6, SS7)
    # print ("SSres: ", SSres)
    # print ("SStot: ", SStot)
    # print("mean mass: ", mean_mass, "sum_mass:", sum_mass)
    stack = stackdata16(cum_time_30min,
                       high_fresh_cum_mass_out_dt, high_aged_cum_mass_out_dt, # leached
                       med12_fresh_cum_mass_out_dt, med12_aged_cum_mass_out_dt,
                       med30_fresh_cum_mass_out_dt, med30_aged_cum_mass_out_dt,
                       low_fresh_cum_mas_out_dt, low_aged_cum_mas_out_dt,
                       high_fresh_overmass_dt, high_aged_overmass_dt,  # ponded
                       med12_fresh_overmass_dt, med12_aged_overmass_dt,
                       med30_fresh_overmass_dt, med30_aged_overmass_dt,
                       low_fresh_overmass_dt, low_aged_overmass_dt)

    return {'Data_stack': stack,
            'R': {'R_fresh': r_squared_fresh, 'R_aged': r_squared_aged},
            'Error_percent': {'high_fresh': high_fresh_error_prc,
                              'med12_fresh': med12_fresh_error_prc,
                              'med30_fresh': med30_fresh_error_prc,
                              'low_fresh': low_fresh_error_prc,
                              'high_aged': high_aged_error_prc,
                              'med12_aged': med12_aged_error_prc,
                              'med30_aged': med30_aged_error_prc,
                              'low_aged': low_aged_error_prc},
            'Koc_fresh': koc_fresh_chosen,
            'Kd_fresh': kd_chosen_fresh/10**3,
            'Koc_aged': koc_aged_chosen,
            'Kd_aged': kd_chosen_aged/10**3}


def pest_Kd_choice(
        Kd,  # Kd coices (length = 2), [fresh, aged]
        x,
        initial_condition,
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

    try:
        mean_mass = sum_mass/masses
    except ZeroDivisionError:
        mean_mass = 0

    SStot_fresh = 0
    SStot_aged = 0
    for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) \
            in enumerate(pest_dict.items()):
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

    # Assign water time series
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

    # define a max error
    error_fresh = 10 ** 9
    error_aged = 10 ** 9
    # Check best Kd fit from selected range

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

        if initial_condition == 'mean':
            mass_tot = mass_i
        elif initial_condition == 'min':
            mass_tot = mass_i - initial_mass_error
        elif initial_condition == 'max':
            mass_tot = mass_i + initial_mass_error

        if index % 2 == 0:
            k = Kd[0]
            # Whelan et al. 1987 (p4.9)
            conc_soil_old = \
                ((mass_tot / soil_m) * pb_i * k/x / (ov + pb_i * (k/x)))  # ug/g
            conc_liq_old = \
                ((mass_tot / soil_m) * pb_i) / (ov + pb_i * k/x)  # ug/g

            r_factor = 1 + (pb_i * k/x) / ovSat
            conc_liq_old = conc_soil_old / (k/x)  # ug/mm3

            #  Compute leaching for each time step under scenario of "index"
            hasleached = 0
            for t in range(len(cum_time_30min)):
                # Update r_factor, based on increasing pb.
                if leached_vol[index][t] > 0:
                    hasleached += 1
                if hasleached > 12:
                    pb_i = pb_f
                    r_factor = 1 + (pb_i * k/x) / ovSat

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

        else:
            k = Kd[1]
            # Whelan et al. 1987 (p4.9)
            conc_soil_old = \
                ((mass_tot / soil_m) * pb_i * k / (ov + pb_i * (k)))  # ug/g
            conc_liq_old = \
                ((mass_tot / soil_m) * pb_i) / (ov + pb_i * k)  # ug/g

            r_factor = 1 + (pb_i * k) / ovSat
            conc_liq_old = conc_soil_old / (k)  # ug/mm3

            #  Compute leaching for each time step under scenario of "index"
            hasleached = 0
            for t in range(len(cum_time_30min)):
                # Update r_factor, based on increasing pb.
                if leached_vol[index][t] > 0:
                    hasleached += 1
                if hasleached > 12:
                    pb_i = pb_f
                    r_factor = 1 + (pb_i * k) / ovSat

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
                    # med12_aged_error_prc = ((cum_mass_out_dt[11] - leach_obs) / leach_obs) * 100
                    # print("Index:", index)
                    # print("Pest dict: ", pest_dict['d_med12_1d'][1])
                    # print("SS3: ", cum_mass_out_dt[11], leach_obs)
                except TypeError:
                    try:
                        med12_aged__error = (mass_overflow_dt[11] - pond_obs) ** 2
                        SS3 = 0.
                        # med12_aged_error_prc = 0.
                    except TypeError:
                        med12_aged__error = 0
                        SS3 = 0.
                        # med12_aged_error_prc = 0.

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
                    # med30_aged_error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                except TypeError:
                    try:
                        med30_aged__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                        SS5 = 0.
                        # med30_aged_error_prc = 0.
                    except TypeError:
                        med30_aged__error = 0.
                        SS5 = 0.
                        # med30_aged_error_prc = 0.

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
                    # low_fresh_error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                except TypeError:
                    try:
                        low_fresh__error = (mass_overflow_dt[-1] - pond_obs) ** 2
                        SS6 = 0.
                        # low_fresh_error_prc = 0.
                    except TypeError:
                        low_fresh__error = 0.
                        SS6 = 0.
                        # low_fresh_error_prc = 0.

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
                    # low_aged_error_prc = ((cum_mass_out_dt[-1] - leach_obs) / leach_obs) * 100
                except TypeError:
                    try:
                        low_aged__error = (mass_overflow_dt[-1] - pond_obs)
                        SS7 = 0.
                        # low_aged_error_prc = 0.
                    except TypeError:
                        low_aged__error = 0.
                        SS7 = 0.
                        # low_aged_error_prc = 0.

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
        try:
            r_squared_fresh = 1 - (SSres_fresh/SStot_fresh)
        except ZeroDivisionError:
            r_squared_fresh = 1
        kd_chosen_fresh = Kd[0]
        x_factor = x
        num_kd_fresh = 0

        # a) Cumulative mass leached
        high_fresh_cum_mass_out_dt = temp_cum_out_high_0d
        try:
            high_fresh_error_prc = ((high_fresh_cum_mass_out_dt[5] - pest_dict['a_high_0d'][1]) / pest_dict['a_high_0d'][1]) * 100
        except TypeError:
            high_fresh_error_prc = "No obs"

        med12_fresh_cum_mass_out_dt = temp_cum_out_med12_0d
        try:
            med12_fresh_error_prc = ((med12_fresh_cum_mass_out_dt[11] - pest_dict['c_med12_0d'][1]) / pest_dict['c_med12_0d'][1]) * 100
        except TypeError:
            med12_fresh_error_prc = "No obs"

        med30_fresh_cum_mass_out_dt = temp_cum_out_med30_0d
        try:
            med30_fresh_error_prc = ((med30_fresh_cum_mass_out_dt[-1] - pest_dict['e_med30_0d'][1]) / pest_dict['e_med30_0d'][1]) * 100
        except TypeError:
            med30_fresh_error_prc = "No obs"

        low_fresh_cum_mas_out_dt = temp_cum_out_low_0d
        try:
            low_fresh_error_prc = ((low_fresh_cum_mas_out_dt[-1] - pest_dict['g_low_0d'][1]) / pest_dict['g_low_0d'][1]) * 100
        except TypeError:
            low_fresh_error_prc = "No obs"

        # b) Mass in ponded water
        high_fresh_overmass_dt = temp_overflow_high_0d
        med12_fresh_overmass_dt = temp_overflow_med12_0d
        med30_fresh_overmass_dt = temp_overflow_med30_0d
        low_fresh_overmass_dt = temp_overflow_low_0d

    if SS_aged < error_aged:
        error_aged = SS_aged
        SSres_aged = SS1 + SS3 + SS5 + SS7
        try:
            r_squared_aged = 1 - (SSres_aged/SStot_aged)
        except ZeroDivisionError:
            r_squared_aged = 1
        kd_chosen_aged = Kd[1]
        num_kd_aged = 1

        # a) Cumulative mass leached
        high_aged_cum_mass_out_dt = temp_cum_out_high_1d
        try:
            high_aged_error_prc = ((high_aged_cum_mass_out_dt[5] - pest_dict['b_high_1d'][1]) /
                                   pest_dict['b_high_1d'][1]) * 100
        except TypeError:
            high_aged_error_prc = "No obs"

        med12_aged_cum_mass_out_dt = temp_cum_out_med12_1d
        try:
            med12_aged_error_prc = ((med12_aged_cum_mass_out_dt[11] - pest_dict['d_med12_1d'][1]) /
                                     pest_dict['d_med12_1d'][1]) * 100
        except TypeError:
            med12_aged_error_prc = "No obs"

        med30_aged_cum_mass_out_dt = temp_cum_out_med30_1d
        try:
            med30_aged_error_prc = ((med30_aged_cum_mass_out_dt[-1] - pest_dict['f_med30_1d'][1]) /
                                     pest_dict['f_med30_1d'][1]) * 100
        except TypeError:
            med30_aged_error_prc = "No obs"

        low_aged_cum_mas_out_dt = temp_cum_out_low_1d
        try:
            low_aged_error_prc = (
                                  (low_aged_cum_mas_out_dt[-1] - pest_dict['h_low_1d'][1]) / pest_dict['h_low_1d'][1]) * 100
        except TypeError:
            low_aged_error_prc = "No obs"

        high_aged_overmass_dt = temp_overflow_high_1d
        med12_aged_overmass_dt = temp_overflow_med12_1d
        med30_aged_overmass_dt = temp_overflow_med30_1d
        low_aged_overmass_dt = temp_overflow_low_1d

    if isFirstCycle:
        print("1st Pulse")
    else:
        print("2nd Pulse")

    if isLiving and soil == 'Alteck':
        print("Living Crop")
        fom_untreat = 5.51 / 100.0
        foc = 0.58 * fom_untreat
    elif not isLiving and soil == 'Alteck':
        print("Sterile Crop")
        fom_sterile = 3.87 / 100.0
        foc = 0.58 * fom_sterile

    elif isLiving and soil == 'Rouff':
        print("Living Vine")
        fom_vine_untreat = 2.93 / 100.0
        foc = 0.58 * fom_vine_untreat
    elif not isLiving and soil == 'Rouff':
        print("Sterile Vine")
        fom_vine_sterile = 3.53 / 100.0
        foc = 0.58 * fom_vine_sterile

    koc_fresh_chosen = (kd_chosen_fresh/x)/10**3/foc
    koc_aged_chosen = kd_chosen_aged/10**3/foc

    print("--------------------------------------------")
    print("Koc tested: ", np.array(Kd)/10**3/foc)  # Convert back to cm3/g
    print("Best Kd (Fresh): ", kd_chosen_fresh/10**3, "cm3/g", "( Num: ", num_kd_fresh + 1, ")", "\n",
          "x factor: ", x, "\n",
          "R2: ", r_squared_fresh, "\n",
          "Effective Koc (fresh) - [cm3/g]:", koc_fresh_chosen)
    print("--------------------------------------------")
    print("Best Kd (Aged): ", kd_chosen_aged / 10 ** 3, "cm3/g", "( Num: ", num_kd_aged + 1, ")", "\n",
          "x factor: No factor considered. ", "\n",
          "R2: ", r_squared_aged, "\n",
          "Effective Koc (aged) [cm3/g]:", koc_aged_chosen)
    print("--------------------------------------------")
    print("--------------------------------------------")
    print("Scenario - modality - Predicted error prcnt (%) | Predicted | Observed |")
    print("--------------------------------------------")
    print("(A) 135 mm/h - Fresh ", high_fresh_error_prc, "|", high_fresh_cum_mass_out_dt[5], "|", pest_dict['a_high_0d'][1])
    print("(A) 135 mm/h - Aged ", high_aged_error_prc, "|", high_aged_cum_mass_out_dt[5], "|", pest_dict['b_high_1d'][1])
    print("(B) 55 mm/h - Fresh ", med12_fresh_error_prc, "|", med12_fresh_cum_mass_out_dt[11], "|", pest_dict['c_med12_0d'][1])
    print("(B) 55 mm/h - Aged ", med12_aged_error_prc, "|", med12_aged_cum_mass_out_dt[11], "|", pest_dict['d_med12_1d'][1])
    print("(C) 55 mm/h - Fresh ", med30_fresh_error_prc, "|", med30_fresh_cum_mass_out_dt[-1], "|", pest_dict['e_med30_0d'][1])
    print("(C) 55 mm/h - Aged ", med30_aged_error_prc, "|", med12_aged_cum_mass_out_dt[-1], "|", pest_dict['f_med30_1d'][1])
    print("(D) 30 mm/h - Fresh ", low_fresh_error_prc, "|", high_fresh_cum_mass_out_dt[-1], "|", pest_dict['g_low_0d'][1])
    print("(D) 30 mm/h - Aged ", low_aged_error_prc, "|", high_aged_cum_mass_out_dt[-1], "|", pest_dict['h_low_1d'][1])

    # print ("SS1 to SS7: ", SS1, SS2, SS3, SS4, SS5, SS6, SS7)
    # print ("SSres: ", SSres)
    # print ("SStot: ", SStot)
    # print("mean mass: ", mean_mass, "sum_mass:", sum_mass)
    stack = stackdata16(cum_time_30min,
                       high_fresh_cum_mass_out_dt, high_aged_cum_mass_out_dt, # leached
                       med12_fresh_cum_mass_out_dt, med12_aged_cum_mass_out_dt,
                       med30_fresh_cum_mass_out_dt, med30_aged_cum_mass_out_dt,
                       low_fresh_cum_mas_out_dt, low_aged_cum_mas_out_dt,
                       high_fresh_overmass_dt, high_aged_overmass_dt,  # ponded
                       med12_fresh_overmass_dt, med12_aged_overmass_dt,
                       med30_fresh_overmass_dt, med30_aged_overmass_dt,
                       low_fresh_overmass_dt, low_aged_overmass_dt)

    return {'Data_stack': stack,
            'R': {'R_fresh': r_squared_fresh, 'R_aged': r_squared_aged},
            'Error_percent': {'high_fresh': high_fresh_error_prc,
                              'med12_fresh': med12_fresh_error_prc,
                              'med30_fresh': med30_fresh_error_prc,
                              'low_fresh': low_fresh_error_prc,
                              'high_aged': high_aged_error_prc,
                              'med12_aged': med12_aged_error_prc,
                              'med30_aged': med30_aged_error_prc,
                              'low_aged': low_aged_error_prc},
            'Koc_fresh': koc_fresh_chosen,
            'Kd_fresh': kd_chosen_fresh/10**3,
            'Koc_aged': koc_aged_chosen,
            'Kd_aged': kd_chosen_aged/10**3}


def extract_goodness_pest(
        pest_data,
        output_pesti_params,
        soil_short,
        isLiving,
        isFirstCycle
):
    if isLiving and isFirstCycle:
        # Living
        output_pesti_params[soil_short]['R1']['LF'] = pest_data['R']['R_fresh']
        output_pesti_params[soil_short]['R1']['LA'] = pest_data['R']['R_aged']

        output_pesti_params[soil_short]['Err1']['LF'].append(pest_data['Error_percent']['high_fresh'])
        output_pesti_params[soil_short]['Err1']['LA'].append(pest_data['Error_percent']['high_aged'])

        output_pesti_params[soil_short]['Err1']['LF'].append(pest_data['Error_percent']['med12_fresh'])
        output_pesti_params[soil_short]['Err1']['LA'].append(pest_data['Error_percent']['med12_aged'])

        output_pesti_params[soil_short]['Err1']['LF'].append(pest_data['Error_percent']['med30_fresh'])
        output_pesti_params[soil_short]['Err1']['LA'].append(pest_data['Error_percent']['med30_aged'])

        output_pesti_params[soil_short]['Err1']['LF'].append(pest_data['Error_percent']['low_fresh'])
        output_pesti_params[soil_short]['Err1']['LA'].append(pest_data['Error_percent']['low_aged'])

        output_pesti_params[soil_short]['Koc1']['LF'] = pest_data['Koc_fresh']
        output_pesti_params[soil_short]['Kd1']['LF'] = pest_data['Kd_fresh']

        output_pesti_params[soil_short]['Koc1']['LA'] = pest_data['Koc_aged']
        output_pesti_params[soil_short]['Kd1']['LA'] = pest_data['Kd_aged']

    elif not isLiving and isFirstCycle:
        # Sterile
        output_pesti_params[soil_short]['R1']['SF'] = pest_data['R']['R_fresh']
        output_pesti_params[soil_short]['R1']['SA'] = pest_data['R']['R_aged']

        output_pesti_params[soil_short]['Err1']['SF'].append(pest_data['Error_percent']['high_fresh'])
        output_pesti_params[soil_short]['Err1']['SA'].append(pest_data['Error_percent']['high_aged'])

        output_pesti_params[soil_short]['Err1']['SF'].append(pest_data['Error_percent']['med12_fresh'])
        output_pesti_params[soil_short]['Err1']['SA'].append(pest_data['Error_percent']['med12_aged'])

        output_pesti_params[soil_short]['Err1']['SF'].append(pest_data['Error_percent']['med30_fresh'])
        output_pesti_params[soil_short]['Err1']['SA'].append(pest_data['Error_percent']['med30_aged'])

        output_pesti_params[soil_short]['Err1']['SF'].append(pest_data['Error_percent']['low_fresh'])
        output_pesti_params[soil_short]['Err1']['SA'].append(pest_data['Error_percent']['low_aged'])

        output_pesti_params[soil_short]['Koc1']['SF'] = pest_data['Koc_fresh']
        output_pesti_params[soil_short]['Kd1']['SF'] = pest_data['Kd_fresh']

        output_pesti_params[soil_short]['Koc1']['SA'] = pest_data['Koc_aged']
        output_pesti_params[soil_short]['Kd1']['SA'] = pest_data['Kd_aged']

    elif isLiving and not isFirstCycle:
        # Living
        output_pesti_params[soil_short]['R2']['LF'] = pest_data['R']['R_fresh']
        output_pesti_params[soil_short]['R2']['LA'] = pest_data['R']['R_aged']

        output_pesti_params[soil_short]['Err2']['LF'].append(pest_data['Error_percent']['high_fresh'])
        output_pesti_params[soil_short]['Err2']['LA'].append(pest_data['Error_percent']['high_aged'])

        output_pesti_params[soil_short]['Err2']['LF'].append(pest_data['Error_percent']['med12_fresh'])
        output_pesti_params[soil_short]['Err2']['LA'].append(pest_data['Error_percent']['med12_aged'])

        output_pesti_params[soil_short]['Err2']['LF'].append(pest_data['Error_percent']['med30_fresh'])
        output_pesti_params[soil_short]['Err2']['LA'].append(pest_data['Error_percent']['med30_aged'])

        output_pesti_params[soil_short]['Err2']['LF'].append(pest_data['Error_percent']['low_fresh'])
        output_pesti_params[soil_short]['Err2']['LA'].append(pest_data['Error_percent']['low_aged'])

        output_pesti_params[soil_short]['Koc2']['LF'] = pest_data['Koc_fresh']
        output_pesti_params[soil_short]['Kd2']['LF'] = pest_data['Kd_fresh']

        output_pesti_params[soil_short]['Koc2']['LA'] = pest_data['Koc_aged']
        output_pesti_params[soil_short]['Kd2']['LA'] = pest_data['Kd_aged']

    elif not isLiving and not isFirstCycle:
        # Sterile
        output_pesti_params[soil_short]['R2']['SF'] = pest_data['R']['R_fresh']
        output_pesti_params[soil_short]['R2']['SA'] = pest_data['R']['R_aged']

        output_pesti_params[soil_short]['Err2']['SF'].append(pest_data['Error_percent']['high_fresh'])
        output_pesti_params[soil_short]['Err2']['SA'].append(pest_data['Error_percent']['high_aged'])

        output_pesti_params[soil_short]['Err2']['SF'].append(pest_data['Error_percent']['med12_fresh'])
        output_pesti_params[soil_short]['Err2']['SA'].append(pest_data['Error_percent']['med12_aged'])

        output_pesti_params[soil_short]['Err2']['SF'].append(pest_data['Error_percent']['med30_fresh'])
        output_pesti_params[soil_short]['Err2']['SA'].append(pest_data['Error_percent']['med30_aged'])

        output_pesti_params[soil_short]['Err2']['SF'].append(pest_data['Error_percent']['low_fresh'])
        output_pesti_params[soil_short]['Err2']['SA'].append(pest_data['Error_percent']['low_aged'])

        output_pesti_params[soil_short]['Koc2']['SF'] = pest_data['Koc_fresh']
        output_pesti_params[soil_short]['Kd2']['SF'] = pest_data['Kd_fresh']

        output_pesti_params[soil_short]['Koc2']['SA'] = pest_data['Koc_aged']
        output_pesti_params[soil_short]['Kd2']['SA'] = pest_data['Kd_aged']

    return output_pesti_params


def extract_goodness_maxmin(
        pest_data,
        initial_condition,
        output_pesti_params,
        soil_short,
        isLiving,
        isFirstCycle
):
    if initial_condition == 'max':
        if isLiving and isFirstCycle:
            # Living
            output_pesti_params[soil_short]['Err1']['LFmax'].append(pest_data['Error_percent']['high_fresh'])
            output_pesti_params[soil_short]['Err1']['LAmax'].append(pest_data['Error_percent']['high_aged'])

            output_pesti_params[soil_short]['Err1']['LFmax'].append(pest_data['Error_percent']['med12_fresh'])
            output_pesti_params[soil_short]['Err1']['LAmax'].append(pest_data['Error_percent']['med12_aged'])

            output_pesti_params[soil_short]['Err1']['LFmax'].append(pest_data['Error_percent']['med30_fresh'])
            output_pesti_params[soil_short]['Err1']['LAmax'].append(pest_data['Error_percent']['med30_aged'])

            output_pesti_params[soil_short]['Err1']['LFmax'].append(pest_data['Error_percent']['low_fresh'])
            output_pesti_params[soil_short]['Err1']['LAmax'].append(pest_data['Error_percent']['low_aged'])

        elif not isLiving and isFirstCycle:
            # Sterile
            output_pesti_params[soil_short]['Err1']['SFmax'].append(pest_data['Error_percent']['high_fresh'])
            output_pesti_params[soil_short]['Err1']['SAmax'].append(pest_data['Error_percent']['high_aged'])

            output_pesti_params[soil_short]['Err1']['SFmax'].append(pest_data['Error_percent']['med12_fresh'])
            output_pesti_params[soil_short]['Err1']['SAmax'].append(pest_data['Error_percent']['med12_aged'])

            output_pesti_params[soil_short]['Err1']['SFmax'].append(pest_data['Error_percent']['med30_fresh'])
            output_pesti_params[soil_short]['Err1']['SAmax'].append(pest_data['Error_percent']['med30_aged'])

            output_pesti_params[soil_short]['Err1']['SFmax'].append(pest_data['Error_percent']['low_fresh'])
            output_pesti_params[soil_short]['Err1']['SAmax'].append(pest_data['Error_percent']['low_aged'])

        elif isLiving and not isFirstCycle:
            # Living
            output_pesti_params[soil_short]['Err2']['LFmax'].append(pest_data['Error_percent']['high_fresh'])
            output_pesti_params[soil_short]['Err2']['LAmax'].append(pest_data['Error_percent']['high_aged'])

            output_pesti_params[soil_short]['Err2']['LFmax'].append(pest_data['Error_percent']['med12_fresh'])
            output_pesti_params[soil_short]['Err2']['LAmax'].append(pest_data['Error_percent']['med12_aged'])

            output_pesti_params[soil_short]['Err2']['LFmax'].append(pest_data['Error_percent']['med30_fresh'])
            output_pesti_params[soil_short]['Err2']['LAmax'].append(pest_data['Error_percent']['med30_aged'])

            output_pesti_params[soil_short]['Err2']['LFmax'].append(pest_data['Error_percent']['low_fresh'])
            output_pesti_params[soil_short]['Err2']['LAmax'].append(pest_data['Error_percent']['low_aged'])

        elif not isLiving and not isFirstCycle:
            # Sterile
            output_pesti_params[soil_short]['Err2']['SFmax'].append(pest_data['Error_percent']['high_fresh'])
            output_pesti_params[soil_short]['Err2']['SAmax'].append(pest_data['Error_percent']['high_aged'])

            output_pesti_params[soil_short]['Err2']['SFmax'].append(pest_data['Error_percent']['med12_fresh'])
            output_pesti_params[soil_short]['Err2']['SAmax'].append(pest_data['Error_percent']['med12_aged'])

            output_pesti_params[soil_short]['Err2']['SFmax'].append(pest_data['Error_percent']['med30_fresh'])
            output_pesti_params[soil_short]['Err2']['SAmax'].append(pest_data['Error_percent']['med30_aged'])

            output_pesti_params[soil_short]['Err2']['SFmax'].append(pest_data['Error_percent']['low_fresh'])
            output_pesti_params[soil_short]['Err2']['SAmax'].append(pest_data['Error_percent']['low_aged'])

    elif initial_condition == 'min':
        if isLiving and isFirstCycle:
            # Living
            output_pesti_params[soil_short]['Err1']['LFmin'].append(pest_data['Error_percent']['high_fresh'])
            output_pesti_params[soil_short]['Err1']['LAmin'].append(pest_data['Error_percent']['high_aged'])

            output_pesti_params[soil_short]['Err1']['LFmin'].append(pest_data['Error_percent']['med12_fresh'])
            output_pesti_params[soil_short]['Err1']['LAmin'].append(pest_data['Error_percent']['med12_aged'])

            output_pesti_params[soil_short]['Err1']['LFmin'].append(pest_data['Error_percent']['med30_fresh'])
            output_pesti_params[soil_short]['Err1']['LAmin'].append(pest_data['Error_percent']['med30_aged'])

            output_pesti_params[soil_short]['Err1']['LFmin'].append(pest_data['Error_percent']['low_fresh'])
            output_pesti_params[soil_short]['Err1']['LAmin'].append(pest_data['Error_percent']['low_aged'])

        elif not isLiving and isFirstCycle:
            # Sterile
            output_pesti_params[soil_short]['Err1']['SFmin'].append(pest_data['Error_percent']['high_fresh'])
            output_pesti_params[soil_short]['Err1']['SAmin'].append(pest_data['Error_percent']['high_aged'])

            output_pesti_params[soil_short]['Err1']['SFmin'].append(pest_data['Error_percent']['med12_fresh'])
            output_pesti_params[soil_short]['Err1']['SAmin'].append(pest_data['Error_percent']['med12_aged'])

            output_pesti_params[soil_short]['Err1']['SFmin'].append(pest_data['Error_percent']['med30_fresh'])
            output_pesti_params[soil_short]['Err1']['SAmin'].append(pest_data['Error_percent']['med30_aged'])

            output_pesti_params[soil_short]['Err1']['SFmin'].append(pest_data['Error_percent']['low_fresh'])
            output_pesti_params[soil_short]['Err1']['SAmin'].append(pest_data['Error_percent']['low_aged'])

        elif isLiving and not isFirstCycle:
            # Living
            output_pesti_params[soil_short]['Err2']['LFmin'].append(pest_data['Error_percent']['high_fresh'])
            output_pesti_params[soil_short]['Err2']['LAmin'].append(pest_data['Error_percent']['high_aged'])

            output_pesti_params[soil_short]['Err2']['LFmin'].append(pest_data['Error_percent']['med12_fresh'])
            output_pesti_params[soil_short]['Err2']['LAmin'].append(pest_data['Error_percent']['med12_aged'])

            output_pesti_params[soil_short]['Err2']['LFmin'].append(pest_data['Error_percent']['med30_fresh'])
            output_pesti_params[soil_short]['Err2']['LAmin'].append(pest_data['Error_percent']['med30_aged'])

            output_pesti_params[soil_short]['Err2']['LFmin'].append(pest_data['Error_percent']['low_fresh'])
            output_pesti_params[soil_short]['Err2']['LAmin'].append(pest_data['Error_percent']['low_aged'])

        elif not isLiving and not isFirstCycle:
            # Sterile
            output_pesti_params[soil_short]['Err2']['SFmin'].append(pest_data['Error_percent']['high_fresh'])
            output_pesti_params[soil_short]['Err2']['SAmin'].append(pest_data['Error_percent']['high_aged'])

            output_pesti_params[soil_short]['Err2']['SFmin'].append(pest_data['Error_percent']['med12_fresh'])
            output_pesti_params[soil_short]['Err2']['SAmin'].append(pest_data['Error_percent']['med12_aged'])

            output_pesti_params[soil_short]['Err2']['SFmin'].append(pest_data['Error_percent']['med30_fresh'])
            output_pesti_params[soil_short]['Err2']['SAmin'].append(pest_data['Error_percent']['med30_aged'])

            output_pesti_params[soil_short]['Err2']['SFmin'].append(pest_data['Error_percent']['low_fresh'])
            output_pesti_params[soil_short]['Err2']['SAmin'].append(pest_data['Error_percent']['low_aged'])

    return output_pesti_params


def extract_pest(
        pesti_data
):
    # Time axis
    cum_time_30min = pesti_data[:, 0]

    # Cumulative leachate sterilized
    high_0d_cum_mass_out_dt = pesti_data[:, 1]
    high_1d_cum_mass_out_dt = pesti_data[:, 2]

    med12_0d_cum_mass_out_dt = pesti_data[:, 3]
    med12_1d_cum_mass_out_dt = pesti_data[:, 4]

    med30_0d_cum_mass_out_dt = pesti_data[:, 5]
    med30_1d_cum_mass_out_dt = pesti_data[:, 6]

    low_0d_cum_mass_out_dt = pesti_data[:, 7]
    low_1d_cum_mass_out_dt = pesti_data[:, 8]

    # Ponded mass
    high_0d_overmass_dt = pesti_data[:, 9]
    high_1d_overmass_dt = pesti_data[:, 10]
    med12_0d_overmass_dt = pesti_data[:, 11]
    med12_1d_overmass_dt = pesti_data[:, 12]

    med30_0d_overmass_dt = pesti_data[:, 13]
    med30_1d_overmass_dt = pesti_data[:, 14]
    low_0d_overmass_dt = pesti_data[:, 15]
    low_1d_overmass_dt = pesti_data[:, 16]

    mass_percol_condition = stackdata8(
        cum_time_30min,
        high_0d_cum_mass_out_dt, high_1d_cum_mass_out_dt,
        med12_0d_cum_mass_out_dt, med12_1d_cum_mass_out_dt,
        med30_0d_cum_mass_out_dt, med30_1d_cum_mass_out_dt,
        low_0d_cum_mass_out_dt, low_1d_cum_mass_out_dt)

    mass_pond_condition = stackdata8(
        cum_time_30min,
        high_0d_overmass_dt, high_1d_overmass_dt,
        med12_0d_overmass_dt, med12_1d_overmass_dt,
        med30_0d_overmass_dt, med30_1d_overmass_dt,
        low_0d_overmass_dt, low_1d_overmass_dt)

    return [mass_percol_condition, mass_pond_condition]


def choose_pesticide(
        pesticide,
        soil
):
    if soil == 'Alteck':
        fom_untreat = 5.51 / 100.0
        foc_living = 0.58 * fom_untreat
        fom_sterile = 3.87 / 100.0
        foc_sterile = 0.58 * fom_sterile

        # Assign dictionaries
        if pesticide == 'Metalaxyl':
            #  Dictionary contains:
            #  Scenario: (initial_mass, leached_mass_observed, ponded_mass_obs, initial_mass_error, error_leach, error_pond)
            pest_dict_S_1st = {
                'a_high_0d': (1818.12, 138.1, 'nan', 755.5, 8, 0.),
                'b_high_1d': (1472.67, 207.1, 'nan', 631.3, 12, 0.),
                'c_med12_0d': (1818.12, 201.0, 'nan', 755.5, 11, 0.),
                'd_med12_1d': (1472.67, 50.4, 'nan', 631.3, 3, 0.),
                'e_med30_0d': (1818.12, 641.8, 'nan', 755.5, 36, 0.),  # Average obs = 1818.12 (+/- : 755.5)
                'f_med30_1d': (1472.67, 356.8, 'nan', 631.3, 20, 0.),
                'g_low_0d': (1818.12, 177.0, 'nan', 755.5, 10, 0.),
                'h_low_1d': (1472.67, 293.5, 'nan', 631.3, 16, 0.)
            }

            pest_dict_L_1st = {
                'a_high_0d': (1518.06, 145.4, 'nan', 648.3, 8, 0.),
                'b_high_1d': (1413.28, 283.5, 'nan', 597.4, 16, 0.),
                'c_med12_0d': (1518.06, 158.4, 'nan', 648.3, 9, 0.),
                'd_med12_1d': (1413.28, 262.3, 'nan', 597.4, 15, 0.),
                'e_med30_0d': (1518.12, 674.9, 'nan', 648.3, 38, 0.),  # Average obs = 1518.12 (+/- : 648.3)
                'f_med30_1d': (1413.28, 360.2, 'nan', 597.4, 20, 0.),
                'g_low_0d': (1518.06, 418.2, 'nan', 648.3, 23, 0.),
                'h_low_1d': (1413.28, 480.9, 'nan', 597.4, 27, 0.)
            }

            pest_dict_S_2nd = {
                'a_high_0d': (1496.75, 8.35, 5.7, 763.2, 0.5, 0.3),
                'b_high_1d': (1127.52, 37.57, 4.1, 642.9, 2.0, 0.2),
                'c_med12_0d': (1440.72, 290.3, 'nan', 766.8, 'nan', 'nan'),
                'd_med12_1d': (1267.11, 'nan', 9.2, 634.1, 'nan', 0.5),
                'e_med30_0d': (1047.95, 93.3, 4.3, 791.4, 5.2, 0.2),
                'f_med30_1d': (994.09, 82.2, 14.0, 651.3, 4.6, 0.8),
                'g_low_0d': (1462.08, 285.3, 'nan', 765.4, 'nan', 'nan'),
                'h_low_1d': (1050.48, 'nan', 12.4, 647.7, 'nan', 0.7)
            }

            pest_dict_L_2nd = {
                'a_high_0d': (1222.86, 175.44, 4.7, 656.4, 9.8, 0.3),
                'b_high_1d': (1006.54, 40.03, 3.2, 613.3, 1.7, 0.2),
                'c_med12_0d': (1211.28, 272.5, 1.8, 657.1, 15.2, 0.1),
                'd_med12_1d': (1025.43, 168.5, 'nan', 612.1, 9.4, 'nan'),
                'e_med30_0d': (751.13, 35.1, 8.9, 686.0, 2.1, 0.5),
                'f_med30_1d': (938.23, 146.1, 0.1, 617.5, 8.2, 0.008),
                'g_low_0d': (979.82, 86.0, 5.8, 671.7, 4.8, 0.3),
                'h_low_1d': (830.68, 76.5, 9.6, 624.3, 4.3, 0.5)
            }
        elif pesticide == 'S-metolachlor':

            #  Dictionary contains:
            #  Scenario:
            # (initial_mass, leached_mass_observed,
            # ponded_mass_obs, initial_mass_error,
            # error_leach, error_pond)
            pest_dict_S_1st = {
                'a_high_0d': (5176.88, 65.7, 'nan', 1657, 3.7, 'nan'),
                'b_high_1d': (4213.1, 77.7, 'nan', 1400, 4.3, 'nan'),
                'c_med12_0d': (5176.88, 79.5, 'nan', 1657, 4.4, 'nan'),
                'd_med12_1d': (4213.1, 16.9, 'nan', 1400, 0.9, 'nan'),
                'e_med30_0d': (5176.88, 327.4, 'nan', 1657, 18.2, 'nan'),  # 5176.88 +/- 1657
                'f_med30_1d': (4213.1, 153.4, 'nan', 1400, 8.5, 'nan'),
                'g_low_0d': (5176.88, 70.9, 'nan', 1657, 3.9, 'nan'),
                'h_low_1d': (4213.1, 110.7, 'nan', 1400, 6.2, 'nan')
            }
            #  Dictionary contains:
            pest_dict_L_1st = {
                'a_high_0d': (3460.85, 53.5, 'nan', 1135, 2.98, 'nan'),
                'b_high_1d': (2832.66, 89.2, 'nan', 946, 4.97, 'nan'),
                'c_med12_0d': (3460.85, 49.4, 'nan', 1135, 2.75, 'nan'),
                'd_med12_1d': (2832.66, 70.5, 'nan', 946, 3.93, 'nan'),
                'e_med30_0d': (3460.85, 281.7, 'nan', 1135, 15.7, 'nan'),  # 3460.85 +/- 1135
                'f_med30_1d': (2832.66, 116.4, 'nan', 946, 6.49, 'nan'),
                'g_low_0d': (3460.85, 142.7, 'nan', 1135, 7.9, 'nan'),
                'h_low_1d': (2832.66, 127.2, 'nan', 946, 7.08, 'nan')
            }

            pest_dict_S_2nd = {
                'a_high_0d': (3698.64, 2.4, 1.9, 1660, 0.1, 0.1),
                'b_high_1d': (2992.53, 16.8, 1.1, 1404, 1.0, 0.1),
                'c_med12_0d': (3688.61, 82.5, 'nan', 1661, 4.6, 'nan'),
                'd_med12_1d': (3036.55, 'nan', 3.4, 1401, 'nan', 0.2),
                'e_med30_0d': (3509.28, 37.9, 2.4, 1675, 2.1, 0.1),
                'f_med30_1d': (2937.71, 24.3, 8.6, 1408, 1.4, 0.5),
                'g_low_0d': (3694.88, 102.7, 'nan', 1661, 5.7, 'nan'),
                'h_low_1d': (2968.65, 'nan', 5.0, 1406, 'nan', 0.3)
            }

            pest_dict_L_2nd = {
                'a_high_0d': (2466.7, 43.7, 1.4, 1138, 2.44, 0.08),
                'b_high_1d': (1985.26, 19.9, 0.4, 951, 1.1, 0.02),
                'c_med12_0d': (2468.61, 66.8, 'nan', 1137, 3.72, 'nan'),
                'd_med12_1d': (1998.26, 37.9, 'nan', 950, 2.1, 'nan'),
                'e_med30_0d': (2300.56, 16.8, 5.6, 1150, 1., 0.31),
                'f_med30_1d': (1965.55, 46.4, 1.4, 952, 2.59, 0.08),
                'g_low_0d': (2401.14, 22.7, 3.4, 1143, 1.26, 0.19),
                'h_low_1d': (1957.78, 16.7, 4.2, 953, 0.93, 0.24)
            }
        elif pesticide == 'Cupper':
            pest_dict_S_1st = {
                'a_high_0d': (1184.44, 8.6, 'nan', 59, 0.4, 'nan'),
                'b_high_1d': (1177.45, 14.0, 'nan', 59, 0.7, 'nan'),
                'c_med12_0d': (1184.44, 9.2, 'nan', 59, 0.5, 'nan'),
                'd_med12_1d': (1177.45, 14.0, 'nan', 59, 0.7, 'nan'),
                'e_med30_0d': (1184.44, 335.6, 'nan', 59, 16.8, 'nan'),
                'f_med30_1d': (1177.45, 278.0, 'nan', 59, 13.9, 'nan'),
                'g_low_0d': (1184.44, 20.5, 'nan', 59, 1.0, 'nan'),
                'h_low_1d': (1177.45, 21.6, 'nan', 59, 1.1, 'nan')
            }
            pest_dict_L_1st = {
                'a_high_0d': (1626.66, 0.7, 'nan', 81, 0.03, 'nan'),
                'b_high_1d': (1106.95, 'nan', 'nan', 55, 'nan', 'nan'),
                'c_med12_0d': (1626.66, 0.4, 'nan', 81, 0.02, 'nan'),
                'd_med12_1d': (1106.95, 'nan', 'nan', 55, 0.02, 'nan'),
                'e_med30_0d': (1626.66, 40.0, 'nan', 81, 2.0, 'nan'),
                'f_med30_1d': (1106.95, 'nan', 'nan', 55, 1.7, 'nan'),
                'g_low_0d': (1626.66, 1.5, 'nan', 81, 0.07, 'nan'),
                'h_low_1d': (1106.95, 'nan', 'nan', 55, 0.04, 'nan')
            }
            pest_dict_S_2nd = {
                'a_high_0d': (1175.8, 4.85, 'nan', 60, 0.24, 'nan'),
                'b_high_1d': (1163.5, 0.74, 'nan', 59, 0.04, 'nan'),
                'c_med12_0d': (1175.2, 6.42, 'nan', 60, 0.32, 'nan'),
                'd_med12_1d': (1163.4, 'nan', 'nan', 59, 'nan', 'nan'),
                'e_med30_0d': (848.9, 4.84, 1.64, 76, 0.24, 0.08),
                'f_med30_1d': (899.4, 0.94, 4.86, 72, 0.05, 0.24),
                'g_low_0d': (1163.9, 8.61, 'nan', 60, 0.43, 'nan'),
                'h_low_1d': (1155.8, 'nan', 'nan', 60, 'nan', 'nan')
            }
            pest_dict_L_2nd = {
                'a_high_0d': (1626.0, 0.59, 'nan', 81.4, 0.3, 'nan'),
                'b_high_1d': (1106.5, 'nan', 'nan', 55.4, 'nan', 'nan'),
                'c_med12_0d': (1626.3, 1.09, 'nan', 81.4, 0.05, 'nan'),
                'd_med12_1d': (1106.5, 'nan', 'nan', 55.4, 'nan', 'nan'),
                'e_med30_0d': (1586.7, 'nan', 'nan', 81.4, 'nan', 'nan'),
                'f_med30_1d': (1072.9, 'nan', 'nan', 55.4, 'nan', 'nan'),
                'g_low_0d': (1625.2, 0.76, 'nan', 81.4, 0.04, 'nan'),
                'h_low_1d': (1106.1, 'nan', 'nan', 55.4, 'nan', 'nan')
            }
        elif pesticide == 'Zinc':
            pest_dict_S_1st = {
                'a_high_0d': (3106.12, 16.8, 'nan', 155, 0.84, 'nan'),
                'b_high_1d': (2594.16, 16.7, 'nan', 130, 0.83, 'nan'),
                'c_med12_0d': (3106.12, 17.7, 'nan', 155, 0.89, 'nan'),
                'd_med12_1d': (2594.16, 16.3, 'nan', 130, 0.81, 'nan'),
                'e_med30_0d': (3106.12, 523.4, 'nan', 155, 26.17, 'nan'),
                'f_med30_1d': (2594.16, 285.1, 'nan', 130, 14.25, 'nan'),
                'g_low_0d': (3106.12, 37.1, 'nan', 155, 1.85, 'nan'),
                'h_low_1d': (2594.16, 23.2, 'nan', 130, 1.16, 'nan')
            }
            pest_dict_L_1st = {
                'a_high_0d': (2636.02, 0.9, 'nan', 132, 0.04, 'nan'),
                'b_high_1d': (2586.08, 3.5, 'nan', 129, 0.176, 'nan'),
                'c_med12_0d': (2636.02, 0.8, 'nan', 132, 0.04, 'nan'),
                'd_med12_1d': (2586.0, 3.2, 'nan', 129, 0.16, 'nan'),
                'e_med30_0d': (2636.02, 17.2, 'nan', 132, 0.86, 'nan'),
                'f_med30_1d': (2586.0, 56.5, 'nan', 129, 2.83, 'nan'),
                'g_low_0d': (2636.02, 1.7, 'nan', 132, 0.08, 'nan'),
                'h_low_1d': (2586.0, 4.5, 'nan', 129, 0.225, 'nan')
            }
            pest_dict_S_2nd = {
                'a_high_0d': (3089.3, 2.5, 'nan', 156, 0.13, 'nan'),
                'b_high_1d': (2577.5, 0.4, 'nan', 131, 0.02, 'nan'),
                'c_med12_0d': (3088.4, 3.6, 'nan', 156, 0.18, 'nan'),
                'd_med12_1d': (2577.9, 'nan', 'nan', 131, 'nan', 'nan'),
                'e_med30_0d': (2582.7, 1.3, 1.5, 181, 0.07, 0.08),
                'f_med30_1d': (2309.1, 0.2, 2.2, 144, 0.01, 0.11),
                'g_low_0d': (3069.0, 4.1, 'nan', 157, 0.2, 'nan'),
                'h_low_1d': (2571.0, 'nan', 'nan', 131, 'nan', 'nan')
            }
            pest_dict_L_2nd = {
                'a_high_0d': (2635.2, 'nan', 'nan', 132, 'nan', 'nan'),
                'b_high_1d': (2582.6, 0.32, 'nan', 129, 0.02, 'nan'),
                'c_med12_0d': (2635.2, 'nan', 'nan', 132, 'nan', 'nan'),
                'd_med12_1d': (2582.9, 0.99, 'nan', 129, 0.05, 'nan'),
                'e_med30_0d': (2618.8, 'nan', 'nan', 132, 'nan', 'nan'),
                'f_med30_1d': (2529.6, 0.45, 'nan', 132, 0.02, 'nan'),
                'g_low_0d': (2634.3, 'nan', 'nan', 132, 'nan', 'nan'),
                'h_low_1d': (2581.6, 0.05, 'nan', 130, 0.002, 'nan')
            }
        else:
            print("Error assigning dictionaries")

    elif soil == 'Rouff':
        fom_untreat = 2.93 / 100.0
        foc_living = 0.58 * fom_untreat
        fom_sterile = 3.53 / 100.0
        foc_sterile = 0.58 * fom_sterile

        # Assign dictionaries
        if pesticide == 'Metalaxyl':
            pest_dict_S_1st = {
                'a_high_0d': (2105.8, 307.0, 'nan', 634.07, 17.2, 'nan'),
                'b_high_1d': (1504.94, 269.5, 'nan', 462.2, 15.1, 'nan'),
                'c_med12_0d': (2105.8, 279.2, 'nan', 634.07, 15.6, 'nan'),
                'd_med12_1d': (1504.94, 294.9, 'nan', 462.2, 16.5, 'nan'),
                'e_med30_0d': (2105.8, 313.8, 'nan', 634.07, 17.5, 'nan'),
                'f_med30_1d': (1504.94, 529.8, 'nan', 462.2, 29.6, 'nan'),
                'g_low_0d': (2105.8, 419.4, 'nan', 634.07, 23.4, 'nan'),
                'h_low_1d': (1504.94, 346.1, 'nan', 462.2, 19.4, 'nan')
            }

            pest_dict_L_1st = {
                'a_high_0d': (1741.66, 271.3, 'nan', 537.1, 15.2, 'nan'),
                'b_high_1d': (1696.73, 116.9, 'nan', 527.4, 6.5, 'nan'),
                'c_med12_0d': (1741.66, 304.3, 'nan', 537.1, 17.0, 'nan'),
                'd_med12_1d': (1696.73, 276.9, 'nan', 527.4, 15.5, 'nan'),
                'e_med30_0d': (1741.66, 616.4, 'nan', 537.1, 34.5, 'nan'),
                'f_med30_1d': (1696.73, 409.1, 'nan', 527.4, 22.9, 'nan'),
                'g_low_0d': (1741.66, 365.0, 'nan', 537.1, 20.4, 'nan'),
                'h_low_1d': (1696.73, 429.5, 'nan', 527.4, 22.0, 'nan')
            }

            pest_dict_S_2nd = {
                'a_high_0d': (1602.5, 63.8, 8.0, 651, 3.6, 0.5),
                'b_high_1d': (1100.6, 100.4, 6.9, 477, 5.6, 0.4),
                'c_med12_0d': (1627.3, 'nan', 18.8, 650, 'nan', 1.1),
                'd_med12_1d': (1078.0, 142.3, 'nan', 479, 8.0, 'nan'),
                'e_med30_0d': (1596.5, 'nan', 16.8, 652, 'nan', 0.9),
                'f_med30_1d': (868.8, 163.8, 'nan', 492, 9.2, 'nan'),
                'g_low_0d': (1502.4, 130.5, 9.9, 658, 7.3, 0.6),
                'h_low_1d': (1032.4, 222.0, 'nan', 482, 12.4, 'nan')
            }
            pest_dict_L_2nd = {
                'a_high_0d': (1309.9, 'nan', 13.2, 552, 'nan', 0.7),
                'b_high_1d': (1407.45, 83.7, 7.5, 534, 4.7, 0.4),
                'c_med12_0d': (1280.5, 'nan', 10.0, 554, 'nan', 0.6),
                'd_med12_1d': (1264.9, 115.4, 7.1, 543, 6.5, 0.4),
                'e_med30_0d': (1002.5, 'nan', 3.6, 572, 'nan', 0.2),
                'f_med30_1d': (1147.2, 92.4, 3.8, 550, 5.2, 0.2),
                'g_low_0d': (1226.4, 241.0, 0.2, 558, 13.5, 0.01),
                'h_low_1d': (1128.9, 15.3, 4.4, 551, 0.9, 0.2)
            }
        elif pesticide == 'S-metolachlor':
            pest_dict_S_1st = {
                'a_high_0d': (7343.7, 148.5, 'nan', 2771, 8, 'nan'),
                'b_high_1d': (3052.6, 111.5, 'nan', 3053, 6, 'nan'),
                'c_med12_0d': (7343.7, 160.1, 'nan', 2771, 9, 'nan'),
                'd_med12_1d': (3052.6, 110.7, 'nan', 3053, 6, 'nan'),
                'e_med30_0d': (7343.7, 188.3, 'nan', 2771, 10, 'nan'),
                'f_med30_1d': (3052.6, 235.8, 'nan', 3053, 13, 'nan'),
                'g_low_0d': (7343.7, 237.4, 'nan', 2771, 13, 'nan'),
                'h_low_1d': (3052.6, 142.8, 'nan', 3053, 8, 'nan')
            }

            pest_dict_L_1st = {
                'a_high_0d': (4038.5, 111.8, 'nan', 4038, 6, 'nan'),
                'b_high_1d': (4904.5, 37.2, 'nan', 4905, 2, 'nan'),
                'c_med12_0d': (4038.5, 111.3, 'nan', 4038, 6, 'nan'),
                'd_med12_1d': (4904.5, 79.3, 'nan', 4905, 4, 'nan'),
                'e_med30_0d': (4038.5, 278.4, 'nan', 4038, 16, 'nan'),
                'f_med30_1d': (4904.5, 156.1, 'nan', 4905, 9, 'nan'),
                'g_low_0d': (4038.5, 142.0, 'nan', 4038, 8, 'nan'),
                'h_low_1d': (4904.5, 134.7, 'nan', 4905, 8, 'nan')
            }

            pest_dict_S_2nd = {
                'a_high_0d': (5206.7, 'nan', 3.5, 2779, 'nan', 0.2),
                'b_high_1d': (2128.2, 29.1, 3.0, 1145, 2, 'nan'),
                'c_med12_0d': (5198.3, 'nan', 9.6, 2779, 'nan', 0.5),
                'd_med12_1d': (2128.8, 41.1, 'nan', 1145, 2, 'nan'),
                'e_med30_0d': (5177.9, 'nan', 9.8, 2781, 'nan', 0.5),
                'f_med30_1d': (2038.3, 74.7, 'nan', 1152, 4, 'nan'),
                'g_low_0d': (5142.4, 40.5, 4.6, 2784, 2., 0.3),
                'h_low_1d': (2105.6, 70.2, 'nan', 1146, 4, 'nan')
            }

            pest_dict_L_2nd = {
                'a_high_0d': (2841.5, 'nan', 4.2, 1507, 'nan', 0.2),
                'b_high_1d': (3522.2, 16.1, 2.4, 1860, 1, 0.1),
                'c_med12_0d': (2841.9, 'nan', 3.0, 1507, 'nan', 0.2),
                'd_med12_1d': (3491.7, 26.3, 2.3, 1862, 1, 0.1),
                'e_med30_0d': (2720.9, 'nan', 1.5, 1517, 'nan', 0.1),
                'f_med30_1d': (3436.1, 28.5, 2.1, 1867, 2, 0.1),
                'g_low_0d': (2819.6, 77.2, 0.1, 1509, 4, 0.01),
                'h_low_1d': (3451.6, 3.4, 2.5, 1865, 0, 0.1)
            }
        elif pesticide == 'Cupper':
            pest_dict_S_1st = {
                'a_high_0d': (3734.8, 23.4, 'nan', 187, 1, 'nan'),
                'b_high_1d': (4008.4, 38.5, 'nan', 200, 2, 'nan'),
                'c_med12_0d': (3734.8, 28.3, 'nan', 187, 1, 'nan'),
                'd_med12_1d': (4008.4, 29.7, 'nan', 200, 1, 'nan'),
                'e_med30_0d': (3734.8, 611.1, 'nan', 187, 31, 'nan'),
                'f_med30_1d': (4008.4, 588.5, 'nan', 200, 29, 'nan'),
                'g_low_0d': (3734.8, 48.4, 'nan', 187, 2, 'nan'),
                'h_low_1d': (4008.4, 51.8, 'nan', 200, 3, 'nan')
            }

            pest_dict_L_1st = {
                'a_high_0d': (4011.0, 1.4, 'nan', 201, 0.07, 'nan'),
                'b_high_1d': (3859.7, 0.6, 'nan', 193, 0.03, 'nan'),
                'c_med12_0d': (4011.0, 1.0, 'nan', 201, 0.05, 'nan'),
                'd_med12_1d': (3859.7, 0.4, 'nan', 193, 0.02, 'nan'),
                'e_med30_0d': (4011.0, 47.3, 'nan', 201, 2.37, 'nan'),
                'f_med30_1d': (3859.7, 33.3, 'nan', 193, 1.67, 'nan'),
                'g_low_0d': (4011.0, 2.7, 'nan', 201, 0.14, 'nan'),
                'h_low_1d': (3859.7, 0.8, 'nan', 193, 0.04, 'nan')
            }

            pest_dict_S_2nd = {
                'a_high_0d': (3711.5, 16.4, 0.4, 188, 1, 0.02),
                'b_high_1d': (3969.9, 12.9, 'nan', 202, 1, 'nan'),
                'c_med12_0d': (3706.5, 'nan', 1.3, 188, 'nan', 0.07),
                'd_med12_1d': (3978.7, 24.4, 'nan', 202, 1, 'nan'),
                'e_med30_0d': (3123.7, 'nan', 16.9, 217, 'nan', 0.84),
                'f_med30_1d': (3419.9, 412.4, 'nan', 230, 21, 'nan'),
                'g_low_0d': (3686.4, 11.5, 0.7, 203, 1, 0.04),
                'h_low_1d': (3956.6, 35.8, 'nan', 203, 2, 'nan')
            }

            pest_dict_L_2nd = {
                'a_high_0d': (4009.6, 'nan', 1.9, 201, 'nan', 0.09),
                'b_high_1d': (3859.1, 0.3, 'nan', 193, 0.01, 'nan'),
                'c_med12_0d': (4010.0, 'nan', 'nan', 201, 'nan', 'nan'),
                'd_med12_1d': (3859.3, 'nan', 'nan', 193, 'nan', 'nan'),
                'e_med30_0d': (3963.6, 1.9, 'nan', 203, 0.09, 'nan'),
                'f_med30_1d': (3826.3, 'nan', 'nan', 195, 'nan', 'nan'),
                'g_low_0d': (4008.3, 1.7, 'nan', 201, 0.08, 'nan'),
                'h_low_1d': (3858.9, 0.3, 'nan', 193, 0.02, 'nan')
            }
        elif pesticide == 'Zinc':
            pest_dict_S_1st = {
                'a_high_0d': (2810.8, 0.6, 'nan', 141, 0.03, 'nan'),
                'b_high_1d': (3121.5, 2.4, 'nan', 156, 0.12, 'nan'),
                'c_med12_0d': (2810.8, 0.8, 'nan', 141, 0.04, 'nan'),
                'd_med12_1d': (3121.5, 'nan', 'nan', 156, 'nan', 'nan'),
                'e_med30_0d': (2810.8, 24.2, 'nan', 141, 1.21, 'nan'),
                'f_med30_1d': (3121.5, 4.3, 'nan', 156, 0.21, 'nan'),
                'g_low_0d': (2810.8, 1.8, 'nan', 141, 0.09, 'nan'),
                'h_low_1d': (3121.5, 'nan', 'nan', 156, 'nan', 'nan')
            }

            pest_dict_L_1st = {
                'a_high_0d': (3216.2, 'nan', 'nan', 161, 'nan', 'nan'),
                'b_high_1d': (3030.5, 'nan', 'nan', 152, 'nan', 'nan'),
                'c_med12_0d': (3216.2, 'nan', 'nan', 161, 'nan', 'nan'),
                'd_med12_1d': (3030.5, 'nan', 'nan', 152, 'nan', 'nan'),
                'e_med30_0d': (3216.2, 0.5, 'nan', 161, 0.02, 'nan'),
                'f_med30_1d': (3030.5, 10.0, 'nan', 152, 0.5, 'nan'),
                'g_low_0d': (3216.2, 'nan', 'nan', 161, 'nan', 'nan'),
                'h_low_1d': (3030.5, 'nan', 'nan', 152, 'nan', 'nan')
            }

            pest_dict_S_2nd = {
                'a_high_0d': (2810.2, 'nan', 'nan', 141, 'nan', 'nan'),
                'b_high_1d': (3119.1, 'nan', 'nan', 156, 'nan', 'nan'),
                'c_med12_0d': (2810.0, 'nan', 'nan', 141, 'nan', 'nan'),
                'd_med12_1d': (3022, 'nan', 'nan', 'nan', 'nan', 'nan'),  # not recorded
                'e_med30_0d': (2786.6, 'nan', 'nan', 142, 'nan', 'nan'),
                'f_med30_1d': (3117.3, 1.6, 'nan', 156, 0.08, 'nan'),
                'g_low_0d': (2809.0, 'nan', 'nan', 141, 'nan', 'nan'),
                'h_low_1d': (3010, 'nan', 'nan', 'nan', 'nan', 'nan')  # not recorded
            }

            pest_dict_L_2nd = {
                'a_high_0d': (3264, 'nan', 'nan', 'nan', 'nan', 'nan'),  # not recorded
                'b_high_1d': (2875, 'nan', 'nan', 'nan', 'nan', 'nan'),  # not recorded
                'c_med12_0d': (3428, 'nan', 'nan', 'nan', 'nan', 'nan'),  # not recorded
                'd_med12_1d': (2890, 'nan', 'nan', 'nan', 'nan', 'nan'),  # not recorded
                'e_med30_0d': (3215.7, 'nan', 'nan', 161, 'nan', 'nan'),
                'f_med30_1d': (3020.6, 'nan', 'nan', 152, 'nan', 'nan'),  # not recorded
                'g_low_0d': (2196, 'nan', 'nan', 'nan', 'nan', 'nan'),  # not recorded
                'h_low_1d': (2860, 'nan', 'nan', 'nan', 'nan', 'nan')  # not recorded
            }
        else:
            print("Error assigning dictionaries")

    # Compute Kd's (foc calculted above)
    if pesticide == 'Metalaxyl':
        # Pesticide Koc
        Koc_mexyl = [100, 75, 50, 30, 28., 25, 17., 15.0, 14., 10.]  # [(a) , (b), (c)] [ml/g]
        # Koc_mexyl = [30, 25, 20, 18.0, 16., 10]
        # Koc_mexyl = [50, 30, 20, 10, 2, 1]
        Koc_mexyl = np.array(Koc_mexyl) * 10 ** 3  # [mm3/g]

        # Kd (a) - NPIC @ http://npic.orst.edu/ingred/ppdmove.htm
        Kd_mexyl1_crop_sterile = Koc_mexyl[0] * foc_sterile
        Kd_mexyl1_crop_untreat = Koc_mexyl[0] * foc_living

        # Kd (b) - PAN @ http://www.pesticideinfo.org/
        Kd_mexyl2_crop_sterile = Koc_mexyl[1] * foc_sterile
        Kd_mexyl2_crop_untreat = Koc_mexyl[1] * foc_living

        # Kd (c) - https://toxnet.nlm.nih.gov/cgi-bin/sis/search/a?dbs+hsdb:@term+@DOCNO+7061
        Kd_mexyl3_crop_sterile = Koc_mexyl[2] * foc_sterile
        Kd_mexyl3_crop_untreat = Koc_mexyl[2] * foc_living

        Kd_mexyl4_crop_sterile = Koc_mexyl[3] * foc_sterile
        Kd_mexyl4_crop_untreat = Koc_mexyl[3] * foc_living

        Kd_mexyl5_crop_sterile = Koc_mexyl[4] * foc_sterile
        Kd_mexyl5_crop_untreat = Koc_mexyl[4] * foc_living

        Kd_mexyl6_crop_sterile = Koc_mexyl[5] * foc_sterile
        Kd_mexyl6_crop_untreat = Koc_mexyl[5] * foc_living

        Kd_mexyl7_crop_sterile = Koc_mexyl[6] * foc_sterile
        Kd_mexyl7_crop_untreat = Koc_mexyl[6] * foc_living

        Kd_mexyl8_crop_sterile = Koc_mexyl[7] * foc_sterile
        Kd_mexyl8_crop_untreat = Koc_mexyl[7] * foc_living

        Kd_mexyl9_crop_sterile = Koc_mexyl[8] * foc_sterile
        Kd_mexyl9_crop_untreat = Koc_mexyl[8] * foc_living

        Kd_pest_sterile = [Kd_mexyl1_crop_sterile,
                           Kd_mexyl2_crop_sterile,
                           Kd_mexyl3_crop_sterile,
                           Kd_mexyl4_crop_sterile,
                           Kd_mexyl5_crop_sterile,
                           Kd_mexyl6_crop_sterile,
                           Kd_mexyl7_crop_sterile,
                           Kd_mexyl8_crop_sterile,
                           Kd_mexyl9_crop_sterile]

        Kd_pest_living = [Kd_mexyl1_crop_untreat,
                          Kd_mexyl2_crop_untreat,
                          Kd_mexyl3_crop_untreat,
                          Kd_mexyl4_crop_untreat,
                          Kd_mexyl5_crop_untreat,
                          Kd_mexyl6_crop_untreat,
                          Kd_mexyl7_crop_untreat,
                          Kd_mexyl8_crop_untreat,
                          Kd_mexyl9_crop_untreat]

    elif pesticide == 'S-metolachlor':
        # Koc Ranges (EU COmmision SANCO/1426/2001 - rev. 3. 4 October 2004)
        # Koc_smeto = [369, 200, 110, 70, 65, 21]  # [ml/g]
        Koc_smeto = [369, 300., 250, 225, 200, 150, 125, 100, 75, 50, 30]  # [ml/g]
        Koc_smeto = np.array(Koc_smeto) * 10 ** 3

        # Kd (S-metolachlor)
        Kd_smeto_crop_sterile1 = Koc_smeto[0] * foc_sterile  # ml/g
        Kd_smeto_crop_untreat1 = Koc_smeto[0] * foc_living
        Kd_smeto_crop_sterile2 = Koc_smeto[1] * foc_sterile  # ml/g
        Kd_smeto_crop_untreat2 = Koc_smeto[1] * foc_living
        Kd_smeto_crop_sterile3 = Koc_smeto[2] * foc_sterile  # ml/g
        Kd_smeto_crop_untreat3 = Koc_smeto[2] * foc_living
        Kd_smeto_crop_sterile4 = Koc_smeto[3] * foc_sterile  # ml/g
        Kd_smeto_crop_untreat4 = Koc_smeto[3] * foc_living
        Kd_smeto_crop_sterile5 = Koc_smeto[4] * foc_sterile  # ml/g
        Kd_smeto_crop_untreat5 = Koc_smeto[4] * foc_living
        Kd_smeto_crop_sterile6 = Koc_smeto[5] * foc_sterile  # ml/g
        Kd_smeto_crop_untreat6 = Koc_smeto[5] * foc_living
        Kd_smeto_crop_sterile7 = Koc_smeto[6] * foc_sterile  # ml/g
        Kd_smeto_crop_untreat7 = Koc_smeto[6] * foc_living
        Kd_smeto_crop_sterile8 = Koc_smeto[7] * foc_sterile  # ml/g
        Kd_smeto_crop_untreat8 = Koc_smeto[7] * foc_living
        Kd_smeto_crop_sterile9 = Koc_smeto[8] * foc_sterile  # ml/g
        Kd_smeto_crop_untreat9 = Koc_smeto[8] * foc_living
        Kd_smeto_crop_sterile10 = Koc_smeto[9] * foc_sterile  # ml/g
        Kd_smeto_crop_untreat10 = Koc_smeto[9] * foc_living
        Kd_smeto_crop_sterile11 = Koc_smeto[10] * foc_sterile  # ml/g
        Kd_smeto_crop_untreat11 = Koc_smeto[10] * foc_living

        Kd_pest_sterile = [Kd_smeto_crop_sterile1,
                           Kd_smeto_crop_sterile2,
                           Kd_smeto_crop_sterile3,
                           Kd_smeto_crop_sterile4,
                           Kd_smeto_crop_sterile5,
                           Kd_smeto_crop_sterile6,
                           Kd_smeto_crop_sterile7,
                           Kd_smeto_crop_sterile8,
                           Kd_smeto_crop_sterile9,
                           Kd_smeto_crop_sterile10,
                           Kd_smeto_crop_sterile11]

        Kd_pest_living = [Kd_smeto_crop_untreat1,
                          Kd_smeto_crop_untreat2,
                          Kd_smeto_crop_untreat3,
                          Kd_smeto_crop_untreat4,
                          Kd_smeto_crop_untreat5,
                          Kd_smeto_crop_untreat6,
                          Kd_smeto_crop_untreat7,
                          Kd_smeto_crop_untreat8,
                          Kd_smeto_crop_untreat9,
                          Kd_smeto_crop_untreat10,
                          Kd_smeto_crop_untreat11]

    elif pesticide == 'Cupper':
        # Kd (Copper). Allison and Allison, 2005 - EPA/600/R-05/074:
        # log(Kd) range: 0.1 - 7.0, max-mean = 5.5
        Kd_copper1 = 10 ** 0.1  # [mL/g] = [L/Kg]
        Kd_copper2 = 10 ** 0.5  # [mL/g] = [L/Kg]
        Kd_copper3 = 10 ** 0.6  # [mL/g] = [L/Kg]
        Kd_copper4 = 10 ** 0.9  # [mL/g] = [L/Kg]
        Kd_copper5 = 10 ** 1.2  # [mL/g] = [L/Kg]
        Kd_copper6 = 10 ** 1.7  # [mL/g] = [L/Kg]
        Kd_copper7 = 10 ** 2.7  # [mL/g] = [L/Kg]
        Kd_copper8 = 10 ** 5.7  # [mL/g] = [L/Kg]

        Kd_copper = [Kd_copper1, Kd_copper2, Kd_copper3, Kd_copper4,
                     Kd_copper5, Kd_copper6, Kd_copper7, Kd_copper8]
        Kd_pest_living = np.array(Kd_copper) * 10 ** 3
        Kd_pest_sterile = np.array(Kd_copper) * 10 ** 3

    elif pesticide == 'Zinc':
        # Kd (Zinc) Allison and Allison, 2005 - EPA/600/R-05/074:
        # log(Kd) range: 1.5 - 6.9
        Kd_zinc1 = 10 ** 0.5  # [mL/g] = [L/Kg]
        Kd_zinc2 = 10 ** 0.8
        Kd_zinc3 = 10 ** 1.0
        Kd_zinc4 = 10 ** 1.2
        Kd_zinc5 = 10 ** 1.5
        Kd_zinc6 = 10 ** 3.5
        Kd_zinc7 = 10 ** 5.5
        Kd_zinc8 = 10 ** 6.9

        Kd_zinc = [Kd_zinc1, Kd_zinc2, Kd_zinc3, Kd_zinc4,
                   Kd_zinc5, Kd_zinc6, Kd_zinc7, Kd_zinc8]
        Kd_pest_living = np.array(Kd_zinc) * 10 ** 3  # mm3/g
        Kd_pest_sterile = np.array(Kd_zinc) * 10 ** 3

    else:
        print("Error assigning Kd's to test")

    return [Kd_pest_sterile, Kd_pest_living,
            pest_dict_S_1st, pest_dict_L_1st, pest_dict_S_2nd, pest_dict_L_2nd]



