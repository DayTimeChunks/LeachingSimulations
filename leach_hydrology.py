from greenampt import *
from hydroplots import *


def leachsim(kSat,  # 2.25 mm/min (13.5 cm/h - Crop Soil) - Alteck
             soil_height, # mm
             soil,
             dtGA = 1,  # Timestep in minute
             StormD = 30,  # Storm duration in min
             intensityM = [2.25, 0.92, 0.5],  # mm/min = [0.225, 0.092, 0.05] cm/min
             numSystems = 3,  # Number of combined replicas
             psi=1100,  # soil suction Alteck mm
             ):
    # Saturated water content (0.61 - Alteck; 0.55 - Rouff)
    if soil == 'Alteck':
        ovSat = 0.61
        ov = 0.20
    elif soil == 'Rouff':
        ovSat = 0.55
        ov = 0.20

    wfs = psi * (ovSat - ov)  # Numerator in Green-Ampt's equation = (suction x distance to saturation)
    """ Microcosm Dimensions """
    zl = soil_height  # Mixing layer depth
    d = (14.93 * 2)  # Diameter of falcon tube (mm)
    area = ((d / 2) ** 2) * 3.1416

    loop = 0
    # Running each of the intesities in one loop
    for i in range(len(intensityM)):
        loop += 1

        """ Initialization """
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
        time_size_dt = []

        """ Start of rainfall event X """
        while cum_time < StormD:
            cum_time += dtGA
            intensity = intensityM[i]
            cum_precip += intensity*dtGA

            """
            Record hydrological data:
            All events are simulated for 30 minutes;
            Simulation result is checked at:
                6min = highest intensity
                12min = medium intesity
                30min = lowest intensity
            """
            green_output = greenampt(ksat=kSat,
                                     wfs=wfs,
                                     dtga=dtGA,
                                     intensity=intensity,
                                     ovsat=ovSat,
                                     ov=ov,
                                     zl=zl,
                                     cum_infilt=cum_inf,
                                     cum_time=cum_time)

            # GA output in (mm)
            cum_inf = green_output[0]
            inf = green_output[1]
            runoff = green_output[2]
            leach = green_output[3]
            time_size = green_output[4]

            # Instantaneous time series
            infil_vol = (inf * area * numSystems)  # mm3
            infil_dt.append(infil_vol)

            runoff_vol = (runoff * area * numSystems)
            runoff_dt.append(runoff_vol)

            leached_vol = (leach * area * numSystems)
            leach_dt.append(leached_vol)

            # Variable time step size (due to ponding in-between dtga)
            time_size_dt.append(time_size)

            # Cummulative time series
            cum_time_dt.append(cum_time)

            precip_vol = (cum_precip * area * numSystems)  # mm3
            cum_precip_dt.append(precip_vol)

            cum_inf_vol = (cum_inf * area * numSystems)
            cum_inf_dt.append(cum_inf_vol)

            cum_runoff += runoff
            cum_runoff_vol = (cum_runoff * area * numSystems)
            cum_runoff_dt.append(cum_runoff_vol)

            cum_leach += leach
            cum_leach_vol = (cum_leach * area * numSystems)
            cum_leach_dt.append(cum_leach_vol)

        # Store hydro. variables  for intesity: 135 mm/h
        if loop == 1:
            cum_time_6min = cum_time_dt
            # print("cum time 6 min: ", cum_time_6min[5])
            # print("cum time 12 min: ", cum_time_6min[11])
            # print("cum time end: ", cum_time_6min[-1])
            # print("cum precip (3 systems): ", cum_precip_dt[-1])
            cum_precip_6min = cum_precip_dt
            cum_runoff_6min = cum_runoff_dt
            cum_leach_6min = cum_leach_dt
            # print(cum_leach_6min)
            cum_inf_6min = cum_inf_dt
            runoff_6min = runoff_dt
            # print("cum runoff end: ", cum_runoff_dt[-1])
            infil_6min = infil_dt
            # print("cum inf. end: ", cum_inf_dt[-1])
            leach_6min = leach_dt
            # print("cum leach end: ", cum_leach_dt[-1])
            time_size_6min = time_size_dt
            mass_bal = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
            print("Mass balance", abs(mass_bal) < 10**(-6))

        # Store hydro. variables  for intesity: 55 mm/h
        elif loop == 2:
            cum_time_12min = cum_time_dt
            cum_precip_12min = cum_precip_dt
            cum_runoff_12min = cum_runoff_dt
            cum_leach_12min = cum_leach_dt
            cum_inf_12min = cum_inf_dt
            runoff_12min = runoff_dt
            infil_12min = infil_dt
            leach_12min = leach_dt
            time_size_12min = time_size_dt
            mass_bal = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
            print("Mass balance", mass_bal)

        # Store hydro. variables  for intesity: 30 mm/h
        elif loop == 3:
            cum_time_30min = cum_time_dt
            cum_precip_30min = cum_precip_dt
            cum_runoff_30min = cum_runoff_dt
            cum_leach_30min = cum_leach_dt
            cum_inf_30min = cum_inf_dt
            runoff_30min = runoff_dt
            infil_30min = infil_dt
            leach_30min = leach_dt
            time_size_30min = time_size_dt
            mass_bal = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
            print("Mass balance", mass_bal)

        else:
            print("Do you have more than 3 rainfall scenarios?")

    # Data in mm3
    return stackdata21(cum_time_30min,
                       leach_6min, leach_12min, leach_30min,
                       cum_inf_6min, cum_inf_12min, cum_inf_30min,
                       cum_leach_6min, cum_leach_12min, cum_leach_30min,
                       runoff_6min, runoff_12min, runoff_30min,
                       cum_runoff_6min, cum_runoff_12min, cum_runoff_30min,
                       infil_6min, infil_12min, infil_30min,
                       time_size_6min, time_size_12min, time_size_30min)


def leachsim2(
        high_leach_obs, med12_leach_obs, med30_leach_obs, low_leach_obs,  # each as list
        kSat,  # mm/min (13.5 cm/h - Crop Soil) - Alteck
        soil_height,  # mm
        soil,
        isAGED,
        isFirstCycle,
        dtGA = 1,  # Timestep in minute
        StormD = 30,  # Storm duration in min
        intensityM = [2.25, 0.92, 0.92, 0.5],  # mm/min = [0.225, 0.092, 0.05] cm/min
        numSystems = 3,  # Number of combined replicas
        psi=1100,  # soil suction Alteck mm
        ):
    """
    finds optimal Ksat for AGED or FRESH Soils only.
    """

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

    wfs = psi * (ovSat - ov)  # Numerator in Green-Ampt's equation = (suction x distance to saturation)
    """ Microcosm Dimensions """
    zl = soil_height  # Mixing layer depth
    d = (14.93 * 2)  # Diameter of falcon tube (mm)
    area = ((d / 2) ** 2) * 3.1416

    intensity_num = 0
    # Running each of the intesities in one loop

    for i in range(len(intensityM)):
        intensity_num += 1
        error = 10**9

        # Test best kSats in range
        for k in kSat:

            """ Initialization """
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
            time_size_dt = []

            """ Start of rainfall event X """
            while cum_time < StormD:
                cum_time += dtGA
                intensity = intensityM[i]
                cum_precip += intensity*dtGA

                """
                Record hydrological data:
                All events are simulated for 30 minutes;
                Simulation result is checked at:
                    6min = highest intensity
                    12min = medium intesity
                    30min = lowest intensity
                """
                green_output = greenampt(ksat=k,
                                         wfs=wfs,
                                         dtga=dtGA,
                                         intensity=intensity,
                                         ovsat=ovSat,
                                         ov=ov,
                                         zl=zl,
                                         cum_infilt=cum_inf,
                                         cum_time=cum_time)

                # GA output in (mm)
                cum_inf = green_output[0]
                inf = green_output[1]
                runoff = green_output[2]
                leach = green_output[3]
                time_size = green_output[4]

                # Instantaneous time series
                infil_vol = (inf * area * numSystems)  # mm3
                infil_dt.append(infil_vol)

                runoff_vol = (runoff * area * numSystems)
                runoff_dt.append(runoff_vol)

                leached_vol = (leach * area * numSystems)
                leach_dt.append(leached_vol)

                # Variable time step size (due to ponding in-between dtga)
                time_size_dt.append(time_size)

                # Cummulative time series
                cum_time_dt.append(cum_time)

                precip_vol = (cum_precip * area * numSystems)  # mm3
                cum_precip_dt.append(precip_vol)

                cum_inf_vol = (cum_inf * area * numSystems)
                cum_inf_dt.append(cum_inf_vol)

                cum_runoff += runoff
                cum_runoff_vol = (cum_runoff * area * numSystems)
                cum_runoff_dt.append(cum_runoff_vol)

                cum_leach += leach
                cum_leach_vol = (cum_leach * area * numSystems)
                cum_leach_dt.append(cum_leach_vol)

            if not isAGED:
                age = 0
            else:
                age = 2

            if intensity_num == 1:
                ksat_error = ((cum_leach_dt[5]/10**3 - high_leach_obs[age]) ** 2 +
                              (cum_leach_dt[5]/10**3 - high_leach_obs[age + 1]) ** 2)/2

            elif intensity_num == 2:
                ksat_error = ((cum_leach_dt[11]/10**3 - med12_leach_obs[age]) ** 2 +
                              (cum_leach_dt[11]/10**3 - med12_leach_obs[age + 1]) ** 2) / 2

                # print("intensity 2nd, error", ksat_error, "sim value at 12: ", cum_leach_dt[11]/10**3)
            elif intensity_num == 3:
                ksat_error = ((cum_leach_dt[-1]/10**3 - med30_leach_obs[age]) ** 2 +
                              (cum_leach_dt[-1]/10**3 - med30_leach_obs[age + 1]) ** 2) / 2

            elif intensity_num == 4:
                ksat_error = ((cum_leach_dt[-1]/10**3 - low_leach_obs[age]) ** 2 +
                              (cum_leach_dt[-1]/10**3 - low_leach_obs[age + 1]) ** 2) / 2

            if ksat_error < error:
                error = ksat_error
                # Store hydro. variables  for intesity: 135 mm/h
                if intensity_num == 1:
                    cum_time_6min = cum_time_dt
                    # print("cum time 6 min: ", cum_time_6min[5])
                    # print("cum time 12 min: ", cum_time_6min[11])
                    # print("cum time end: ", cum_time_6min[-1])
                    # print("cum precip (3 systems): ", cum_precip_dt[-1])
                    cum_precip_6min = cum_precip_dt
                    cum_runoff_6min = cum_runoff_dt
                    cum_leach_6min = cum_leach_dt
                    cum_inf_6min = cum_inf_dt
                    runoff_6min = runoff_dt
                    # print("cum runoff end: ", cum_runoff_dt[-1])
                    infil_6min = infil_dt
                    # print("cum inf. end: ", cum_inf_dt[-1])
                    leach_6min = leach_dt
                    # print("cum leach end: ", cum_leach_dt[-1])
                    time_size_6min = time_size_dt
                    mass_bal_1 = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                    k_high = k
                    high_error_prc = (((cum_leach_dt[5] / 10 ** 3 - high_leach_obs[age]) / high_leach_obs[age] +
                                       (cum_leach_dt[5] / 10 ** 3 - high_leach_obs[age + 1]) / high_leach_obs[
                                           age + 1]) / 2) * 100

                    # print("Sim:", cum_leach_dt[5] / 10 ** 3)
                    # print("Obs:", high_leach_obs[age], high_leach_obs[age+1])
                    # print(cum_leach_dt)
                    SS1 = ((cum_leach_dt[5] / 10 ** 3 - high_leach_obs[age]) ** 2 +
                           (cum_leach_dt[5] / 10 ** 3 - high_leach_obs[age + 1]) ** 2)

                # Store hydro. variables  for intesity: 55 mm/h
                elif intensity_num == 2:
                    cum_time_med_12min = cum_time_dt
                    cum_precip_med_12min = cum_precip_dt
                    cum_runoff_med_12min = cum_runoff_dt
                    cum_leach_med_12min = cum_leach_dt
                    cum_inf_med_12min = cum_inf_dt
                    runoff_med_12min = runoff_dt
                    infil_med_12min = infil_dt
                    leach_med_12min = leach_dt
                    time_size_med_12min = time_size_dt
                    mass_bal_2 = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                    k_med12 = k
                    med12_error_prc = ((((cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[age]) / med12_leach_obs[age]) +
                                        ((cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[age + 1]) / med12_leach_obs[
                                            age + 1])
                                        ) / 2
                                       ) * 100
                    SS2 = ((cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[age]) ** 2 +
                           (cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[age + 1]) ** 2)

                elif intensity_num == 3:
                    cum_time_med_30min = cum_time_dt
                    cum_precip_med_30min = cum_precip_dt
                    cum_runoff_med_30min = cum_runoff_dt
                    cum_leach_med_30min = cum_leach_dt
                    cum_inf_med_30min = cum_inf_dt
                    runoff_med_30min = runoff_dt
                    infil_med_30min = infil_dt
                    leach_med_30min = leach_dt
                    time_size_med_30min = time_size_dt
                    mass_bal_3 = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                    k_med30 = k

                    med30_error_prc = (((cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[age]) / med30_leach_obs[age] +
                                        (cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[age + 1]) / med30_leach_obs[
                                            age + 1]) / 2) * 100
                    SS3 = ((cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[age]) ** 2 +
                           (cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[age + 1]) ** 2)

                # Store hydro. variables  for intesity: 30 mm/h
                elif intensity_num == 4:
                    cum_time_30min = cum_time_dt
                    cum_precip_low_30min = cum_precip_dt
                    cum_runoff_low_30min = cum_runoff_dt
                    cum_leach_low_30min = cum_leach_dt
                    cum_inf_low_30min = cum_inf_dt
                    runoff_low_30min = runoff_dt
                    infil_low_30min = infil_dt
                    leach_low_30min = leach_dt
                    time_size_low_30min = time_size_dt
                    mass_bal_4 = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                    k_low = k
                    low_error_prc = (((cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[age]) / low_leach_obs[age] +
                                      (cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[age + 1]) / low_leach_obs[
                                          age + 1]) / 2) * 100
                    SS4 = ((cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[age]) ** 2 +
                           (cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[age + 1]) ** 2)

                else:
                    print("Do you have more than 3 rainfall scenarios?")

            else:
                continue
    mean_obs = (np.mean(high_leach_obs) + np.mean(med12_leach_obs) + np.mean(med30_leach_obs) + np.mean(low_leach_obs))/4
    SStot = 0
    for i in high_leach_obs:
        SStot += (i-mean_obs)**2
    for i in med12_leach_obs:
        SStot += (i-mean_obs)**2
    for i in med30_leach_obs:
        SStot += (i-mean_obs)**2
    for i in low_leach_obs:
        SStot += (i-mean_obs)**2

    SSres = SS1 + SS2 + SS3 + SS4
    r_squared = 1 - (SSres / SStot)

    if isAGED:
        print("AGED")
    else:
        print("FRESH")
    print("--------------------------------------------")
    print("ksat high: ", k_high/10*60, "cm/h")
    print("ksat med12: ", k_med12/10*60, "cm/h")
    print("ksat med30: ", k_med30/10*60, "cm/h")
    print("ksat low: ", k_low/10*60, "cm/h")
    print("--------------------------------------------")
    print("R2: ", r_squared)
    print("--------------------------------------------")
    print("Simulation error percent (%), by modality")
    print("--------------------------------------------")
    print("135 mm/h - 6min ", high_error_prc)
    print("55 mm/h - 12min ", med12_error_prc)
    print("55 mm/h - 30min ", med30_error_prc)
    print("30 mm/h - 30min ", low_error_prc)
    print("--------------------------------------------")
    print("Mass balance: ", abs(mass_bal_1) < 10 ** (-6), abs(mass_bal_2) < 10 ** (-6), abs(mass_bal_3) < 10 ** (-6),
          abs(mass_bal_4) < 10 ** (-6))

    # Data in mm3
    return stackdata28(cum_time_30min,
                       leach_6min, leach_med_12min, leach_med_30min, leach_low_30min,
                       cum_inf_6min, cum_inf_med_12min, cum_inf_med_30min, cum_inf_low_30min,
                       cum_leach_6min, cum_leach_med_12min, cum_leach_med_30min, cum_leach_low_30min,
                       runoff_6min, runoff_med_12min, runoff_med_30min, runoff_low_30min,
                       cum_runoff_6min, cum_runoff_med_12min, cum_runoff_med_30min, cum_runoff_low_30min,
                       infil_6min, infil_med_12min, infil_med_30min, infil_low_30min,
                       time_size_6min, time_size_med_12min, time_size_med_30min, time_size_low_30min)


def leachsim3(
        high_leach_obs, med12_leach_obs, med30_leach_obs, low_leach_obs,  # each as list
        kSat,  # [list of Ksats to test]
        soil_height,  # mm
        soil,
        isFirstCycle,
        dtGA = 1,  # Timestep in minute
        StormD = 30,  # Storm duration in min
        intensityM = [2.25, 0.92, 0.92, 0.5],  # mm/min = [0.225, 0.092, 0.05] cm/min
        numSystems = 3,  # Number of combined replicas
        psi=1100,  # soil suction Alteck mm
        ):
    """
    Finds optimal Ksats for each soil.
    """

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

    wfs = psi * (ovSat - ov)  # Numerator in Green-Ampt's equation = (suction x distance to saturation)
    """ Microcosm Dimensions """
    zl = soil_height  # Mixing layer depth
    d = (14.93 * 2)  # Diameter of falcon tube (mm)
    area = ((d / 2) ** 2) * 3.1416

    # Running each of the intesities in one loop
    for soil in range(4):  # [sterile, untreat, sterile_aged, untreat_aged]
        for ii in range(len(intensityM)):
            error = 10 ** 9
            for k in kSat:  # Test best kSat in list
                """ Initialization """
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
                time_size_dt = []

                """ Start of rainfall event X """
                while cum_time < StormD:
                    cum_time += dtGA
                    intensity = intensityM[ii]
                    cum_precip += intensity*dtGA

                    """
                    Record hydrological data:
                    All events are simulated for 30 minutes;
                    Simulation result is checked at:
                        6min = highest intensity
                        12min = medium intesity
                        30min = lowest intensity
                    """
                    green_output = greenampt(ksat=k,
                                             wfs=wfs,
                                             dtga=dtGA,
                                             intensity=intensity,
                                             ovsat=ovSat,
                                             ov=ov,
                                             zl=zl,
                                             cum_infilt=cum_inf,
                                             cum_time=cum_time)

                    # GA output in (mm)
                    cum_inf = green_output[0]
                    inf = green_output[1]
                    runoff = green_output[2]
                    leach = green_output[3]
                    time_size = green_output[4]

                    # Instantaneous time series
                    infil_vol = (inf * area * numSystems)  # mm3
                    infil_dt.append(infil_vol)

                    runoff_vol = (runoff * area * numSystems)
                    runoff_dt.append(runoff_vol)

                    leached_vol = (leach * area * numSystems)
                    leach_dt.append(leached_vol)

                    # Variable time step size (due to ponding in-between dtga)
                    time_size_dt.append(time_size)

                    # Cummulative time series
                    cum_time_dt.append(cum_time)

                    precip_vol = (cum_precip * area * numSystems)  # mm3
                    cum_precip_dt.append(precip_vol)

                    cum_inf_vol = (cum_inf * area * numSystems)
                    cum_inf_dt.append(cum_inf_vol)

                    cum_runoff += runoff
                    cum_runoff_vol = (cum_runoff * area * numSystems)
                    cum_runoff_dt.append(cum_runoff_vol)

                    cum_leach += leach
                    cum_leach_vol = (cum_leach * area * numSystems)
                    cum_leach_dt.append(cum_leach_vol)

                if ii == 0:
                    ksat_error = (cum_leach_dt[5]/10**3 - high_leach_obs[soil]) ** 2

                elif ii == 1:
                    ksat_error = (cum_leach_dt[11]/10**3 - med12_leach_obs[soil]) ** 2

                elif ii == 2:
                    ksat_error = (cum_leach_dt[-1]/10**3 - med30_leach_obs[soil]) ** 2

                elif ii == 3:
                    ksat_error = (cum_leach_dt[-1]/10**3 - low_leach_obs[soil]) ** 2

                if ksat_error < error:
                    error = ksat_error
                    # s index -> [sterile, untreat, sterile_aged, untreat_aged] -> [SF, LF, SA, LA]
                    ####################################################
                    # Store hydro. variables  for intesity: 135 mm/h
                    if ii == 0 and soil == 0:
                        cum_time_6min_SF = cum_time_dt
                        cum_precip_6min_SF = cum_precip_dt
                        cum_runoff_6min_SF = cum_runoff_dt
                        cum_leach_6min_SF = cum_leach_dt
                        cum_inf_6min_SF = cum_inf_dt
                        runoff_6min_SF = runoff_dt
                        # print("cum runoff end: ", cum_runoff_dt[-1])
                        infil_6min_SF = infil_dt
                        # print("cum inf. end: ", cum_inf_dt[-1])
                        leach_6min_SF = leach_dt
                        # print("cum leach end: ", cum_leach_dt[-1])
                        time_size_6min_SF = time_size_dt
                        mass_bal_1_SF = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_high_SF = k
                        high_SF_error_prc = ((cum_leach_dt[5] / 10 ** 3 - high_leach_obs[soil]) / high_leach_obs[soil]) * 100

                        # print("Sim:", cum_leach_dt[5] / 10 ** 3)
                        # print("Obs:", high_leach_obs[s])
                        # print(cum_leach_dt)
                        SS1_SF = (cum_leach_dt[5] / 10 ** 3 - high_leach_obs[soil]) ** 2

                    elif ii == 0 and soil == 1:
                        cum_time_6min_LF = cum_time_dt
                        cum_precip_6min_LF = cum_precip_dt
                        cum_runoff_6min_LF = cum_runoff_dt
                        cum_leach_6min_LF = cum_leach_dt
                        cum_inf_6min_LF = cum_inf_dt
                        runoff_6min_LF = runoff_dt
                        # print("cum runoff end: ", cum_runoff_dt[-1])
                        infil_6min_LF = infil_dt
                        # print("cum inf. end: ", cum_inf_dt[-1])
                        leach_6min_LF = leach_dt
                        # print("cum leach end: ", cum_leach_dt[-1])
                        time_size_6min_LF = time_size_dt
                        mass_bal_1_LF = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_high_LF = k
                        high_LF_error_prc = ((cum_leach_dt[5] / 10 ** 3 - high_leach_obs[soil]) / high_leach_obs[soil]) * 100

                        # print("Sim:", cum_leach_dt[5] / 10 ** 3)
                        # print("Obs:", high_leach_obs[s])
                        # print(cum_leach_dt)
                        SS1_LF = (cum_leach_dt[5] / 10 ** 3 - high_leach_obs[soil]) ** 2

                    elif ii == 0 and soil == 2:
                        cum_time_6min_SA = cum_time_dt
                        cum_precip_6min_SA = cum_precip_dt
                        cum_runoff_6min_SA = cum_runoff_dt
                        cum_leach_6min_SA = cum_leach_dt
                        cum_inf_6min_SA = cum_inf_dt
                        runoff_6min_SA = runoff_dt
                        # print("cum runoff end: ", cum_runoff_dt[-1])
                        infil_6min_SA = infil_dt
                        # print("cum inf. end: ", cum_inf_dt[-1])
                        leach_6min_SA = leach_dt
                        # print("cum leach end: ", cum_leach_dt[-1])
                        time_size_6min_SA = time_size_dt
                        mass_bal_1_SA = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_high_SA = k
                        high_SA_error_prc = ((cum_leach_dt[5] / 10 ** 3 - high_leach_obs[soil]) / high_leach_obs[
                            soil]) * 100

                        # print("Sim:", cum_leach_dt[5] / 10 ** 3)
                        # print("Obs:", high_leach_obs[s])
                        # print(cum_leach_dt)
                        SS1_SA = (cum_leach_dt[5] / 10 ** 3 - high_leach_obs[soil]) ** 2

                    elif ii == 0 and soil == 3:
                        cum_time_6min_LA = cum_time_dt
                        cum_precip_6min_LA = cum_precip_dt
                        cum_runoff_6min_LA = cum_runoff_dt
                        cum_leach_6min_LA = cum_leach_dt
                        cum_inf_6min_LA = cum_inf_dt
                        runoff_6min_LA = runoff_dt
                        # print("cum runoff end: ", cum_runoff_dt[-1])
                        infil_6min_LA = infil_dt
                        # print("cum inf. end: ", cum_inf_dt[-1])
                        leach_6min_LA = leach_dt
                        # print("cum leach end: ", cum_leach_dt[-1])
                        time_size_6min_LA = time_size_dt
                        mass_bal_1_LA = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_high_LA = k
                        high_LA_error_prc = ((cum_leach_dt[5] / 10 ** 3 - high_leach_obs[soil]) / high_leach_obs[
                            soil]) * 100

                        # print("Sim:", cum_leach_dt[5] / 10 ** 3)
                        # print("Obs:", high_leach_obs[s])
                        # print(cum_leach_dt)
                        SS1_LA = (cum_leach_dt[5] / 10 ** 3 - high_leach_obs[soil]) ** 2

                    ####################################################
                    # Store hydro. variables  for intesity: 55 mm/h @ 12 min
                    elif ii == 1 and soil == 0:
                        cum_time_med_12min_SF = cum_time_dt
                        cum_precip_med_12min_SF = cum_precip_dt
                        cum_runoff_med_12min_SF = cum_runoff_dt
                        cum_leach_med_12min_SF = cum_leach_dt
                        cum_inf_med_12min_SF = cum_inf_dt
                        runoff_med_12min_SF = runoff_dt
                        infil_med_12min_SF = infil_dt
                        leach_med_12min_SF = leach_dt
                        time_size_med_12min_SF = time_size_dt
                        mass_bal_2_SF = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_med12_SF = k
                        med12_SF_error_prc = ((cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[soil]) / med12_leach_obs[soil]) * 100
                        SS2_SF = (cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[soil]) ** 2

                    elif ii == 1 and soil == 1:
                        cum_time_med_12min_LF = cum_time_dt
                        cum_precip_med_12min_LF = cum_precip_dt
                        cum_runoff_med_12min_LF = cum_runoff_dt
                        cum_leach_med_12min_LF = cum_leach_dt
                        cum_inf_med_12min_LF = cum_inf_dt
                        runoff_med_12min_LF = runoff_dt
                        infil_med_12min_LF = infil_dt
                        leach_med_12min_LF = leach_dt
                        time_size_med_12min_LF = time_size_dt
                        mass_bal_2_LF = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_med12_LF = k
                        med12_LF_error_prc = ((cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[soil]) / med12_leach_obs[soil]) * 100
                        SS2_LF = (cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[soil]) ** 2

                    elif ii == 1 and soil == 2:
                        cum_time_med_12min_SA = cum_time_dt
                        cum_precip_med_12min_SA = cum_precip_dt
                        cum_runoff_med_12min_SA = cum_runoff_dt
                        cum_leach_med_12min_SA = cum_leach_dt
                        cum_inf_med_12min_SA = cum_inf_dt
                        runoff_med_12min_SA = runoff_dt
                        infil_med_12min_SA = infil_dt
                        leach_med_12min_SA = leach_dt
                        time_size_med_12min_SA = time_size_dt
                        mass_bal_2_SA = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_med12_SA = k
                        med12_SA_error_prc = ((cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[soil]) / med12_leach_obs[soil]) * 100
                        SS2_SA = (cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[soil]) ** 2

                    elif ii == 1 and soil == 3:
                        cum_time_med_12min_LA = cum_time_dt
                        cum_precip_med_12min_LA = cum_precip_dt
                        cum_runoff_med_12min_LA = cum_runoff_dt
                        cum_leach_med_12min_LA = cum_leach_dt
                        cum_inf_med_12min_LA = cum_inf_dt
                        runoff_med_12min_LA = runoff_dt
                        infil_med_12min_LA = infil_dt
                        leach_med_12min_LA = leach_dt
                        time_size_med_12min_LA = time_size_dt
                        mass_bal_2_LA = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_med12_LA = k
                        med12_LA_error_prc = ((cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[soil]) / med12_leach_obs[soil]) * 100
                        SS2_LA = (cum_leach_dt[11] / 10 ** 3 - med12_leach_obs[soil]) ** 2

                    ####################################################
                    # Store hydro. variables  for intesity: 55 mm/h @ 30min
                    elif ii == 2 and soil == 0:
                        cum_time_med_30min_SF = cum_time_dt
                        cum_precip_med_30min_SF = cum_precip_dt
                        cum_runoff_med_30min_SF = cum_runoff_dt
                        cum_leach_med_30min_SF = cum_leach_dt
                        cum_inf_med_30min_SF = cum_inf_dt
                        runoff_med_30min_SF = runoff_dt
                        infil_med_30min_SF = infil_dt
                        leach_med_30min_SF = leach_dt
                        time_size_med_30min_SF = time_size_dt
                        mass_bal_3_SF = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_med30_SF = k

                        med30_SF_error_prc = ((cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[soil]) / med30_leach_obs[soil]) * 100
                        SS3_SF = (cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[soil]) ** 2

                    elif ii == 2 and soil == 1:
                        cum_time_med_30min_LF = cum_time_dt
                        cum_precip_med_30min_LF = cum_precip_dt
                        cum_runoff_med_30min_LF = cum_runoff_dt
                        cum_leach_med_30min_LF = cum_leach_dt
                        cum_inf_med_30min_LF = cum_inf_dt
                        runoff_med_30min_LF = runoff_dt
                        infil_med_30min_LF = infil_dt
                        leach_med_30min_LF = leach_dt
                        time_size_med_30min_LF = time_size_dt
                        mass_bal_3_LF = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_med30_LF = k

                        med30_LF_error_prc = ((cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[soil]) / med30_leach_obs[soil]) * 100
                        SS3_LF = (cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[soil]) ** 2

                    elif ii == 2 and soil == 2:
                        cum_time_med_30min_SA = cum_time_dt
                        cum_precip_med_30min_SA = cum_precip_dt
                        cum_runoff_med_30min_SA = cum_runoff_dt
                        cum_leach_med_30min_SA = cum_leach_dt
                        cum_inf_med_30min_SA = cum_inf_dt
                        runoff_med_30min_SA = runoff_dt
                        infil_med_30min_SA = infil_dt
                        leach_med_30min_SA = leach_dt
                        time_size_med_30min_SA = time_size_dt
                        mass_bal_3_SA = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_med30_SA = k

                        med30_SA_error_prc = ((cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[soil]) / med30_leach_obs[soil]) * 100
                        SS3_SA = (cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[soil]) ** 2

                    elif ii == 2 and soil == 3:
                        cum_time_med_30min_LA = cum_time_dt
                        cum_precip_med_30min_LA = cum_precip_dt
                        cum_runoff_med_30min_LA = cum_runoff_dt
                        cum_leach_med_30min_LA = cum_leach_dt
                        cum_inf_med_30min_LA = cum_inf_dt
                        runoff_med_30min_LA = runoff_dt
                        infil_med_30min_LA = infil_dt
                        leach_med_30min_LA = leach_dt
                        time_size_med_30min_LA = time_size_dt
                        mass_bal_3_LA = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_med30_LA = k

                        med30_LA_error_prc = ((cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[soil]) / med30_leach_obs[soil]) * 100
                        SS3_LA = (cum_leach_dt[-1] / 10 ** 3 - med30_leach_obs[soil]) ** 2

                    ####################################################
                    # Store hydro. variables  for intesity: 30 mm/h
                    elif ii == 3 and soil == 0:
                        cum_time_30min_SF = cum_time_dt
                        cum_precip_low_30min_SF = cum_precip_dt
                        cum_runoff_low_30min_SF = cum_runoff_dt
                        cum_leach_low_30min_SF = cum_leach_dt
                        cum_inf_low_30min_SF = cum_inf_dt
                        runoff_low_30min_SF = runoff_dt
                        infil_low_30min_SF = infil_dt
                        leach_low_30min_SF = leach_dt
                        time_size_low_30min_SF = time_size_dt
                        mass_bal_4_SF = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_low_SF = k
                        low_SF_error_prc = ((cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[soil]) / low_leach_obs[soil]) * 100
                        SS4_SF = (cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[soil]) ** 2

                    elif ii == 3 and soil == 1:
                        cum_time_30min_LF = cum_time_dt
                        cum_precip_low_30min_LF = cum_precip_dt
                        cum_runoff_low_30min_LF = cum_runoff_dt
                        cum_leach_low_30min_LF = cum_leach_dt
                        cum_inf_low_30min_LF = cum_inf_dt
                        runoff_low_30min_LF = runoff_dt
                        infil_low_30min_LF = infil_dt
                        leach_low_30min_LF = leach_dt
                        time_size_low_30min_LF = time_size_dt
                        mass_bal_4_LF = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_low_LF = k
                        low_LF_error_prc = ((cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[soil]) / low_leach_obs[soil]) * 100
                        SS4_LF = (cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[soil]) ** 2

                    elif ii == 3 and soil == 2:
                        cum_time_30min_SA = cum_time_dt
                        cum_precip_low_30min_SA = cum_precip_dt
                        cum_runoff_low_30min_SA = cum_runoff_dt
                        cum_leach_low_30min_SA = cum_leach_dt
                        cum_inf_low_30min_SA = cum_inf_dt
                        runoff_low_30min_SA = runoff_dt
                        infil_low_30min_SA = infil_dt
                        leach_low_30min_SA = leach_dt
                        time_size_low_30min_SA = time_size_dt
                        mass_bal_4_SA = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_low_SA = k
                        low_SA_error_prc = ((cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[soil]) / low_leach_obs[soil]) * 100
                        SS4_SA = (cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[soil]) ** 2

                    elif ii == 3 and soil == 3:
                        cum_time_30min_LA = cum_time_dt
                        cum_precip_low_30min_LA = cum_precip_dt
                        cum_runoff_low_30min_LA = cum_runoff_dt
                        cum_leach_low_30min_LA = cum_leach_dt
                        cum_inf_low_30min_LA = cum_inf_dt
                        runoff_low_30min_LA = runoff_dt
                        infil_low_30min_LA = infil_dt
                        leach_low_30min_LA = leach_dt
                        time_size_low_30min_LA = time_size_dt
                        mass_bal_4_LA = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
                        k_low_LA = k
                        low_LA_error_prc = ((cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[soil]) / low_leach_obs[soil]) * 100
                        SS4_LA = (cum_leach_dt[-1] / 10 ** 3 - low_leach_obs[soil]) ** 2

                    else:
                        print("Do you have more than 3 rainfall scenarios?")

                else:
                    continue

    mean_obs = (np.mean(high_leach_obs) + np.mean(med12_leach_obs) + np.mean(med30_leach_obs) + np.mean(low_leach_obs))/4
    SStot = 0
    for obs in high_leach_obs:
        SStot += (obs-mean_obs)**2
    for obs in med12_leach_obs:
        SStot += (obs-mean_obs)**2
    for obs in med30_leach_obs:
        SStot += (obs-mean_obs)**2
    for obs in low_leach_obs:
        SStot += (obs-mean_obs)**2

    SSres = SS1_SF + SS2_SF + SS3_SF + SS4_SF + \
            SS1_LF + SS2_LF + SS3_LF + SS4_LF + \
            SS1_SA + SS2_SA + SS3_SA + SS4_SA + \
            SS1_LA + SS2_LA + SS3_LA + SS4_LA

    r_squared = 1 - (SSres / SStot)
    print("R2: ", r_squared)
    print("--------------------------------------------")
    print("Sterile, Fresh Soil")
    print("ksat high: ", k_high_SF/10*60, "cm/h | Error %", high_SF_error_prc)
    print("ksat med12: ", k_med12_SF/10*60, "cm/h | Error %", med12_SF_error_prc)
    print("ksat med30: ", k_med30_SF/10*60, "cm/h | Error %", med30_SF_error_prc)
    print("ksat low: ", k_low_SF/10*60, "cm/h | Error %", low_SF_error_prc)
    print("--------------------------------------------")
    print("Sterile, Aged Soil")
    print("ksat high: ", k_high_SA/10*60, "cm/h | Error %", high_SA_error_prc)
    print("ksat med12: ", k_med12_SA/10*60, "cm/h | Error %", med12_SA_error_prc)
    print("ksat med30: ", k_med30_SA/10*60, "cm/h | Error %", med30_SA_error_prc)
    print("ksat low: ", k_low_SA/10*60, "cm/h | Error %", low_SA_error_prc)
    print("--------------------------------------------")
    print("Living, Fresh Soil")
    print("ksat high: ", k_high_LF / 10 * 60, "cm/h | Error %", high_LF_error_prc)
    print("ksat med12: ", k_med12_LF / 10 * 60, "cm/h | Error %", med12_LF_error_prc)
    print("ksat med30: ", k_med30_LF / 10 * 60, "cm/h | Error %", med30_LF_error_prc)
    print("ksat low: ", k_low_LF / 10 * 60, "cm/h | Error %", low_LF_error_prc)
    print("--------------------------------------------")
    print("Living, Fresh Aged")
    print("ksat high: ", k_high_LA / 10 * 60, "cm/h | Error %", high_LA_error_prc)
    print("ksat med12: ", k_med12_LA / 10 * 60, "cm/h | Error %", med12_LA_error_prc)
    print("ksat med30: ", k_med30_LA / 10 * 60, "cm/h | Error %", med30_LA_error_prc)
    print("ksat low: ", k_low_LA / 10 * 60, "cm/h | Error %", low_LA_error_prc)
    print("--------------------------------------------")

    # print("Mass balance: ", abs(mass_bal_1) < 10 ** (-6), abs(mass_bal_2) < 10 ** (-6), abs(mass_bal_3) < 10 ** (-6),
    #      abs(mass_bal_4) < 10 ** (-6))

    # Data in mm3
    # Available for each modality with eg. "_SF" subscript but not returned are also:
    # leach_6min, leach_med_12min_SF, leach_med_30min, leach_low_30min,
    # cum_inf_6min, cum_inf_med_12min_SF, cum_inf_med_30min, cum_inf_low_30min,
    # infil_6min, infil_med_12min_SF, infil_med_30min, infil_low_30min,
    return stackdata36(
        cum_time_30min_SF,  # 0
        cum_runoff_6min_SF, cum_leach_6min_SF,  # 2
        cum_runoff_6min_SA, cum_leach_6min_SA,  # 4
        cum_runoff_6min_LF, cum_leach_6min_LF,
        cum_runoff_6min_LA, cum_leach_6min_LA,  # 8
        cum_runoff_med_12min_SF, cum_leach_med_12min_SF,
        cum_runoff_med_12min_SA, cum_leach_med_12min_SA,
        cum_runoff_med_12min_LF, cum_leach_med_12min_LF,
        cum_runoff_med_12min_LA, cum_leach_med_12min_LA,  # 16
        cum_runoff_med_30min_SF, cum_leach_med_30min_SF,
        cum_runoff_med_30min_SA, cum_leach_med_30min_SA,
        cum_runoff_med_30min_LF, cum_leach_med_30min_LF,
        cum_runoff_med_30min_LA, cum_leach_med_30min_LA,
        cum_runoff_low_30min_SF, cum_leach_low_30min_SF,
        cum_runoff_low_30min_SA, cum_leach_low_30min_SA,  # 28
        cum_runoff_low_30min_LF, cum_leach_low_30min_LF,
        cum_runoff_low_30min_LA, cum_leach_low_30min_LA,
        time_size_6min_SF, time_size_med_12min_SF, time_size_med_30min_SF, time_size_low_30min_SF)  # 36
