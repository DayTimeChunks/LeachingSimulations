from greenampt import *
from hydroplots import *


def leachsim(dtGA = 1,  # Timestep in minute
             StormD = 30,  # Storm duration in min
             intensityM = [2.25, 0.92, 0.5],  # mm/min = [0.225, 0.092, 0.05] cm/min
             numSystems = 3,  # Number of combined replicas
             ov=0.25,   # Initial water content m3. m-3
             ovSat = 0.45,  # Saturated water content (0.61 - Alteck; 0.55 - Rouff)
             kSat = 2.25,  # mm/min (13.5 cm/h - Crop Soil) - Alteck
             psi=1100,  # soil suction Alteck mm
             soil_height=20  # mm
             ):
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
            cum_inf_6min = cum_inf_dt
            runoff_6min = runoff_dt
            # print("cum runoff end: ", cum_runoff_dt[-1])
            infil_6min = infil_dt
            # print("cum inf. end: ", cum_inf_dt[-1])
            leach_6min = leach_dt
            # print("cum leach end: ", cum_leach_dt[-1])
            time_size_6min = time_size_dt
            mass_bal = cum_precip_dt[-1] - (cum_runoff_dt[-1] + cum_inf_dt[-1])
            print("Mass balance", mass_bal)

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
        dtGA = 1,  # Timestep in minute
        StormD = 30,  # Storm duration in min
        intensityM = [2.25, 0.92, 0.92, 0.5],  # mm/min = [0.225, 0.092, 0.05] cm/min
        numSystems = 3,  # Number of combined replicas
        psi=1100,  # soil suction Alteck mm
        AGED = False):

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

            if not AGED:
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

                else:
                    print("Do you have more than 3 rainfall scenarios?")

            else:
                continue
    print(mass_bal_1, mass_bal_2, mass_bal_3, mass_bal_4)
    print("ksat high: ", k_high)
    print("ksat med12: ", k_med12)
    print("ksat med30: ", k_med30)
    print("ksat low: ", k_low)

    # Data in mm3
    return stackdata28(cum_time_30min,
                       leach_6min, leach_med_12min, leach_med_30min, leach_low_30min,
                       cum_inf_6min, cum_inf_med_12min, cum_inf_med_30min, cum_inf_low_30min,
                       cum_leach_6min, cum_leach_med_12min, cum_leach_med_30min, cum_leach_low_30min,
                       runoff_6min, runoff_med_12min, runoff_med_30min, runoff_low_30min,
                       cum_runoff_6min, cum_runoff_med_12min, cum_runoff_med_30min, cum_runoff_low_30min,
                       infil_6min, infil_med_12min, infil_med_30min, infil_low_30min,
                       time_size_6min, time_size_med_12min, time_size_med_30min, time_size_low_30min)