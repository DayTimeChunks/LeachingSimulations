from greenampt import *
from hydroplots import *


def leachsim(dtGA = 0.25,  # Timestep in minute fraction
             StormD = 30,  # Storm duration in minutes
             intensityM = [0.225, 0.092, 0.05],  # cm/min = [135, 55, 30] mm/h
             numSystems = 3,  # Number of combined replicas
             ov=0.25,   # Initial water content m3. m-3
             ovSat = 0.45,  # Saturated water content (assumed)
             kSat = 0.225,  # cm/min (13.5 cm/h - Crop Soil) - Alteck
             psi=110,  # soil suction Alteck
             soil_height=3.0  # cm
             ):
    wfs = psi * (ovSat - ov)  # Numerator in Green-Ampt's equation = (suction x distance to saturation)
    """ Microcosm Dimensions """
    zl = soil_height  # Mixing layer depth in cm
    d = (1.493 * 2)  # Diameter of falcon tube (cm)
    area = ((d / 2) ** 2) * 3.1416  # (cm2)



    loop = 0
    # Running each of the intesities in one loop
    for i in range(0, len(intensityM)):
        loop += 1
        """ Reynolds parameters"""
        v_ro = intensityM[i] - kSat  # Runoff velocity (cm/min)

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

        """ Start of rainfall event X """
        while cum_time < StormD:
            cum_time += dtGA
            p = intensityM[i] * dtGA  # Precipitation (cm) (in one time step)
            cum_precip += p

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
                                     precip=p,
                                     ovsat=ovSat,
                                     ov=ov,
                                     zl=zl,
                                     cum_infilt=cum_inf)

            cum_inf = green_output[0]
            inf = green_output[1]
            runoff = green_output[2]
            leach = green_output[3]

            leached_vol = leach * area * numSystems

            infil_dt.append(inf * area * numSystems)
            runoff_dt.append(runoff * area * numSystems)
            leach_dt.append(leached_vol)

            cum_time_dt.append(cum_time)

            cum_precip_dt.append(cum_precip * area * numSystems)
            cum_inf_dt.append(cum_inf * area * numSystems)

            cum_runoff += runoff
            cum_runoff_dt.append(cum_runoff * area * numSystems)

            cum_leach += leach
            cum_leach_dt.append(cum_leach * area * numSystems)

        # Store hydro. variables  for intesity: 135 mm/h
        if loop == 1:
            cum_time_6min = cum_time_dt
            cum_precip_6min = cum_precip_dt
            cum_runoff_6min = cum_runoff_dt
            cum_leach_6min = cum_leach_dt
            cum_inf_6min = cum_inf_dt
            runoff_6min = runoff_dt
            infil_6min = infil_dt
            leach_6min = leach_dt

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

        else:
            print("Do you have more than 3 rainfall scenarios?")

    return stackdata9(cum_time_30min,
                     leach_6min, leach_12min, leach_30min,
                     cum_inf_6min, cum_inf_12min, cum_inf_30min,
                     cum_leach_6min, cum_leach_12min, cum_leach_30min)
