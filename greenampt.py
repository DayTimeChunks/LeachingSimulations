from math import log


def greenampt(ksat, wfs, dtga, precip, ovsat, ov, zl, cum_infilt):

    """
    Function returns a list: [F2, inf, runoff, leach].
    Where,
        F2 = cumulative infiltration;
        inf = F2 - F1 (infiltration at time t)
        runoff = precip - inf
        leach = amount leached out of the unit volume at time t

    Calculation of Green-Ampt Infiltration in X steps (i.e. if->while sequences).
    1) Define infiltration capacity (f)
    2) Iterate to calculate cumulative infiltration during time step size (F2)
        2.1) Test if ponding time "tp" occurs in the middle of the time step (tp)

    
    wfs = suction pressure * change in moisture
    """

    """ 1) """
    # Infiltration Capacity "f"
    F1 = cum_infilt

    if F1 > 0.0:  # replaced (5.0 * 1e-3) with 0.
        f = ksat * (1 + wfs / F1)
    else:
        f = 100.0 * ksat

    # print(f)
    # print(ksat)

    """ 2) Cumulative infiltration - Green-Ampt infiltration calculation """
    if precip > f:
        Fold2 = 0
        error = 1
        cc = 1

        while (error > 0.00001) and (cc < 1000):
            Fnew = F1 + (ksat * dtga) + wfs * log((Fold2 + wfs) / (F1 + wfs))
            error = abs((Fnew - Fold2) / Fnew)
            Fold2 = Fnew
            cc += 1
            F2 = Fnew
            runoff = precip * dtga - (F2 - F1)  # runoff (cm) = precip - infil

            if cc >= 1000:
                print("GAInflt recursion error 1")
                break

            # Test for large error in iteration
            if F2 < F1:
                print("Error F2 < F1, you have to decrease the error in G&A iteration")
                break

                # If no error, calculation of Fnew has converged.
    else:
        F2 = F1 + precip

        """ 2.1) Check if infiltration cap. has been exceeded within time step - TestCondition 1 """
        ff2 = ksat * (1 + wfs / F2)  # ff2 = new infiltration capacity

        if precip < ff2:
            runoff = 0
        else:
            # "potential ponding during time step"
            Fp = ksat * wfs / (precip - ksat)
            tp = (Fp - F1) / precip  # calculates time at ponding

            if tp > dtga:
                runoff = 0
            else:
                # Calculate F2 iteratively again, but with F1 being updated until ponding occurs
                print("G&A ponding during time step")
                # TestCondition = 3
                F1 = F1 + precip * tp
                Fold2 = 0.0
                error = 1
                cc = 1

                while (error > 0.00001) and (cc < 1000):
                    Fnew = F1 + ksat * (dtga - tp) + wfs * log((Fold2 + wfs) / (F1 + wfs))
                    error = (Fnew - Fold2) / Fnew
                    Fold2 = Fnew
                    cc += 1
                    # Partial changes: F2 & runoff:

                    F2 = Fnew  # Cumul. infiltration (mm)
                    runoff = precip * (dtga - tp) - (F2 - F1)  # Runoff (mm)

                    if cc >= 1000:
                        print("GAInflt recursion error 2")
                        break
                    elif F2 < F1:
                        print("Error F2 < F1, you have to decrease the error in G&A iteration")
                        break

                # print("runoff from intermediate dt =", runoff)

    waterfront = F2 / (ovsat - ov)

    inf = F2 - F1

    if waterfront > zl:
        leach = (F2 - F1)
    else:
        leach = 0

    cumInfilt = F2

    return [F2, inf, runoff, leach]
