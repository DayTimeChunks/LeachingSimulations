from math import log


def greenampt(ksat, wfs, dtga, intensity, ovsat, ov, zl, cum_infilt, cum_time):

    """
    Function returns a list: [F2, inf, ponding, leach].
    Where,
        F2 = cumulative infiltration;
        inf = F2 - F1 (infiltration at time t)
        ponding = precip - inf
        leach = amount leached out of the unit volume at time t

    Calculation of Green-Ampt Infiltration in X steps (i.e. if->while sequences).
    1) Define infiltration capacity (f)
    2) Iterate to calculate cumulative infiltration during time step size (F2)
        2.1) Test if ponding time "tp" occurs in the middle of the time step (tp)

    
    wfs = suction pressure * change in moisture
    """

    """ 1) """

    F1 = cum_infilt
    # print("F1 - in:", F1)

    # Infiltration Capacity "f"
    if F1 > 0:  # could replace (5.0 * 1e-3) with 0.
        f = ksat * (1 + (wfs / F1))
    else:
        f = 1000.0 * ksat  # infinite

    # print(f)
    # print(ksat)

    """ 2.1) Cumulative infiltration - Green-Ampt infiltration calculation """
    if intensity > f:
        # print("G&A ponding ongoing with cum infil (F1): ", F1, "at Cum time: ", cum_time)

        tp = 0
        Fold2 = 0
        error = 1
        cc = 1

        while (error > 0.000001) and (cc < 1000):
            Fnew = F1 + (ksat * dtga) + wfs * log( (Fold2 + wfs) / (F1 + wfs) )
            error = abs((Fnew - Fold2) / Fnew)
            Fold2 = Fnew
            cc += 1

            if cc >= 1000:
                print("GAInflt recursion error 1")
                break

        F2 = Fnew

        # Test for large error in iteration
        if F2 < F1:
            print("Error 2.1, F2 < F1, you have to decrease the error in G&A iteration")
            print("F2: ", F2, "F1: ", F1, "cc: ", cc, "error:", error)

        # If no error, calculation of Fnew has converged.
        ponding = intensity * dtga - (F2 - F1)  # ponding (mm) = precip - infil
        # print("ponding in 2.1", ponding, "F2 at ponding: ", F2)

    else:
        F2 = F1 + intensity * dtga
        ponding = 0

        """ 2.2) Check if infiltration cap. has been exceeded within time step - TestCondition 1 """
        ff2 = ksat * (1 + wfs / F2)  # ff2 = new infiltration capacity

        if intensity < ff2:
            ponding = 0
            # print("F2 at false inter-ponding 1: ", F2)
        else:
            # "potential ponding during time step"
            Fp = ksat * wfs / (intensity - ksat)  # [mm]
            tp = (Fp - F1) / intensity  # calculates time at ponding [min]

            if tp > dtga:
                ponding = 0
                # print("F2 at false inter-ponding 2: ", F2)
            else:
                # Calculate F2 iteratively again, but with F1 being updated until ponding occurs
                print("G&A ponding during time step, with size:", tp)
                # TestCondition = 3
                F1 = F1 + intensity * tp
                Fold2 = 0.0
                error = 1
                cc = 1

                while (error > 0.00001) and (cc < 1000):
                    Fnew = F1 + ksat * (dtga - tp) + wfs * log((Fold2 + wfs) / (F1 + wfs))
                    error = (Fnew - Fold2) / Fnew
                    Fold2 = Fnew
                    cc += 1
                    # Partial changes: F2 & ponding:

                    if cc >= 1000:
                        print("GAInflt recursion error 2")
                        break

                F2 = Fnew  # Cumul. infiltration (mm)
                # print("F2 after TRUE inter-ponding: ", F2)

                if F2 < F1:
                    print("Error 2.2, F2 < F1, you have to decrease the error in G&A iteration")
                    print("Size of error: ", F2 - F1)

                ponding = intensity * (dtga - tp) - (F2 - F1)  # Runoff (mm)
                # print("ponding depth at intermediate timesetp (tp) = ", ponding)

    waterfront = F2 / (ovsat - ov)

    inf = F2 - cum_infilt

    if waterfront > zl:
        leach = (F2 - cum_infilt)
    else:
        leach = 0

    # print("Final F2", F2, "inf(t): ", inf, "leach(t): ", leach)
    # print("intes*dtGA", intensity * dtga)

    return [F2, inf, ponding, leach]
