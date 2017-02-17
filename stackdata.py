import numpy as np


def stackdata36(time,
                a, b, c, d, e, f, g, h, i, j,
                k, l, m, n, o, p, q, r, s, t,
                u, v, w, x, y, z, aa, ab, ac, ad,
                af, ag, ah, ai, aj, ak):
    data = np.column_stack((
        time,
        a, b, c, d, e, f, g, h, i, j, k, l, m, n,
        o, p, q, r, s, t, u, v, w, x, y, z, aa, ab, ac, ad,
        af, ag, ah, ai, aj, ak))
    return data


def stackdata30(time,
                a, b, c, d,
                e, f, g, h,
                i, j, k, l,
                m, n, o, p,
                q, r, s, t,
                u, v, w, x,
                y, z, aa, ab,
                ac, ad):
    data = np.column_stack((
        time,
        a, b, c, d, e, f, g, h, i, j, k, l, m, n,
        o, p, q, r, s, t, u, v, w, x, y, z, aa, ab, ac, ad ))
    return data


def stackdata28(time,
                a, b, c, d,
                e, f, g, h,
                i, j, k, l,
                m, n, o, p,
                q, r, s, t,
                u, v, w, x,
                y, z, aa, ab):
    data = np.column_stack((time,
                            a, b, c, d, e, f, g, h, i, j, k, l, m, n,
                            o, p, q, r, s, t, u, v, w, x, y, z, aa, ab))
    return data

def stackdata21(time,
                a, b, c, d,
                e, f, g, h,
                i, j, k, l,
                m, n, o,
                p, q, r,
                s, t, u):
    data = np.column_stack((time, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u))
    return data


def stackdata18(time,
                a, b, c, d,
                e, f, g, h,
                i, j, k, l,
                m, n, o,
                p, q, r):
    data = np.column_stack((time, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r))
    return data


def stackdata16(time,
                a, b, c, d,
                e, f, g, h,
                i, j, k, l,
                m, n, o, p):
    data = np.column_stack((time, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p))
    return data


def stackdata15(time,
                a, b, c, d,
                e, f, g, h,
                i, j, k, l,
                m, n, o):
    data = np.column_stack((time, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o))
    return data


def stackdata9(time, l1, l2, l3, l4, l5, l6, l7, l8, l9):
    data = np.column_stack((time, l1, l2, l3, l4, l5, l6, l7, l8, l9))
    return data


def stackdata8(time, l1, l2, l3, l4, l5, l6, l7, l8):
    data = np.column_stack((time, l1, l2, l3, l4, l5, l6, l7, l8))
    return data


def stackdata6(time, l1, l2, l3, l4, l5, l6):
    data = np.column_stack((time, l1, l2, l3, l4, l5, l6))
    return data


def stackdata4(time, a, b, c, d):
    data = np.column_stack((time, a, b, c, d))
    return data


def stackdata3(time, a, b, c):
    data = np.column_stack((time, a, b, c))
    return data