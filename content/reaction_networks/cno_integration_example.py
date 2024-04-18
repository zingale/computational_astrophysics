import numba
import numpy as np
from scipy import constants
from numba.experimental import jitclass

from pynucastro.rates import TableIndex, TableInterpolator, TabularRate, Tfactors
from pynucastro.screening import PlasmaState, ScreenFactors

jp = 0
jhe4 = 1
jc12 = 2
jc13 = 3
jn13 = 4
jn14 = 5
jn15 = 6
jo14 = 7
jo15 = 8
jo16 = 9
jo17 = 10
jo18 = 11
jf17 = 12
jf18 = 13
jf19 = 14
jne18 = 15
nnuc = 16

A = np.zeros((nnuc), dtype=np.int32)

A[jp] = 1
A[jhe4] = 4
A[jc12] = 12
A[jc13] = 13
A[jn13] = 13
A[jn14] = 14
A[jn15] = 15
A[jo14] = 14
A[jo15] = 15
A[jo16] = 16
A[jo17] = 17
A[jo18] = 18
A[jf17] = 17
A[jf18] = 18
A[jf19] = 19
A[jne18] = 18

Z = np.zeros((nnuc), dtype=np.int32)

Z[jp] = 1
Z[jhe4] = 2
Z[jc12] = 6
Z[jc13] = 6
Z[jn13] = 7
Z[jn14] = 7
Z[jn15] = 7
Z[jo14] = 8
Z[jo15] = 8
Z[jo16] = 8
Z[jo17] = 8
Z[jo18] = 8
Z[jf17] = 9
Z[jf18] = 9
Z[jf19] = 9
Z[jne18] = 10

# masses in ergs
mass = np.zeros((nnuc), dtype=np.float64)

mass[jp] = 0.0015040963030260536
mass[jhe4] = 0.0059735574925878256
mass[jc12] = 0.017909017027273523
mass[jc13] = 0.019406441930882663
mass[jn13] = 0.01940999951603316
mass[jn14] = 0.020898440903103103
mass[jn15] = 0.0223864338056853
mass[jo14] = 0.020906683076491985
mass[jo15] = 0.02239084645968795
mass[jo16] = 0.023871099858982767
mass[jo17] = 0.02536981167252093
mass[jo18] = 0.02686227133140636
mass[jf17] = 0.025374234423440733
mass[jf18] = 0.026864924401329426
mass[jf19] = 0.028353560468882166
mass[jne18] = 0.026872045275379234

names = []
names.append("H1")
names.append("He4")
names.append("C12")
names.append("C13")
names.append("N13")
names.append("N14")
names.append("N15")
names.append("O14")
names.append("O15")
names.append("O16")
names.append("O17")
names.append("O18")
names.append("F17")
names.append("F18")
names.append("F19")
names.append("Ne18")

def to_composition(Y):
    """Convert an array of molar fractions to a Composition object."""
    from pynucastro import Composition, Nucleus
    nuclei = [Nucleus.from_cache(name) for name in names]
    comp = Composition(nuclei)
    for i, nuc in enumerate(nuclei):
        comp.X[nuc] = Y[i] * A[i]
    return comp


def energy_release(dY):
    """return the energy release in erg/g (/s if dY is actually dY/dt)"""
    enuc = 0.0
    for i, y in enumerate(dY):
        enuc += y * mass[i]
    enuc *= -1*constants.Avogadro
    return enuc

@jitclass([
    ("p_C12__N13", numba.float64),
    ("p_C13__N14", numba.float64),
    ("N13__C13__weak__wc12", numba.float64),
    ("p_N13__O14", numba.float64),
    ("p_N14__O15", numba.float64),
    ("p_N15__He4_C12", numba.float64),
    ("O14__N14__weak__wc12", numba.float64),
    ("O15__N15__weak__wc12", numba.float64),
    ("p_N15__O16", numba.float64),
    ("p_O16__F17", numba.float64),
    ("F17__O17__weak__wc12", numba.float64),
    ("p_O17__He4_N14", numba.float64),
    ("p_O17__F18", numba.float64),
    ("F18__O18__weak__wc12", numba.float64),
    ("p_O18__He4_N15", numba.float64),
    ("p_O18__F19", numba.float64),
    ("p_F19__He4_O16", numba.float64),
    ("He4_O14__p_F17", numba.float64),
    ("p_F17__Ne18", numba.float64),
    ("Ne18__F18__weak__wc12", numba.float64),
    ("p_F18__He4_O15", numba.float64),
])
class RateEval:
    def __init__(self):
        self.p_C12__N13 = np.nan
        self.p_C13__N14 = np.nan
        self.N13__C13__weak__wc12 = np.nan
        self.p_N13__O14 = np.nan
        self.p_N14__O15 = np.nan
        self.p_N15__He4_C12 = np.nan
        self.O14__N14__weak__wc12 = np.nan
        self.O15__N15__weak__wc12 = np.nan
        self.p_N15__O16 = np.nan
        self.p_O16__F17 = np.nan
        self.F17__O17__weak__wc12 = np.nan
        self.p_O17__He4_N14 = np.nan
        self.p_O17__F18 = np.nan
        self.F18__O18__weak__wc12 = np.nan
        self.p_O18__He4_N15 = np.nan
        self.p_O18__F19 = np.nan
        self.p_F19__He4_O16 = np.nan
        self.He4_O14__p_F17 = np.nan
        self.p_F17__Ne18 = np.nan
        self.Ne18__F18__weak__wc12 = np.nan
        self.p_F18__He4_O15 = np.nan

@numba.njit()
def ye(Y):
    return np.sum(Z * Y)/np.sum(A * Y)

@numba.njit()
def p_C12__N13(rate_eval, tf):
    # C12 + p --> N13
    rate = 0.0

    # ls09n
    rate += np.exp(  17.1482 + -13.692*tf.T913i + -0.230881*tf.T913
                  + 4.44362*tf.T9 + -3.15898*tf.T953 + -0.666667*tf.lnT9)
    # ls09r
    rate += np.exp(  17.5428 + -3.77849*tf.T9i + -5.10735*tf.T913i + -2.24111*tf.T913
                  + 0.148883*tf.T9 + -1.5*tf.lnT9)

    rate_eval.p_C12__N13 = rate

@numba.njit()
def p_C13__N14(rate_eval, tf):
    # C13 + p --> N14
    rate = 0.0

    # nacrr
    rate += np.exp(  15.1825 + -13.5543*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  18.5155 + -13.72*tf.T913i + -0.450018*tf.T913
                  + 3.70823*tf.T9 + -1.70545*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  13.9637 + -5.78147*tf.T9i + -0.196703*tf.T913
                  + 0.142126*tf.T9 + -0.0238912*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_C13__N14 = rate

@numba.njit()
def N13__C13__weak__wc12(rate_eval, tf):
    # N13 --> C13
    rate = 0.0

    # wc12w
    rate += np.exp(  -6.7601)

    rate_eval.N13__C13__weak__wc12 = rate

@numba.njit()
def p_N13__O14(rate_eval, tf):
    # N13 + p --> O14
    rate = 0.0

    # lg06r
    rate += np.exp(  10.9971 + -6.12602*tf.T9i + 1.57122*tf.T913i
                  + -1.5*tf.lnT9)
    # lg06n
    rate += np.exp(  18.1356 + -15.1676*tf.T913i + 0.0955166*tf.T913
                  + 3.0659*tf.T9 + -0.507339*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_N13__O14 = rate

@numba.njit()
def p_N14__O15(rate_eval, tf):
    # N14 + p --> O15
    rate = 0.0

    # im05n
    rate += np.exp(  17.01 + -15.193*tf.T913i + -0.161954*tf.T913
                  + -7.52123*tf.T9 + -0.987565*tf.T953 + -0.666667*tf.lnT9)
    # im05r
    rate += np.exp(  6.73578 + -4.891*tf.T9i
                  + 0.0682*tf.lnT9)
    # im05r
    rate += np.exp(  7.65444 + -2.998*tf.T9i
                  + -1.5*tf.lnT9)
    # im05n
    rate += np.exp(  20.1169 + -15.193*tf.T913i + -4.63975*tf.T913
                  + 9.73458*tf.T9 + -9.55051*tf.T953 + 0.333333*tf.lnT9)

    rate_eval.p_N14__O15 = rate

@numba.njit()
def p_N15__He4_C12(rate_eval, tf):
    # N15 + p --> He4 + C12
    rate = 0.0

    # nacrn
    rate += np.exp(  27.4764 + -15.253*tf.T913i + 1.59318*tf.T913
                  + 2.4479*tf.T9 + -2.19708*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  -6.57522 + -1.1638*tf.T9i + 22.7105*tf.T913
                  + -2.90707*tf.T9 + 0.205754*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  20.8972 + -7.406*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  -4.87347 + -2.02117*tf.T9i + 30.8497*tf.T913
                  + -8.50433*tf.T9 + -1.54426*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_N15__He4_C12 = rate

@numba.njit()
def O14__N14__weak__wc12(rate_eval, tf):
    # O14 --> N14
    rate = 0.0

    # wc12w
    rate += np.exp(  -4.62354)

    rate_eval.O14__N14__weak__wc12 = rate

@numba.njit()
def O15__N15__weak__wc12(rate_eval, tf):
    # O15 --> N15
    rate = 0.0

    # wc12w
    rate += np.exp(  -5.17053)

    rate_eval.O15__N15__weak__wc12 = rate

@numba.njit()
def p_N15__O16(rate_eval, tf):
    # N15 + p --> O16
    rate = 0.0

    # li10n
    rate += np.exp(  20.0176 + -15.24*tf.T913i + 0.334926*tf.T913
                  + 4.59088*tf.T9 + -4.78468*tf.T953 + -0.666667*tf.lnT9)
    # li10r
    rate += np.exp(  14.5444 + -10.2295*tf.T9i
                  + 0.0459037*tf.T9 + -1.5*tf.lnT9)
    # li10r
    rate += np.exp(  6.59056 + -2.92315*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_N15__O16 = rate

@numba.njit()
def p_O16__F17(rate_eval, tf):
    # O16 + p --> F17
    rate = 0.0

    # ia08n
    rate += np.exp(  19.0904 + -16.696*tf.T913i + -1.16252*tf.T913
                  + 0.267703*tf.T9 + -0.0338411*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_O16__F17 = rate

@numba.njit()
def F17__O17__weak__wc12(rate_eval, tf):
    # F17 --> O17
    rate = 0.0

    # wc12w
    rate += np.exp(  -4.53318)

    rate_eval.F17__O17__weak__wc12 = rate

@numba.njit()
def p_O17__He4_N14(rate_eval, tf):
    # O17 + p --> He4 + N14
    rate = 0.0

    # il10r
    rate += np.exp(  5.5336 + -2.11477*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -7.20763 + -0.753395*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.579 + -16.9078*tf.T913i
                  + -2.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  10.174 + -4.95865*tf.T9i + 5.10182*tf.T913
                  + 0.379373*tf.T9 + -0.0672515*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_O17__He4_N14 = rate

@numba.njit()
def p_O17__F18(rate_eval, tf):
    # O17 + p --> F18
    rate = 0.0

    # il10n
    rate += np.exp(  15.8929 + -16.4035*tf.T913i + 4.31885*tf.T913
                  + -0.709921*tf.T9 + -2.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  9.39048 + -6.22828*tf.T9i + 2.31435*tf.T913
                  + -0.302835*tf.T9 + 0.020133*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -13.077 + -0.746296*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_O17__F18 = rate

@numba.njit()
def F18__O18__weak__wc12(rate_eval, tf):
    # F18 --> O18
    rate = 0.0

    # wc12w
    rate += np.exp(  -9.15982)

    rate_eval.F18__O18__weak__wc12 = rate

@numba.njit()
def p_O18__He4_N15(rate_eval, tf):
    # O18 + p --> He4 + N15
    rate = 0.0

    # il10r
    rate += np.exp(  10.2725 + -1.663*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -27.9044 + -0.245884*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  26.9671 + -16.6979*tf.T913i
                  + -3.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  8.94352 + -5.32335*tf.T9i + 11.6568*tf.T913
                  + -2.16303*tf.T9 + 0.209965*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_O18__He4_N15 = rate

@numba.njit()
def p_O18__F19(rate_eval, tf):
    # O18 + p --> F19
    rate = 0.0

    # il10r
    rate += np.exp(  -35.0079 + -0.244743*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.917 + -16.7246*tf.T913i
                  + -3.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  7.26876 + -6.7253*tf.T9i + 3.99059*tf.T913
                  + -0.593127*tf.T9 + 0.0877534*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  5.07648 + -1.65681*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_O18__F19 = rate

@numba.njit()
def p_F19__He4_O16(rate_eval, tf):
    # F19 + p --> He4 + O16
    rate = 0.0

    # nacr 
    rate += np.exp(  8.239 + -2.46828*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  -52.7043 + -0.12765*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  26.2916 + -18.116*tf.T913i
                  + 1.86674*tf.T9 + -7.5666*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  14.3586 + -3.286*tf.T9i
                  + -0.21103*tf.T9 + 2.87702*tf.lnT9)
    # nacr 
    rate += np.exp(  15.1955 + -3.75185*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_F19__He4_O16 = rate

@numba.njit()
def He4_O14__p_F17(rate_eval, tf):
    # O14 + He4 --> p + F17
    rate = 0.0

    # Ha96n
    rate += np.exp(  40.8358 + -39.388*tf.T913i + -17.4673*tf.T913
                  + 35.3029*tf.T9 + -24.8162*tf.T953 + -0.666667*tf.lnT9)
    # Ha96r
    rate += np.exp(  16.3087 + -22.51*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  11.1184 + -13.6*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  -106.091 + -0.453036*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  12.1289 + -12.0223*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  18.6518 + -26.0*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_O14__p_F17 = rate

@numba.njit()
def p_F17__Ne18(rate_eval, tf):
    # F17 + p --> Ne18
    rate = 0.0

    # cb09 
    rate += np.exp(  -7.84708 + -0.0323504*tf.T9i + -14.2191*tf.T913i + 34.0647*tf.T913
                  + -16.5698*tf.T9 + 2.48116*tf.T953 + -2.13376*tf.lnT9)
    # cb09 
    rate += np.exp(  27.5778 + -4.95969*tf.T9i + -21.3249*tf.T913i + -0.230774*tf.T913
                  + 0.917931*tf.T9 + -0.0440377*tf.T953 + -7.36014*tf.lnT9)

    rate_eval.p_F17__Ne18 = rate

@numba.njit()
def Ne18__F18__weak__wc12(rate_eval, tf):
    # Ne18 --> F18
    rate = 0.0

    # wc12w
    rate += np.exp(  -0.879336)

    rate_eval.Ne18__F18__weak__wc12 = rate

@numba.njit()
def p_F18__He4_O15(rate_eval, tf):
    # F18 + p --> He4 + O15
    rate = 0.0

    # il10n
    rate += np.exp(  62.0058 + -21.4023*tf.T913i + -80.8891*tf.T913
                  + 134.6*tf.T9 + -126.504*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  1.75704 + -3.01675*tf.T9i + 13.3223*tf.T913
                  + -1.36696*tf.T9 + 0.0757363*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -31.7388 + -0.376432*tf.T9i + 61.738*tf.T913
                  + -108.29*tf.T9 + -34.2365*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_F18__He4_O15 = rate

def rhs(t, Y, rho, T, screen_func=None):
    return rhs_eq(t, Y, rho, T, screen_func)

@numba.njit()
def rhs_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    p_C12__N13(rate_eval, tf)
    p_C13__N14(rate_eval, tf)
    N13__C13__weak__wc12(rate_eval, tf)
    p_N13__O14(rate_eval, tf)
    p_N14__O15(rate_eval, tf)
    p_N15__He4_C12(rate_eval, tf)
    O14__N14__weak__wc12(rate_eval, tf)
    O15__N15__weak__wc12(rate_eval, tf)
    p_N15__O16(rate_eval, tf)
    p_O16__F17(rate_eval, tf)
    F17__O17__weak__wc12(rate_eval, tf)
    p_O17__He4_N14(rate_eval, tf)
    p_O17__F18(rate_eval, tf)
    F18__O18__weak__wc12(rate_eval, tf)
    p_O18__He4_N15(rate_eval, tf)
    p_O18__F19(rate_eval, tf)
    p_F19__He4_O16(rate_eval, tf)
    He4_O14__p_F17(rate_eval, tf)
    p_F17__Ne18(rate_eval, tf)
    Ne18__F18__weak__wc12(rate_eval, tf)
    p_F18__He4_O15(rate_eval, tf)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 1, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C12__N13 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C13__N14 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N13__O14 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N14__O15 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N15__He4_C12 *= scor
        rate_eval.p_N15__O16 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O16__F17 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O17__He4_N14 *= scor
        rate_eval.p_O17__F18 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O18__He4_N15 *= scor
        rate_eval.p_O18__F19 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F19__He4_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O14__p_F17 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F17__Ne18 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F18__He4_O15 *= scor

    dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[jp] = (
       -rho*Y[jp]*Y[jc12]*rate_eval.p_C12__N13
       -rho*Y[jp]*Y[jc13]*rate_eval.p_C13__N14
       -rho*Y[jp]*Y[jn13]*rate_eval.p_N13__O14
       -rho*Y[jp]*Y[jn14]*rate_eval.p_N14__O15
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17
       -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14
       -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       -rho*Y[jp]*Y[jf17]*rate_eval.p_F17__Ne18
       -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15
       +rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17
       )

    dYdt[jhe4] = (
       -rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17
       +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14
       +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15
       )

    dYdt[jc12] = (
       -rho*Y[jp]*Y[jc12]*rate_eval.p_C12__N13
       +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       )

    dYdt[jc13] = (
       -rho*Y[jp]*Y[jc13]*rate_eval.p_C13__N14
       +Y[jn13]*rate_eval.N13__C13__weak__wc12
       )

    dYdt[jn13] = (
       -Y[jn13]*rate_eval.N13__C13__weak__wc12
       -rho*Y[jp]*Y[jn13]*rate_eval.p_N13__O14
       +rho*Y[jp]*Y[jc12]*rate_eval.p_C12__N13
       )

    dYdt[jn14] = (
       -rho*Y[jp]*Y[jn14]*rate_eval.p_N14__O15
       +rho*Y[jp]*Y[jc13]*rate_eval.p_C13__N14
       +Y[jo14]*rate_eval.O14__N14__weak__wc12
       +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14
       )

    dYdt[jn15] = (
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16
       +Y[jo15]*rate_eval.O15__N15__weak__wc12
       +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       )

    dYdt[jo14] = (
       -Y[jo14]*rate_eval.O14__N14__weak__wc12
       -rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17
       +rho*Y[jp]*Y[jn13]*rate_eval.p_N13__O14
       )

    dYdt[jo15] = (
       -Y[jo15]*rate_eval.O15__N15__weak__wc12
       +rho*Y[jp]*Y[jn14]*rate_eval.p_N14__O15
       +rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15
       )

    dYdt[jo16] = (
       -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17
       +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16
       +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       )

    dYdt[jo17] = (
       -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14
       -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18
       +Y[jf17]*rate_eval.F17__O17__weak__wc12
       )

    dYdt[jo18] = (
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       +Y[jf18]*rate_eval.F18__O18__weak__wc12
       )

    dYdt[jf17] = (
       -Y[jf17]*rate_eval.F17__O17__weak__wc12
       -rho*Y[jp]*Y[jf17]*rate_eval.p_F17__Ne18
       +rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17
       +rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17
       )

    dYdt[jf18] = (
       -Y[jf18]*rate_eval.F18__O18__weak__wc12
       -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15
       +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18
       +Y[jne18]*rate_eval.Ne18__F18__weak__wc12
       )

    dYdt[jf19] = (
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       )

    dYdt[jne18] = (
       -Y[jne18]*rate_eval.Ne18__F18__weak__wc12
       +rho*Y[jp]*Y[jf17]*rate_eval.p_F17__Ne18
       )

    return dYdt

def jacobian(t, Y, rho, T, screen_func=None):
    return jacobian_eq(t, Y, rho, T, screen_func)

@numba.njit()
def jacobian_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    p_C12__N13(rate_eval, tf)
    p_C13__N14(rate_eval, tf)
    N13__C13__weak__wc12(rate_eval, tf)
    p_N13__O14(rate_eval, tf)
    p_N14__O15(rate_eval, tf)
    p_N15__He4_C12(rate_eval, tf)
    O14__N14__weak__wc12(rate_eval, tf)
    O15__N15__weak__wc12(rate_eval, tf)
    p_N15__O16(rate_eval, tf)
    p_O16__F17(rate_eval, tf)
    F17__O17__weak__wc12(rate_eval, tf)
    p_O17__He4_N14(rate_eval, tf)
    p_O17__F18(rate_eval, tf)
    F18__O18__weak__wc12(rate_eval, tf)
    p_O18__He4_N15(rate_eval, tf)
    p_O18__F19(rate_eval, tf)
    p_F19__He4_O16(rate_eval, tf)
    He4_O14__p_F17(rate_eval, tf)
    p_F17__Ne18(rate_eval, tf)
    Ne18__F18__weak__wc12(rate_eval, tf)
    p_F18__He4_O15(rate_eval, tf)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 1, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C12__N13 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C13__N14 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N13__O14 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N14__O15 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N15__He4_C12 *= scor
        rate_eval.p_N15__O16 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O16__F17 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O17__He4_N14 *= scor
        rate_eval.p_O17__F18 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O18__He4_N15 *= scor
        rate_eval.p_O18__F19 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F19__He4_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O14__p_F17 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F17__Ne18 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F18__He4_O15 *= scor

    jac = np.zeros((nnuc, nnuc), dtype=np.float64)

    jac[jp, jp] = (
       -rho*Y[jc12]*rate_eval.p_C12__N13
       -rho*Y[jc13]*rate_eval.p_C13__N14
       -rho*Y[jn13]*rate_eval.p_N13__O14
       -rho*Y[jn14]*rate_eval.p_N14__O15
       -rho*Y[jn15]*rate_eval.p_N15__He4_C12
       -rho*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jo16]*rate_eval.p_O16__F17
       -rho*Y[jo17]*rate_eval.p_O17__He4_N14
       -rho*Y[jo17]*rate_eval.p_O17__F18
       -rho*Y[jo18]*rate_eval.p_O18__He4_N15
       -rho*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jf19]*rate_eval.p_F19__He4_O16
       -rho*Y[jf17]*rate_eval.p_F17__Ne18
       -rho*Y[jf18]*rate_eval.p_F18__He4_O15
       )

    jac[jp, jhe4] = (
       +rho*Y[jo14]*rate_eval.He4_O14__p_F17
       )

    jac[jp, jc12] = (
       -rho*Y[jp]*rate_eval.p_C12__N13
       )

    jac[jp, jc13] = (
       -rho*Y[jp]*rate_eval.p_C13__N14
       )

    jac[jp, jn13] = (
       -rho*Y[jp]*rate_eval.p_N13__O14
       )

    jac[jp, jn14] = (
       -rho*Y[jp]*rate_eval.p_N14__O15
       )

    jac[jp, jn15] = (
       -rho*Y[jp]*rate_eval.p_N15__He4_C12
       -rho*Y[jp]*rate_eval.p_N15__O16
       )

    jac[jp, jo14] = (
       +rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       )

    jac[jp, jo16] = (
       -rho*Y[jp]*rate_eval.p_O16__F17
       )

    jac[jp, jo17] = (
       -rho*Y[jp]*rate_eval.p_O17__He4_N14
       -rho*Y[jp]*rate_eval.p_O17__F18
       )

    jac[jp, jo18] = (
       -rho*Y[jp]*rate_eval.p_O18__He4_N15
       -rho*Y[jp]*rate_eval.p_O18__F19
       )

    jac[jp, jf17] = (
       -rho*Y[jp]*rate_eval.p_F17__Ne18
       )

    jac[jp, jf18] = (
       -rho*Y[jp]*rate_eval.p_F18__He4_O15
       )

    jac[jp, jf19] = (
       -rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jhe4, jp] = (
       +rho*Y[jn15]*rate_eval.p_N15__He4_C12
       +rho*Y[jo17]*rate_eval.p_O17__He4_N14
       +rho*Y[jo18]*rate_eval.p_O18__He4_N15
       +rho*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jf18]*rate_eval.p_F18__He4_O15
       )

    jac[jhe4, jhe4] = (
       -rho*Y[jo14]*rate_eval.He4_O14__p_F17
       )

    jac[jhe4, jn15] = (
       +rho*Y[jp]*rate_eval.p_N15__He4_C12
       )

    jac[jhe4, jo14] = (
       -rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       )

    jac[jhe4, jo17] = (
       +rho*Y[jp]*rate_eval.p_O17__He4_N14
       )

    jac[jhe4, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jhe4, jf18] = (
       +rho*Y[jp]*rate_eval.p_F18__He4_O15
       )

    jac[jhe4, jf19] = (
       +rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jc12, jp] = (
       -rho*Y[jc12]*rate_eval.p_C12__N13
       +rho*Y[jn15]*rate_eval.p_N15__He4_C12
       )

    jac[jc12, jc12] = (
       -rho*Y[jp]*rate_eval.p_C12__N13
       )

    jac[jc12, jn15] = (
       +rho*Y[jp]*rate_eval.p_N15__He4_C12
       )

    jac[jc13, jp] = (
       -rho*Y[jc13]*rate_eval.p_C13__N14
       )

    jac[jc13, jc13] = (
       -rho*Y[jp]*rate_eval.p_C13__N14
       )

    jac[jc13, jn13] = (
       +rate_eval.N13__C13__weak__wc12
       )

    jac[jn13, jp] = (
       -rho*Y[jn13]*rate_eval.p_N13__O14
       +rho*Y[jc12]*rate_eval.p_C12__N13
       )

    jac[jn13, jc12] = (
       +rho*Y[jp]*rate_eval.p_C12__N13
       )

    jac[jn13, jn13] = (
       -rate_eval.N13__C13__weak__wc12
       -rho*Y[jp]*rate_eval.p_N13__O14
       )

    jac[jn14, jp] = (
       -rho*Y[jn14]*rate_eval.p_N14__O15
       +rho*Y[jc13]*rate_eval.p_C13__N14
       +rho*Y[jo17]*rate_eval.p_O17__He4_N14
       )

    jac[jn14, jc13] = (
       +rho*Y[jp]*rate_eval.p_C13__N14
       )

    jac[jn14, jn14] = (
       -rho*Y[jp]*rate_eval.p_N14__O15
       )

    jac[jn14, jo14] = (
       +rate_eval.O14__N14__weak__wc12
       )

    jac[jn14, jo17] = (
       +rho*Y[jp]*rate_eval.p_O17__He4_N14
       )

    jac[jn15, jp] = (
       -rho*Y[jn15]*rate_eval.p_N15__He4_C12
       -rho*Y[jn15]*rate_eval.p_N15__O16
       +rho*Y[jo18]*rate_eval.p_O18__He4_N15
       )

    jac[jn15, jn15] = (
       -rho*Y[jp]*rate_eval.p_N15__He4_C12
       -rho*Y[jp]*rate_eval.p_N15__O16
       )

    jac[jn15, jo15] = (
       +rate_eval.O15__N15__weak__wc12
       )

    jac[jn15, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jo14, jp] = (
       +rho*Y[jn13]*rate_eval.p_N13__O14
       )

    jac[jo14, jhe4] = (
       -rho*Y[jo14]*rate_eval.He4_O14__p_F17
       )

    jac[jo14, jn13] = (
       +rho*Y[jp]*rate_eval.p_N13__O14
       )

    jac[jo14, jo14] = (
       -rate_eval.O14__N14__weak__wc12
       -rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       )

    jac[jo15, jp] = (
       +rho*Y[jn14]*rate_eval.p_N14__O15
       +rho*Y[jf18]*rate_eval.p_F18__He4_O15
       )

    jac[jo15, jn14] = (
       +rho*Y[jp]*rate_eval.p_N14__O15
       )

    jac[jo15, jo15] = (
       -rate_eval.O15__N15__weak__wc12
       )

    jac[jo15, jf18] = (
       +rho*Y[jp]*rate_eval.p_F18__He4_O15
       )

    jac[jo16, jp] = (
       -rho*Y[jo16]*rate_eval.p_O16__F17
       +rho*Y[jn15]*rate_eval.p_N15__O16
       +rho*Y[jf19]*rate_eval.p_F19__He4_O16
       )

    jac[jo16, jn15] = (
       +rho*Y[jp]*rate_eval.p_N15__O16
       )

    jac[jo16, jo16] = (
       -rho*Y[jp]*rate_eval.p_O16__F17
       )

    jac[jo16, jf19] = (
       +rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jo17, jp] = (
       -rho*Y[jo17]*rate_eval.p_O17__He4_N14
       -rho*Y[jo17]*rate_eval.p_O17__F18
       )

    jac[jo17, jo17] = (
       -rho*Y[jp]*rate_eval.p_O17__He4_N14
       -rho*Y[jp]*rate_eval.p_O17__F18
       )

    jac[jo17, jf17] = (
       +rate_eval.F17__O17__weak__wc12
       )

    jac[jo18, jp] = (
       -rho*Y[jo18]*rate_eval.p_O18__He4_N15
       -rho*Y[jo18]*rate_eval.p_O18__F19
       )

    jac[jo18, jo18] = (
       -rho*Y[jp]*rate_eval.p_O18__He4_N15
       -rho*Y[jp]*rate_eval.p_O18__F19
       )

    jac[jo18, jf18] = (
       +rate_eval.F18__O18__weak__wc12
       )

    jac[jf17, jp] = (
       -rho*Y[jf17]*rate_eval.p_F17__Ne18
       +rho*Y[jo16]*rate_eval.p_O16__F17
       )

    jac[jf17, jhe4] = (
       +rho*Y[jo14]*rate_eval.He4_O14__p_F17
       )

    jac[jf17, jo14] = (
       +rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       )

    jac[jf17, jo16] = (
       +rho*Y[jp]*rate_eval.p_O16__F17
       )

    jac[jf17, jf17] = (
       -rate_eval.F17__O17__weak__wc12
       -rho*Y[jp]*rate_eval.p_F17__Ne18
       )

    jac[jf18, jp] = (
       -rho*Y[jf18]*rate_eval.p_F18__He4_O15
       +rho*Y[jo17]*rate_eval.p_O17__F18
       )

    jac[jf18, jo17] = (
       +rho*Y[jp]*rate_eval.p_O17__F18
       )

    jac[jf18, jf18] = (
       -rate_eval.F18__O18__weak__wc12
       -rho*Y[jp]*rate_eval.p_F18__He4_O15
       )

    jac[jf18, jne18] = (
       +rate_eval.Ne18__F18__weak__wc12
       )

    jac[jf19, jp] = (
       -rho*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jo18]*rate_eval.p_O18__F19
       )

    jac[jf19, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__F19
       )

    jac[jf19, jf19] = (
       -rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jne18, jp] = (
       +rho*Y[jf17]*rate_eval.p_F17__Ne18
       )

    jac[jne18, jf17] = (
       +rho*Y[jp]*rate_eval.p_F17__Ne18
       )

    jac[jne18, jne18] = (
       -rate_eval.Ne18__F18__weak__wc12
       )

    return jac
