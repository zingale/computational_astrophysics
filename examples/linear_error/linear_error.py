import numpy as np

def f(x):
    return x * np.sin(x)

def sample_func(f, npts, xmin=0.0, xmax=5.0):
    x = np.linspace(xmin, xmax, npts)
    return x, f(x)

def interp(x0, x, f):

    # find the index where x0 lies
    idx = max(0, np.argwhere(x > x0)[0][0] - 1)

    slope = (f[idx+1] - f[idx]) / (x[idx+1] - x[idx])

    return slope * (x0 - x[idx]) + f[idx]

x0 = 1.0

npts = 20
xv, fv = sample_func(f, npts)

print(f(x0), interp(x0, xv, fv))

