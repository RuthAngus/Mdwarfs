import numpy as np
import matplotlib.pyplot as plt
from gatspy.periodic import LombScargle
import glob

def Xmatch(fnames):
    names = np.genfromtxt("../data/RVplanetHIP.txt")
    names = np.array([int(names[i]) for i in range(len(names))])
    match = []
    for i, fname in enumerate(fnames):
        l = names==fname.isdigit
        if sum(l) > 0:
            match.append(names[names==fname.isdigit])
    print match

def periodogram(fs, t, y, dy):
    model = LombScargle().fit(t, y, dy)
    period = 1. / fs
    power = model.periodogram(period)
    return power

def peak_detect(x, y):
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    l = y[peaks] == max(y[peaks])
    return x[peaks], y[peaks], x[peaks][l], y[peaks][l]

"""
For each star without a planet they compute, on a fine grid of orbital
periods, the maximum mass planet compatible with the RV measurements.
What does that actually mean? Do they inject and recover?
Do they compute a periodogram?
"""

if __name__ == "__main__":

    fnames = glob.glob("../data/vels/*")

    fs = np.arange(1./100, 1/.1, .01)
    print len(fs)

    for fname in fnames:
        print fname
        t, y, dy = np.genfromtxt(fname).T
        pgram = periodogram(fs, t, y, dy)

        fpeaks, ppeaks, mf, mp = peak_detect(fs, pgram)

        plt.clf()
        plt.plot(fs, pgram)
        plt.plot(mf, mp, "ro")
        plt.show()
        raw_input('enter')
