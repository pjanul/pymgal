import numpy as np
__ver__ = '1.0'


class charlot_fall(object):
    """ callable-object implementation of the Charlot and Fall (2000) dust law """
    tau1 = 0.0
    tau2 = 0.0
    tbreak = 0.0

    def __init__(self, tau1=1.0, tau2=0.5, tbreak=0.01):
        """ dust_obj = charlot_fall(tau1=1.0, tau2=0.3, tbreak=0.01)
        Return a callable object for returning the dimming factor as a function of age
        for a Charlot and Fall (2000) dust law.  The dimming is:

        np.exp(-1*Tau(t)(lambda/5500angstroms))

        Where Tau(t) = `tau1` for t < `tbreak` (in gyrs) and `tau2` otherwise. """

        self.tau1 = tau1
        self.tau2 = tau2
        self.tbreak = tbreak

    def __call__(self, ts, ls):
        ls1 = np.copy(ls)
        ts1 = np.copy(ts)
        ls1 = np.asarray(ls1)
        ts1 = np.asarray(ts1)
        ls1.shape = (ls1.size, 1)
        ts1.shape = (1, ts1.size)

        taus = np.asarray([self.tau1] * ts1.size)
        m = (ts1 > self.tbreak).ravel()
        if m.sum():
            taus[m] = self.tau2

        return np.exp(-1.0 * taus * (ls1 / 5500.0)**-0.7)
