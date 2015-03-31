import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
import scipy.stats as sps

def load_data():
    msini_obs = np.genfromtxt("observed.txt").T
    r_sim, p_sim, m_star, inc = np.genfromtxt("sim.txt").T
    return msini_obs, r_sim, p_sim, m_star, inc

# convert the radii of the imaginary planets into masses
def radius2mass(theta, radius):
    m, c = theta
    return m * radius + c

# convert the radii of the imaginary planets into masses
def mass2radius(theta, mass):
    m, c = theta
    return (mass - c) / m

def rv_K(P, m_star, m_planet, inc):
    m_star_kg = m_star * 1.9891e30
    m_planet_kg = m_planet * 5.972e24
    """
    Calculate RV semi-amplitude of body 1 for orbital period P and mass
    function f, = m1 (sin(incl))**2 / (m1+m2)**3, where m1 and m2 are
    the masses of bodies 1 and 2, respectively, and incl is the
    inclination of the orbital plane (incl = pi/2 for a system seen
    edge-on). Takes period in days, m_star in solar masses, m_planet in
    Earth masses.
    """
    f = m_planet_kg * np.sin(inc) / (m_star_kg + m_planet)**(2/3.)
    return (2 * np.pi * 6.67e-11 / (P*86400))**(1/3.) * f

# create a population of the synthesised planets that would be detected
def populationize(theta, sim_data):
    r_sim, p_sim, m_star, inc = sim_data
    m_planet = radius2mass(theta, r_sim)  # calculate mass
    K = rv_K(p_sim, m_star, m_planet, inc)  # calculate semi-amplitude

    # detected?
    l = (K > .1) * (p_sim < 100)
    return m_planet[l] * np.sin(inc[l])  # return msinis

# compute the sufficient statistic, the KS statistic of the mass distributions
def KS_statistic(m_obs, m_sim):
    return sps.ks_2samp(m_obs, m_sim)[0]  # ks_stat, pvalue

# decide whether to accept or reject
def accept_reject(theta, msini_obs, sim_data, tolerance):
    msini_sim = populationize(theta, sim_data)  # simulated msinis

    # accept/reject
    if KS_statistic(msini_obs, msini_sim) < tolerance:
        return theta
    else:
        return np.ones_like(theta) * np.inf

# generate ABC posterior
# wrapper: generate parameters, compute KS, ask if rejected, keep if accepted
def ABC_sampling(msini_obs, sim_data, tolerance, nsamp):

    # generate parameter samples
    thetas = np.zeros((2, nsamp))
#     thetas[0, :] = np.random.randn(nsamp)*.01 + 1  # slope
#     thetas[1, :] = np.random.randn(nsamp)*.5  # intercept
    thetas[0, :] = np.random.uniform(0., 3, nsamp)  # slope
    thetas[1, :] = np.random.uniform(-1., 3, nsamp)  # intercept
    keepsies = np.zeros_like(thetas)

    # determine whether samples are accepted or rejected
    for i in range(nsamp):
        keepsies[:, i] = accept_reject(thetas[:, i], msini_obs, sim_data,
                                       tolerance)
    l = np.isfinite(keepsies[0, :])
    parameters = np.zeros((2, len(keepsies[0, :][l])))
    parameters[0] = keepsies[0, :][l]
    parameters[1] = keepsies[1, :][l]

    plt.clf()
    plt.subplot(2, 1, 1)
#     plt.hist(thetas[0, :], alpha=.5, color="r")
    plt.hist(parameters[0], alpha=.5, color="b")
    plt.subplot(2, 1, 2)
#     plt.hist(thetas[1, :], alpha=.5, color="r")
    plt.hist(parameters[1], alpha=.5, color="b")
    plt.show()
    return parameters

if __name__ == "__main__":

    # load data
    msini_obs, r_sim, p_sim, m_star, inc = load_data()
    sim_data = np.vstack((r_sim, p_sim, m_star, inc))

    # do some ABC
    tolerance, nsamp = .2, 10000
    parameters = ABC_sampling(msini_obs, sim_data, tolerance, nsamp)
