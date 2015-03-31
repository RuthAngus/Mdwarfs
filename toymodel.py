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

def rv_K(P, m_star, m_planet, inc):
    m_star *= const.M_sun
    m_planet *= const.M_earth
    """
    Calculate RV semi-amplitude of body 1 for orbital period P and mass
    function f, = m1 (sin(incl))**2 / (m1+m2)**3, where m1 and m2 are
    the masses of bodies 1 and 2, respectively, and incl is the
    inclination of the orbital plane (incl = pi/2 for a system seen
    edge-on). Takes period in days, m_star in solar masses, m_planet in
    Earth masses.
    """
    f = m_star * np.sin(inc) ** 2 / (m_star + m_planet) ** 3
    return (2 * np.pi * 6.67e-11 / (P*86400)**2)**(1/3.) * f

# Given the K of each planet, calculate whether it would be
# detected by an RV survey (simplified for now)
def detected(K, P):
    l = (K > 10) * (P < 200)
    return l

# create a population of the synthesised planets that would be detected
def populationize(theta, sim_data):
    r_sim, p_sim, m_star, inc = sim_data
    m_planet = radius2mass(theta, r_sim)  # calculate mass
    K = rv_K(p_sim, m_star, m_planet, inc)  # calculate semi-amplitude
    l = detected(K, p_sim)  # throw away undetectable planets
    return m_planet[l] * np.sin(inc[l])  # return msinis

# compute the sufficient statistic, the KS statistic of the mass distributions
def KS_statistic(m_obs, m_sim):
    return sps.ks_2samp(m_obs, m_sim)[0]  # ks_stat, pvalue

# decide whether to accept or reject
def accept_reject(theta, msini_obs, sim_data, tolerance):
    msini_sim = populationize(theta, sim_data)  # simulated msinis
    if KS_statistic(msini_obs, msini_sim) < tolerance:
        return theta
    else:
        return np.ones_like(theta) * np.inf

# generate ABC posterior
# wrapper: generate parameters, compute KS, ask if rejected, keep if accepted
def ABC_sampling(msini_obs, sim_data, tolerance, nsamp):

    # generate parameter samples
    thetas = np.zeros((2, nsamp))
    thetas[0, :] = 2*np.random.randn(nsamp)*.1  # gradient
    thetas[1, :] = .5*np.random.randn(nsamp)*.1  # intercept
    keepsies = np.zeros_like(thetas)

    # determine whether samples are accepted or rejected
    for i in range(nsamp):
        keepsies[:, i] = accept_reject(thetas[:, i], msini_obs, sim_data,
                                       tolerance)
    l = np.isfinite(keepsies[0, :])
    parameters = np.zeros((2, len(keepsies[0, :][l])))
    parameters[0] = keepsies[0, :][l]
    parameters[1] = keepsies[1, :][l]
    return parameters

# approximate the posterior from the pooled parameters
def ABC_posterior(parameters):
    plt.clf()
    plt.subplot(2, 1, 1)
    plt.hist(parameters[0])
    plt.subplot(2, 1, 2)
    plt.hist(parameters[1])
    plt.savefig("ABC_posterior")

if __name__ == "__main__":
    msini_obs, r_sim, p_sim, m_star, inc = load_data()
    plt.clf()
    plt.hist(msini_obs)
    plt.show()
    assert 0
    sim_data = np.vstack((r_sim, p_sim, m_star, inc))
    ABC_sampling(msini_obs, sim_data, tolerance, nsamp)
