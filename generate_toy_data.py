import numpy as np
import matplotlib.pyplot as plt
from toymodel import mass2radius

def generate_toy_data(nobs, nsim):

    # observed data, or "truth" - these will be the observed msinis
    m_obs = np.random.randn(nobs)*.1 + 2
    inclinations = np.random.uniform(0, np.pi, nobs)
    np.savetxt("observed.txt", np.transpose((m_obs*np.sin(inclinations))))

    # data for simulations (in reality, these will be radius and period
    # measurements taken from Kepler. The inclinations will be random.
    truth = [1.4, 1.5]  # slope, intercept
    radii = mass2radius(truth, m_obs)
    periods = np.random.uniform(1., 100., nsim)  # days
    star_masses = .01*np.random.randn(nsim) + .5  # M_sun
    data = np.vstack((radii, periods, star_masses, inclinations))
    np.savetxt("sim.txt", data.T)

if __name__ == "__main__":
    generate_toy_data(10000, 10000)
