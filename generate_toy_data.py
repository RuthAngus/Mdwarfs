import numpy as np
import matplotlib.pyplot as plt
from toymodel import mass2radius

def generate_toy_data(nobs, nsim):

    # observed data, or "truth"
    m_obs = np.random.randn(nobs)*.1 + 2
    i = np.random.uniform(0, np.pi, nobs)
    np.savetxt("observed.txt", np.transpose((m_obs*np.sin(i))))

    # data for simulations (in reality, taken from Kepler)
    m_for_sim = np.random.randn(nsim)*.1 + 2
    truth = [.9, .2]
    radii = mass2radius(truth, m_for_sim)
    periods = np.random.uniform(1., 100., nsim)  # days
    star_masses = .01*np.random.randn(nsim) + .5  # M_sun
    inclinations = np.random.uniform(0, np.pi, nsim)
    data = np.vstack((radii, periods, star_masses, inclinations))
    np.savetxt("sim.txt", data.T)

if __name__ == "__main__":
    generate_toy_data(1000, 1000)
