import numpy as np
import matplotlib.pyplot as plt

def generate_toy_data(nobs, nsim):

    # observed data, or "truth"
    msini_obs = np.random.randn(nobs)*.5
    np.savetxt("observed.txt", np.transpose((msini_obs)))

    # data for simulations (in reality, taken from Kepler)
    radii = .5*np.random.randn(nsim)  # 1000 planets with mean = 1 earth radius
    periods = np.random.uniform(1., 100., nsim)  # days
    star_masses = .5*np.random.randn(nsim) + .3  # M_sun
    inclinations = np.random.uniform(0, 360, nsim)
    data = np.vstack((radii, periods, star_masses, inclinations))
    np.savetxt("sim.txt", data.T)

if __name__ == "__main__":
    generate_toy_data(100, 1000)
