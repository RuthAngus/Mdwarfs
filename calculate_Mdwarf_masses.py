import numpy as np
import matplotlib.pyplot as plt
from isochrones.dartmouth import Dartmouth_Isochrone
from isochrones import StarModel
import time
import h5py

# find a better mass for an individual star
def colour2mass(B, J, K, R, V, MCMC=False):
    mags = {'B': (B, 0.08),
            'J': (J, 0.08),
            'K': (K, 0.08),
            'R': (R, 0.08),
            'V': (V, 0.08)}
    # find the mass and radius using isochrones
    dar = Dartmouth_Isochrone()  # Tim's isochrone instance
    smod_phot = StarModel(dar, **mags)
    if MCMC:
        result = smod_phot.fit_mcmc()
        # save samples and triangle plot
        smod_phot.triangle_plots("test", format="pdf")
        smod_phot.save_hdf("samples.h5")
    else:
        result = smod_phot.maxlike()  # calculate maximum likelihood
    return np.array(result)

def get_stellar_parameters(min_distance=0, max_distance=10):
    # load data
    B, R, J, K, V, J_V, plx = np.genfromtxt("data/superblink.txt",
                                          skip_header=50, skip_footer=3,
                                          invalid_raise=False,
                                          usecols=(0, 1, 2, 3, 4, 5, 6)).T
    st = np.genfromtxt("data/superblink.txt", skip_header=50, skip_footer=3,
                       invalid_raise=False, usecols=(7), dtype=str).T

    D = 1./plx
    l = (min_distance < D) * (D < max_distance)  # volume limited sample

    start = time.time()
    star_params = np.zeros((len(B[l]), 5))  # mass, age, feh, distance A_V
    print "time expected ~ %.2f" % (25*len(B[l])/60.), "minutes"
    for i in range(len(B[l])):
        print i, "of", len(B[l])
        star_params[i, :] = colour2mass(B[l][i], J[l][i], K[l][i], R[l][i],
                                        V[l][i])

    end = time.time()
    print "time taken = %.2f" % ((end-start)/60.), "minutes"

    print "saving results"
    f = h5py.File("Mdwarf_parameters.h5", "w")
    data = f.create_dataset("params", (np.shape(star_params)))
    data[:, 0] = star_params[:, 0]
    data[:, 1] = star_params[:, 1]
    data[:, 2] = star_params[:, 2]
    data[:, 3] = star_params[:, 3]
    data[:, 4] = star_params[:, 4]
    f.close()

if __name__ == "__main__":

    get_stellar_parameters(max_distance=10)
#     with h5py.File("Mdwarf_parameters.h5", "r") as f:
#         parameters = f["params"][:, :]
#     masses = parameters[:, 0]
#     plt.clf()
#     plt.hist(masses, 50)
#     plt.show()
