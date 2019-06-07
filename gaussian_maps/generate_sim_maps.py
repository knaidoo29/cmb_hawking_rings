import numpy as np
import healpy as hp

# Calculated from command line running of CAMB-0.1.7 the default base_planck_lowl_lowLike.ini file
data = np.loadtxt('cls/test_lensedCls.dat', unpack=True)
l = data[0]
cl = data[1]

# cls are multiplied by l(l+1)/(2pi) so need to get rid of this constant [units are muK^2]
cl *= 2.*np.pi/(l*(l+1.))

# need to add cls for l = 0, 1 (i.e. l=0 => constant background and l=1 => dipole)
l = np.concatenate([np.array([0., 1.]), l])
cl = np.concatenate([np.zeros(2), cl])

nside, lmax = 2048, 1500

for i in range(0, 1000):
    map_sim = hp.synfast(cl, nside, lmax=lmax)
    filename = 'maps/map_' + str(i+1) + '.fits'
    hp.write_map(filename, map_sim)
