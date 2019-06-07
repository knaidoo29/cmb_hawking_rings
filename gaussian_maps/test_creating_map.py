import numpy as np
import matplotlib.pylab as plt
import healpy as hp

# Calculated from command line running of CAMB-0.1.7 the default base_planck_lowl_lowLike.ini file
data = np.loadtxt('cls/test_lensedCls.dat', unpack=True)
l = data[0]
cl = data[1]

# Sanity check...
plt.figure(figsize=(8, 6))
plt.plot(l, cl)
plt.xlabel(r'$l$', fontsize=16)
plt.ylabel(r'$l(l+1)C_{l}/(2\pi)$ $[(\mu K)^{2}]$', fontsize=16)
plt.show()

# cls are multiplied by l(l+1)/(2pi) so need to get rid of this constant [units are muK^2]
cl *= 2.*np.pi/(l*(l+1.))

# need to add cls for l = 0, 1 (i.e. l=0 => constant background and l=1 => dipole)
l = np.concatenate([np.array([0., 1.]), l])
cl = np.concatenate([np.zeros(2), cl])

# see https://healpy.readthedocs.io/en/latest/generated/healpy.sphtfunc.synfast.html#healpy.sphtfunc.synfast
# for additional parameters (such as adding a Gaussian smoothing, etc)

nside, lmax = 2048, 1500
map_sim = hp.synfast(cl, nside, lmax=lmax)

# saving maps as a fits file using healpy function, of course you could save this in anyway you want. filsize = 192Mb
filename = 'maps/test_map.fits'
hp.write_map(filename, map_sim)

# opening maps (checking this works)
map_sim = hp.read_map(filename)

# plot simulated map
hp.mollview(map_sim)
plt.show()

# measure cls from the simulated map
cl_sim = hp.anafast(map_sim, lmax=lmax)
cl_sim *= (l[:len(cl_sim)]*(l[:len(cl_sim)]+1.))/(2.*np.pi)
cl *= (l*(l+1.))/(2.*np.pi)

# plot input cls vs cls from simulated map
plt.figure(figsize=(8, 6))
plt.plot(l, cl)
plt.plot(l[:len(cl_sim)], cl_sim, color='grey', alpha=0.5)
plt.xlabel(r'$l$', fontsize=16)
plt.ylabel(r'$l(l+1)C_{l}/(2\pi)$ $[(\mu K)^{2}]$', fontsize=16)
plt.show()
