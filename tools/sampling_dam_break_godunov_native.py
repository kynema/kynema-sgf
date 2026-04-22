#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
KYNEMA_SGF_PATH = '.'
sys.path.append(KYNEMA_SGF_PATH+'/tools/')
import amrex_particle
from amrex_particle import AmrexParticleFile

loc_dir = "."
pp_dir = loc_dir + "/post_processing"

pfile = AmrexParticleFile(loc_dir)

nt = 5 + 1
out_int = 30

plt.figure()
base_color = 0.8

for n in range(nt):
    pt = pfile.load(n * out_int, root_dir = pp_dir)
    pt.parse_header()
    pt.load_binary_data()
    data = pt.df
    x_oo = data.xco
    z_oo = data.zco
    u_oo = data.velocityx
    v_oo = data.velocityy
    w_oo = data.velocityz
    # Reorder arrays
    ind_0 = np.argsort(x_oo) 
    x = np.zeros(np.shape(x_oo))
    z = np.zeros(np.shape(x_oo))
    Vmag = np.zeros(np.shape(x_oo))
    flag = np.zeros(np.shape(x_oo))
    for i in range(0, len(ind_0)): 
        x[i]= x_oo[ind_0[i]]
        z[i]= z_oo[ind_0[i]]
        Vmag[i] = np.sqrt(u_oo[ind_0[i]]**2+v_oo[ind_0[i]]**2+w_oo[ind_0[i]]**2)
        # When interface is not present, z location is set to zlo
        flag[i] = 1.0 if z[i] > 1e-8 else 0.
    # Shorten arrays (exclude points where interface not detected)
    nvalid = int(flag.sum())
    xshort = x[0:nvalid]
    zshort = z[0:nvalid]
    vshort = Vmag[0:nvalid]

    color = base_color - (base_color) * (n/nt)
    cstr = str(color)

    plt.plot(xshort,zshort,color=cstr)
    plt.scatter(xshort,zshort,c=vshort,cmap="jet",vmin=0.,vmax=2.)

plt.ylabel(r'$z$',fontsize=16)
plt.xlabel(r'$x$',fontsize=16)
plt.colorbar()
plt.savefig('plot_sampling_native.pdf',format='pdf',dpi=300,bbox_inches="tight")