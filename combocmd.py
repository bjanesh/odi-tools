#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

title_string = 'test'

r_mag, r_ierr, gmr = np.loadtxt('calibrated_mags_long.dat', usecols=(-3,-2,-1), unpack=True)
r_mag_sh, r_ierr_sh, gmr_sh = np.loadtxt('calibrated_mags.dat', usecols=(-3,-2,-1), unpack=True)


plt.clf()

fig = plt.figure(figsize=(10, 8)) 
fig.subplots_adjust(hspace=0)
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
ax0 = plt.subplot(gs[0])

ax0.scatter(gmr, r_mag, s=2, color='black', marker='o', edgecolors='none')
ax0.scatter(gmr_sh, r_mag_sh, s=2, color='red', marker='o', edgecolors='none')

ax0.set_ylabel('$r$')
ax0.set_xlabel('$(g-r)$')
ax0.set_ylim(22,12)
ax0.set_xlim(-1,2)

ax1 = plt.subplot(gs[1])

ax1.scatter(r_ierr, r_mag, s=2, color='black', marker='o', edgecolors='none')
ax1.scatter(r_ierr_sh, r_mag_sh, s=2, color='red', marker='o', edgecolors='none')
ax1.vlines([0.01, 0.02,0.03,0.04],12,22, linestyle='--', colors='gray')

# ax1.set_ylabel('$r$')
ax1.set_xlabel('inst $r$ err.')
ax1.set_ylim(22,12)
ax1.set_xlim(-0.002,0.05)
plt.setp(ax1.get_yticklabels(), visible=False)

plt.tight_layout()
plt.savefig(title_string+"_CMD.pdf")