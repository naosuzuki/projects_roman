import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions, CircleSkyRegion, RegionMeta, RegionVisual
from astropy.coordinates import SkyCoord
import astropy.units as u


wcs = WCS(fits.getheader('/Users/conor/Desktop/IRAC/NEP.linear.1.mosaic.fits'))
data = fits.getdata('/Users/conor/Desktop/IRAC/NEP.linear.1.mosaic.fits')

center = SkyCoord('17:58:55.9 +66:01:03.7', unit='hourangle,deg')
Xc,Yc = wcs.world_to_pixel(center)
print(Xc, Yc)
sX, sY = wcs.proj_plane_pixel_scales()

survey_radius = np.sqrt(10./np.pi)
extended_survey_radius = np.sqrt(20./np.pi)

radius = 1.1 * u.degree
leaves = 6
Theta = np.linspace(0, -2*np.pi, leaves+1)[:-1]
Theta += np.pi/2. # Rotate leaves by 90 deg

ext_radius = 2.05 * u.degree # fiducial
ext_leaves = 9
ext_Theta = np.linspace(0, -2*np.pi, ext_leaves+1)[:-1] + np.radians(10)

EDFN_radius = np.sqrt(20/np.pi)/wcs.pixel_scale_matrix[1,1]

innner_leaf_offsets = [(radius*np.cos(t), radius*np.sin(t)) for t in Theta]
outer_leaf_offsets = [(ext_radius*np.cos(t), ext_radius*np.sin(t)) for t in ext_Theta]
t_offset = 2*np.pi / ext_leaves / 2
outer_leaf_offsets += [(ext_radius*np.cos(t+t_offset), ext_radius*np.sin(t+t_offset)) for t in ext_Theta]


inner_leaf_pixcoords = [((x/sX).value + Xc, (y/sY).value + Yc) for x,y in innner_leaf_offsets]
outer_leaf_pixcoords = [((x/sX).value + Xc, (y/sY).value + Yc) for x,y in outer_leaf_offsets]


leaf_pixcoords = inner_leaf_pixcoords + outer_leaf_pixcoords
#print(leaf_pixcoords)

#exit()
leaf_pos = wcs.all_pix2world(leaf_pixcoords, 0)
leaf_coords = SkyCoord(leaf_pos[:,0], leaf_pos[:,1], unit='deg')
# print(leaf_coords)

vis = lambda i: RegionVisual(color='Lime' if i < 6 else 'Red', width=2, fontsize=20)
meta = lambda i: RegionMeta(text=str(i if i <= 6 else i-6))
regs = Regions( [CircleSkyRegion(center, radius=0.75*u.degree, meta=meta(0), visual=vis(0))] +
    [CircleSkyRegion(c, radius=0.75*u.degree, 
                     meta=meta(i+1),
                     visual=vis(i)) 
     for i, c in enumerate(leaf_coords)]
)

# H20_selfcal = SkyCoord(268.5295950, 65.1203827, unit='deg')
# regs.append(
#     CircleSkyRegion(
#         H20_selfcal, 0.75*u.degree, 
#         meta=RegionMeta(text='H20\nSelfcal'), 
#         visual=RegionVisual(color='Cyan', width=2)
#     )
# )
H20_selfcal = SkyCoord(268.813,	+65.29, unit='deg')
regs.append(
    CircleSkyRegion(
        H20_selfcal, np.sqrt(2.5/np.pi)*u.degree, 
        meta=RegionMeta(text=r'Selfcal'), 
        visual=RegionVisual(color='Cyan', width=2)
    )
)

regs.write('EDF-N_extended_pointings.reg', format='ds9', overwrite=True)


leaf_pixcoords = np.array(leaf_pixcoords)

import matplotlib.pyplot as plt
from astropy.visualization import simple_norm

plt.rcParams['font.family'] = 'serif'
plt.rcParams['text.usetex'] = True
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13

hsc_pixrad = (0.75/sX).value

norm = simple_norm(data, stretch='log', min_cut=-0, max_cut=5)

box = dict(boxstyle='circle', facecolor='w', alpha=0.3, pad=0.1)

plt.figure(figsize=(8,8))
plt.subplot(111,projection=wcs)
#plt.imshow(data, norm=norm, cmap='Greys_r')
for i,r in enumerate(regs):
    pixr = r.to_pixel(wcs)
    if r != regs[-1]:
        pixr.plot()
        if i > 6:
            plt.text(*pixr.center.xy, r'\textbf{'+r.meta['text']+'}', fontsize=40, 
                    ha='center', va='center', color=r.visual['color'])
    else:
        continue
        pixr.plot(facecolor='w', alpha=0.5, fill=True)
        pixr.plot()
        x,y = pixr.center.xy
        plt.text(x, y+1300, r'\textbf{Euclid}', fontsize=40, 
                ha='center', va='center', color=r.visual['color'], zorder=10)
        plt.text(x, y-1300, r'\textbf{Selfcal}', fontsize=40, 
                ha='center', va='center', color=r.visual['color'], zorder=10)

t = np.linspace(0,2*np.pi,100)
plt.plot(Xc + EDFN_radius*np.cos(t), Yc + EDFN_radius*np.sin(t), 'k--', lw=2)

plt.xlim(leaf_pixcoords[:,0].min() - hsc_pixrad, leaf_pixcoords[:,0].max() + hsc_pixrad)
plt.ylim(leaf_pixcoords[:,1].min() - hsc_pixrad, leaf_pixcoords[:,1].max() + hsc_pixrad)
plt.xlabel('RA', fontsize=20)
plt.ylabel('Dec', fontsize=20)
#plt.axis('off')

plt.tight_layout()
plt.show()