# A dumb test of everything stupid and uncertain.
# Later this will be shaped into a good simulator of optics system. 
from zen_calc import *
from dask.optimize import unwrap_partial
from skimage.restoration import unwrap_phase
import glob
import numpy as np

path1 = '/home/sillycat/Documents/Light_sheet/Data/Jun09/'
path2 = '/home/sillycat/Documents/Light_sheet/Data/May17/'

# PSF1 = np.load(path1+'T0_modnone.npy')
# PSF2 = np.load(path2+ 'psfm01.npy')

class SPIM_specs():
    

    dx = 0.097
    dz = 0.30
    l = 0.550   
    n = 1.33
    NA = 1.00
    f = 9000
    nIt = 4

# PSF = PSF2
path = path1

ii=1
plt.close('all')


for pfile in PSF_list30:
    PSF = np.load(pfile)
    PF = retrievePF(PSF, dx, dz, l, n, NA, f, nIt)
    PF = unwrap_phase(PF)
    plt.figure(figsize=(4.5,4))
#     plt.imshow(PF, cmap = 'RdBu', )
    im = plt.imshow(PF, cmap = plt.cm.RdBu, extent=(-2,2,2,-2))
    plt.tick_params(
        axis = 'both',
        which = 'both', 
        bottom = 'off',
        top = 'off',
        right = 'off',
        left = 'off',
        labelleft='off',
        labelbottom = 'off')


    plt.tight_layout()
    plt.savefig(pfile[:-4]+'pf')



plt.show()

