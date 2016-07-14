import numpy as np
import libtim.zern
import matplotlib.pyplot as plt
from scipy import optimize
import pupil2device as pupil
from numpy.lib.scimath import sqrt as _msqrt

# single zernike mode 
def singleZernike(mode, amp, radius=None, useMask=True):
        # modes: a scalar
    modes = np.zeros((mode))
    modes[mode-1]=amp
    zernike = libtim.zern.calc_zernike(modes, radius, mask=useMask,
                                                zern_data = {})
    return zernike
#     calc_zern_basis(nmodes, rad, modestart=1, calc_covmat=False):

def multiZernike(amps, rad, useMask = True):
        #amps: parameter as modes in calc_zernike
        zernike = libtim.zern.calc_zernike(amps, rad, mask = useMask, zern_data={})
        return zernike



def plotZernike(zernike):
    plt.imshow(zernike, interpolation='nearest')
    plt.show()
    
def retrievePF(PSF, dx, dz, l, n, NA, f, nIt, bscale=1.0, neglect_defocus=True,
                   invert=False, wavelengths=1, resetAmp=False,
                   symmeterize=False):

    z_offset = 0
    if neglect_defocus:
            # We try to estimate the axial position of the emitter
            # The coordinates of the brightest pixel
        cz, cx, cy = np.unravel_index(PSF.argmax(), PSF.shape)
            # Intensity trace along z
        i = PSF[:,cx,cy]
            # Get z positions
        nz = PSF.shape[0]
        upper = 0.5*(nz-1)*dz
        z = np.linspace(-upper, upper, nz)
            # Initial fit parameters
        b = np.mean((i[0],i[-1]))
        a = i.max() - b
        w = l/3.2
        p0 = (a,0,w,b)
        def gaussian(z, a, z0, w, b):
            return a * np.exp(-(z-z0)**2/w) + b
            # Fit gaussian to axial intensity trace
        popt, pcov = optimize.curve_fit(gaussian, z, i, p0)
            # Where we think the emitter is axially located:
        z_offset = -1.0*popt[1] #Added on April 3, 2015
#         plt.plot(z, i)
#         plt.plot(z, gaussian(z,*popt))
#         plt.savefig('z_fit.png')
    
    nx,ny = PSF.shape[1:3]
    PF = pupil.Simulation(nx,dx,l,n,NA,f,wavelengths=wavelengths)
#         self, nx=256, dx=0.1, l=0.68, n=1.33, NA=1.27, f=3333.33, wavelengths=10, wavelength_halfmax=0.005
    A = PF.plane
    
    Mx, My = np.meshgrid(np.arange(nx)-nx/2., np.arange(nx)-nx/2.)
    r_pxl = _msqrt(Mx**2 + My**2)
    
    hcyl = np.array(nz*[np.logical_and(r_pxl>=65, r_pxl<76)])
    background = np.mean(PSF[hcyl])*bscale
    print "Finding PF..."
    print "   Using parameters:"
    print "   dz = ", dz
    print "   background = ", background
    print "   z_offset = ", z_offset
    complex_PF = PF.psf2pf(PSF, dz, background, A, nIt, z_offset,
                                        resetAmp=resetAmp,symmeterize=symmeterize)
#       def psf2pf(self, PSF, dz, mu, A, nIterations=10, z_offset=0, use_pyfftw=True, resetAmp=True,
#                symmeterize=False):

    if invert:
        complex_PF = abs(complex_PF) * np.exp(-1*1j*np.angle(complex_PF))

    Pupil_final = _PupilFunction(complex_PF, PF)

    nz = PSF.shape[0]
    upper = 0.5*(nz-1)*dz
    z = np.linspace(-upper, upper, nz)
    psft = PF.pf2psf(complex_PF, z)
#         def psf2pf(self, PSF, dz, mu, A, nIterations=10, z_offset=0, use_pyfftw=True, resetAmp=True,
#                symmeterize=False):
    np.save('retrieved_psf', psft)

    return Pupil_final.phase


class _PupilFunction(object):
    '''
    A pupil function that keeps track when when either complex or amplitude/phase
    representation is changed.
    '''
    def __init__(self, complex, geometry):
        self.complex = complex
        self._geometry = geometry

    @property
    def complex(self):
        return self._complex

    @complex.setter
    def complex(self, new):
        self._complex = new
        self._amplitude = abs(new)
        self._phase = np.angle(new)

    @property
    def amplitude(self):
        return self._amplitude

    @amplitude.setter
    def amplitude(self, new):
        self._amplitude = new
        self._complex = new * np.exp(1j*self._phase)

    @property
    def phase(self):
        return self._phase

    @phase.setter
    def phase(self, new):
        self._phase = new
        self._complex = self._amplitude * np.exp(1j*new)

    @property
    def zernike_coefficients(self):
        return self._zernike_coefficients

    @zernike_coefficients.setter
    def zernike_coefficients(self, new):
        self._zernike_coefficients = new
        #self._zernike = zernike.basic_set(new, self._geometry.r, self._geometry.theta)
        self._zernike = libtim.zern.calc_zernike(new, self._geometry.nx/2.0)

    @property
    def zernike(self):
        return self._zernike