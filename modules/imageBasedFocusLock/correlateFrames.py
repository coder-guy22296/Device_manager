import scipy, pylab
from scipy import signal
from scipy.ndimage import center_of_mass
from scipy.fftpack import fft2,ifft2,fftshift
from Utilities import IO
from Utilities import gaussfitter
from Utilities import scipy_gaussfitter



def corr(frame1, frame2, one_is_ffted=False,
         two_is_ffted=False):
    """
    Parameters
    ----------


    """
    frame_size_x, frame_size_y = frame2.shape

    if not one_is_ffted:
        fft_1 = fft2(frame1)/(frame_size_x*frame_size_y)
    else:
        fft_1 = frame1
    if not two_is_ffted:
        fft_2 = fft2(frame2)/(frame_size_x*frame_size_y)
    else:
        fft_2 = frame2

    correlation = fftshift(ifft2(fft_1*scipy.conj(fft_2)))**2

    #pylab.imshow(abs(correlation))

    return correlation

def corr_high_freq(frame1, frame2, cut_off, one_is_ffted=False):
    """
    Parameters
    ----------


    """
    frame_size_x, frame_size_y = frame2.shape

    x,y = scipy.meshgrid(scipy.arange(-frame_size_x/2,frame_size_x/2),
                         scipy.arange(-frame_size_y/2,frame_size_y/2))

    freqs = scipy.sqrt(x**2 + y**2)

    frame2 = frame2 - frame2.mean()

    if not one_is_ffted:
        fft_1 = fft2(frame1)/(frame_size_x*frame_size_y)
    else:
        fft_1 = frame1
    fft_2 = fft2(frame2)/(frame_size_x*frame_size_y)

    correlation = fftshift(ifft2(fft_1*scipy.conj(fft_2)*(freqs>cut_off)))**2

    return correlation

def fitCorrGaussian(corr, center_region):
    frame_dim_x,frame_dim_y = corr.shape

    center_guess = frame_dim_x/2.0
    height_guess = corr.max()
    noise_guess = corr[1:frame_dim_x/10,1:frame_dim_y/10].mean()

    guess_params = [noise_guess,height_guess-noise_guess,
                    center_region, center_region, 5., 5.]#, 0.]
    

    max_height = height_guess*1.1
    max_width = center_region*1.1
    max_params = [0,max_height,0,0,max_width,max_width]#,360]
    use_max = [False, True, False, False, True, True]#, True]
    corr_region = corr[center_guess-center_region:center_guess+center_region,
                       center_guess-center_region:center_guess+center_region]
    guess_params[2:4] = scipy.unravel_index(corr_region.argmax(), corr_region.shape)
    fits = gaussfitter.gaussfit(corr_region,
                                params=guess_params,
                                maxpars=max_params,
                                limitedmax=use_max,
                                rotate = 0)
    central_val = corr[center_guess-center_region+fits[2]-1:center_guess-center_region+fits[2]+2,
                        center_guess-center_region+fits[3]-1:center_guess-center_region+fits[3]+2].mean()

    return fits, central_val

def fitGaussian1D(vector):
    length = len(vector)

    center_guess = scipy.argmax(vector)
    height_guess = vector.max()
    noise_guess = vector.min()

    guess_params = [noise_guess,height_guess-noise_guess,
                    center_guess, length/6.0]
    

    max_height = height_guess*1.1
    #max_width = center_region*1.1
    #max_params = [0,max_height,0,0,max_width,max_width]#,360]
    #use_max = [False, True, False, False, True, True]#, True]

    fits = gaussfitter.onedgaussfit(scipy.arange(0,length),vector,params=guess_params)

    return fits[0], fits[1]

def fitCorrGaussian_v2(corr, center_region):
    frame_dim_x,frame_dim_y = corr.shape

    center_guess = frame_dim_x/2.0

    corr_region = corr[center_guess-center_region:center_guess+center_region,
                       center_guess-center_region:center_guess+center_region]
    fits = scipy_gaussfitter.fitgaussian(corr_region)
    central_val = corr[center_guess-center_region+fits[1]-1:center_guess-center_region+fits[1]+2,
                        center_guess-center_region+fits[2]-1:center_guess-center_region+fits[2]+2].mean()

    return fits, central_val

def cent(corr, center_region):
    frame_dim_x,frame_dim_y = corr.shape

    center_guess = frame_dim_x/2.0

    corr_region = corr[center_guess-center_region:center_guess+center_region,
                       center_guess-center_region:center_guess+center_region]

    xc,yc = center_of_mass(corr_region - corr_region.min())
    peak = corr_region[xc-2:xc+2, yc-2:yc+2].mean()

    results = scipy.array([peak, xc, yc])

    return results, 0

def polyFitToCorr(zs, corrs, pnum):
    pfit = scipy.polyfit(zs, corrs, pnum)
    evalpfit = scipy.polyval(pfit, zs)
    newzs = scipy.linspace(zs[0], zs[-1], len(zs)*10)
    neweval = scipy.polyval(pfit, newzs)
    w = scipy.argmax(neweval)
    return evalpfit, newzs[w]
