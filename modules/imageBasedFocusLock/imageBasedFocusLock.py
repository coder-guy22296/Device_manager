import inLib
import time
import scipy,pylab,Image
import scipy.misc.pilutil
from scipy import ndimage
import numpy as np
from scipy.ndimage import uniform_filter as uf
from scipy.ndimage import center_of_mass as c_o_m
from scipy.fftpack import fft2,ifft2,fftshift
import os, glob
import correlateFrames
from multiprocessing import Process, Value, Array
import threading
    
gs_ht = 0
xval = 1
yval = 2

#for back of microscope:
xval = 2
yval = 1

reflection_settings = {'zprop': -0.008,
                       'intz': -0.0008,
                       'rolling': 4,
                       'thresh': 0.001}
                       

# For reflection of excitation laser:
#{'zprop': -0.001,
#                       'intz': -0.0001,
#                       'rolling': 4}
image_settings = {'zprop': 0.08,
                  'intz': 0.0001,
                  'rolling': 20}



class Control(inLib.Module):

    def __init__(self, control, settings):
        print "Initializing imageBasedFocusLock..."
        inLib.Module.__init__(self, control, settings)
        print "imageBasedFocusLock initialized."
        
        self.setSettings(settings)
        self.reInit()

        self.acquiredFlatField = False
        
        self.addedSig = 0
        self.lastzs = 1
        self.movedLastTime = [False]

        if settings['camera'] == 'ccd':
            self.useFly = False
        else:
            self.useFly = True
        if self.useFly:
            self.normedMax = 65535.0#65535.0
        else:
            self.normedMax = 255.0 #or 65535

        self.rot90 = False
        
        self.conversions = []
        self.fr_change = 0
        self.cameraCaptureDelay = 0.5 #used in calibration
        self.targetAcquireDelay = 1.0 #used only in getting targets and calibration
        self.monitorDelay = 0.00 #used during focus lock
        self.no_z = False
        self.no_xy = False
        self.number_exposures = 1
        self.exposure_time = 0.2
        self.do_uniform_filter = False
        self.multiplane_step = 0.2
        self.use_multiplane = False
        self.newplane = 1000
        self.zoffset = 0
        self.active = False
        self.quit=False
        self.laserShutOff=False

        self.reflectionMin = 10000

        self.rollingAvXY = False

        self.xystage = self._control.stage
        self.piezo = self._control.piezo
        if self.useFly:
            self.camera = self._control.fly
        else:
            self.camera = self._control.ccd
        self.lasers = self._control.lasers

        self.use_marz=True

        self.micronperpixel_x = -0.01
        self.micronperpixel_y = -0.01
        self.xyprop = 0.09
        self.xyint = 0.0
        self.xyderiv = 0.0

        self.conversion = 0.08
        self.integrated_conv = 0
        self.derivative_conv = 0
        self.threshold_diff = 0.001

        self.conv_reflect = reflection_settings['zprop']
        self.int_reflect = reflection_settings['intz']
        self.thresh_reflect = reflection_settings['thresh']
        self.rolling_reflect = reflection_settings['rolling']

        self.integrateStart = 0
        self.scaleFactor = 1

        self.xythreshold_pixels = 0.05

        self.maybeTroubleThreshold = 25

        self.x = 128
        self.y = 128

        #Upper right corner of ROI
        self.x0 = 0
        self.y0 = 0

        self.num_refs = 3 #could also be 5

        self.below_im = None
        self.target_im = None
        self.above_im = None

        self.above_im2 = None
        self.below_im2 = None

        self.fft_below = None
        self.fft_target = None
        self.fft_above = None
        self.acquiredTargets = False

        self.fft_above2 = None
        self.fft_below2 = None

        self.save_images = True
        self.filenamebase = 'D:\\Data\\Ryan\\image'

        self.dry_run = False
        self.dry_run_z = False

        self.move_every_other = False

        self.saveAllImages = False

        self.findCOM = False


    def reInit(self):
        self.target_signal = []
        self.defocus_signal = []
        self.defocus_signal2 = []
        self.diff_signal = []
        self.previousInitialTargets = []
        self.xdrift = []
        self.ydrift = []
        self.stage_zs = []
        self.stage_xs = []
        self.stage_ys = []
        self.movedz = []
        self.movedx = []
        self.movedy = []
        self.times = []
        self.times2 = []
        self.modelx = []
        self.modely = []
        self.sumx = []
        self.centroidX = []
        return self

    def setReflectionDefaultParams(self):
        self.conversion = self.conv_reflect
        self.integrated_conv = self.int_reflect
        self.threshold_diff = self.thresh_reflect
        self.rolling = self.rolling_reflect
        return [self.conversion, self.integrated_conv,
                self.threshold_diff, self.rolling]

    def getTS(self):
        return self.target_signal

    def flatField(self):
        self.flatfield = self.camera.flatFieldCorrect()
        print self.flatfield
        self.acquiredFlatField = True

    def stopFlatField(self):
        self.acquiredFlatField = False

    def setSettings(self, settings):
        self.gaussFitRegion = settings['xyFitRegion']
        self.target = settings['target']
        self.below = settings['below']
        self.above = settings['above']
        self.rolling = settings['rollingAv']
        self.initialFrames = settings['initialFrames']
        self.unactive_until = settings['unactiveUntil']
        self.frames = settings['frames']
        self.normalize = settings['normalize']
        self.centroid = settings['centroid']

    def setTargets(self, target=None, above=None, below=None):
        if target is not None:
            self.target = target
        if above is not None:
            self.above = above
        if below is not None:
            self.below = below

    def setXYProp(self, xprop, yprop):
        self.micronperpixel_x = xprop
        self.micronperpixel_y = yprop

    def setXYThresh(self, xythresh):
        self.xythreshold_pixels = xythresh

    def setReflectionMin(self, minValue):
        self.reflectionMin = minValue

    def setRotate90(self, state):
        self.rot90 = state

    def setMarz(self, useMarz):
        self.use_marz = useMarz

    def setSaveAll(self, save):
        self.saveAllImages = save

    def setZProp(self, zprop):
        self.conversion = zprop

    def setZInt(self, zint):
        self.integrated_conv = zint

    def setZThresh(self, zthresh):
        self.threshold_diff = zthresh

    def setRolling(self, rolling):
        self.rolling = rolling

    def setRollingXY(self, rolling):
        self.rollingAvXY = rolling

    def setZOffset(self, zoff):
        self.zoffset = zoff

    def getZOffset(self):
        return self.zoffset

    def setFrames(self, frames):
        self.frames = frames

    def _setDryRun(self, dryrun):
        self.dry_run = dryrun

    def setQuit(self, quitState):
        self.quit = quitState

    def setExposureTime(self, exposure):
        self.camera.setExposureTime(exposure)

    def setNoXY(self, state):
        self.no_xy = state

    def setXYDim(self, x, y):
        self.x = x
        self.y = y
        self.camera.newSettings(xdim=x, ydim=y)

    def setX0Y0(self, x0, y0):
        self.camera.newSettings(x0=x0,y0=y0)

    def setSecondImageSize(self, x0,y0,xdim,ydim):
        self.camera.newSettings(xdim_2=xdim,
                                ydim_2=ydim,
                                x0_2=x0,
                                y0_2=y0)

    def setTargetToCurrent(self, zsep):
        zpos = self.piezo.getPosition(3)
        self.target = zpos
        self.above = zpos+zsep
        self.below = zpos-zsep
        return zpos

    def changeCamCaptDelay(self, delay):
        self.cameraCaptureDelay = delay

    def changeMonitorDelay(self, delay):
        self.monitorDelay = delay
        
    def shutDown(self):
        self.quit = True

    def readyCam(self):
        if not self.useFly:
            self.camera.readyCam()

    def setFindCOM(self, state):
        self.findCOM = state

    def getImage(self, secondImage=False):
        if self.useFly:
            im, imB = self.camera.getImage()
        else:
            im,imB = self.camera.getImage(secondary=secondImage)
        if self.acquiredFlatField:
            im = 1.0*im - self.flatfield
        if self.findCOM:
            self.centroidX.append(c_o_m(im)[1])
            #print "Centroid_x: ", self.centroidX[-1]
        if self.saveAllImages:
            num = 1 + len(glob.glob('image*.bmp'))
            scipy.misc.pilutil.imsave('image_'+str(num)+'.bmp', im)
        if self.normalize:
            im2 = (1.0*im)-im.min()
            im = im2 * (self.normedMax/im2.max())
        if self.do_uniform_filter:
            im = uf(im*1.0,[5,5])
        if self.rot90:
            im = np.rot90(im)
        if secondImage:
            if self.normalize:
                im2B = (1.0*imB)-imB.min()
                imB = im2B * (self.normedMax/im2B.max())
            return (1.0*im, 1.0*imB)
        else:
            return 1.0*im

    def moveRelative(self,x,y):
        if self.use_marz:
            self.marz.goRelative(x, y)
        else:
            current_x = self.stage.getPosition(1)
            current_y = self.stage.getPosition(2)
            self.stage.moveTo(1,current_x+x)
            self.stage.moveTo(2,current_y+y)

    def saveImage(self, im, string):
        im_to_save = (1.0*im)-im.min()
        im = im_to_save * (255/(im_to_save.max()))
        im = Image.fromarray(im)
        im = im.convert('L')
        im.save(self.filenamebase+string+'.bmp')
            

    def acquireTargets(self, print_mesg=True):


        if self.no_z:
            time.sleep(5)
            self.target_im = self.getImage()

        else:
        
            #move to below:
            self.piezo.moveTo(3,self.below,waitForConvergence=True)
            time.sleep(self.targetAcquireDelay)
            self.below_im = self.getImage()

            if print_mesg:
                print "Current position",self.piezo.getPosition(3)

            if self.save_images:
                self.saveImage(self.below_im, "_below")

            #move to target
            self.piezo.moveTo(3,self.target,waitForConvergence=True)
            time.sleep(self.targetAcquireDelay)
            self.target_im = self.getImage()

            if print_mesg:
                print "Current position",self.piezo.getPosition(3)

            if self.save_images:
                self.saveImage(self.target_im, "_target")
         
            #move to above
            self.piezo.moveTo(3,self.above,waitForConvergence=True)
            time.sleep(self.targetAcquireDelay)
            self.above_im = self.getImage()

            if print_mesg:
                print "Current position",self.piezo.getPosition(3)

            if self.save_images:
                self.saveImage(self.above_im, "_above")

            self.piezo.moveTo(3,self.target,waitForConvergence=True)
            time.sleep(self.targetAcquireDelay)
            self.piezo.moveTo(3,self.target,waitForConvergence=True)
            time.sleep(self.targetAcquireDelay)

            self.acquiredTargets = True



    def acquireNewTarget(self, newpos, relative=False):
        current_z = self.piezo.getPosition(3)
        if relative:
            self.piezo.moveTo(3, newpos+current_z,waitForConvergence=True)
        else:
            self.piezo.moveTo(3, newpos,waitForConvergence=True)
        time.sleep(self.targetAcquireDelay)
        current_z = self.piezo.getPosition(3)
        self.target_im = self.getImage()
        norm = self.x*self.y
        self.target_im = self.target_im - self.target_im.mean()
        self.fft_target = fft2(self.target_im[0:self.x,0:self.y])/norm
        return current_z
        

    def fftReferences(self):
        norm = self.x*self.y
        
        #Subtract off mean:
        self.target_im = self.target_im - self.target_im.mean()
        if not self.no_z:
            self.below_im = self.below_im - self.below_im.mean()
            self.above_im = self.above_im - self.above_im.mean()
            if self.below_im2 is not None:
                self.below_im2 = self.below_im2 - self.below_im2.mean()
                self.above_im2 = self.above_im2 - self.above_im2.mean()

        self.fft_target = fft2(self.target_im[0:self.x,0:self.y])/norm
        if not self.no_z:
            self.fft_below = fft2(self.below_im[0:self.x,0:self.y])/norm
            self.fft_above = fft2(self.above_im[0:self.x,0:self.y])/norm
            if self.below_im2 is not None:
                self.fft_below2 = fft2(self.below_im2[0:self.x,0:self.y])/norm
                self.fft_above2 = fft2(self.above_im2[0:self.x,0:self.y])/norm

    def saveimage(self,image,num):
        im = scipy.misc.pilutil.toimage(image)
        im.save(self.filenamebase+str(num).zfill(4)+'.bmp')
        

    def getImage_doCorrs(self, num=0):

        if self.centroid==True:
            fit_func = correlateFrames.cent
        else:
            fit_func = correlateFrames.fitCorrGaussian_v2
        corr_func = correlateFrames.corr
        
        xdim_fft, ydim_fft = self.fft_target.shape
        
        current_image = self.getImage()

        xdim, ydim = current_image.shape

        if xdim>xdim_fft or ydim>ydim_fft:
            current_image = current_image[0:xdim_fft,0:ydim_fft]

        current_image = current_image - current_image.mean()
        fft_current_im = fft2(current_image)/(xdim_fft*ydim_fft)
        
        corr_below = corr_func((self.fft_below), fft_current_im, one_is_ffted=True, two_is_ffted=True)
        corr_target = corr_func((self.fft_target), fft_current_im, one_is_ffted=True, two_is_ffted=True)
        corr_above = corr_func(self.fft_above, fft_current_im, one_is_ffted=True, two_is_ffted=True)

        if self.fft_below2 is not None:
            corr_below2 = corr_func(self.fft_below2, fft_current_im, one_is_ffted=True, two_is_ffted=True)
            corr_above2 = corr_func(self.fft_above2, fft_current_im, one_is_ffted=True, two_is_ffted=True)
        
        fit_below,cv = fit_func(abs(corr_below),self.gaussFitRegion)
        fit_target,cv = fit_func(abs(corr_target),self.gaussFitRegion)
        fit_above,cv = fit_func(abs(corr_above),self.gaussFitRegion)

        if self.fft_below2 is not None:
            fit_below2,cv = fit_func(abs(corr_below2), self.gaussFitRegion)
            fit_above2,cv = fit_func(abs(corr_above2), self.gaussFitRegion)
        

        if self.save_images:
            self.saveImage(scipy.real(corr_target),"_corrtarg_"+str(num).zfill(4))
            self.saveImage(current_image, "_current_im_"+str(num).zfill(4))


        if self.num_refs == 3:
            return fit_target[gs_ht], fit_below[gs_ht]-fit_above[gs_ht], fit_target[xval], fit_target[yval]
        elif self.num_refs == 5:
            return fit_target[gs_ht], fit_below[gs_ht]-fit_above[gs_ht], fit_below2[gs_ht]-fit_above2[gs_ht], fit_target[xval], fit_target[yval]

    def doCorrAndFitInThread(self, fit_func, corr_func, ffted_arr1, ffted_arr2, gaussFitRegion, result):
        corr = corr_func(ffted_arr1, ffted_arr2, one_is_ffted=True, two_is_ffted=True)
        fit, cv = fit_func(abs(corr), gaussFitRegion)
        for i in range(0,len(result)):
            result[i] = fit[i]

    def getImage_reflectionBased(self, mask=None):
        current_image = self.getImage()
        #reflectionMin = 12000
        reflectionMin = self.reflectionMin
        if mask is not None:
            current_image = current_image * mask
        sum_x = current_image.sum(axis=0)
        if self.findCOM:
            self.centroidX.append(c_o_m(current_image)[1])
            #print self.centroidX[-1]
            return self.centroidX[-1], 0, 0
        if (sum_x.max()-sum_x.min())>reflectionMin:
            fitresult_x, model_x = correlateFrames.fitGaussian1D(sum_x)
            return fitresult_x[2], model_x, sum_x
        else:
            return 0,0,0

    def justFitPeak(self, sum_x):
        fitresult_x, model_x = correlateFrames.fitGaussian1D(sum_x)
        return fitresult_x[2], model_x

    def getCOMs(self):
        return self.centroidX

    def resetCOMs(self):
        self.centroidX = []

    def captureStack(self, numImages=100):
        im = self.getImage()
        ims = np.zeros((numImages, im.shape[0], im.shape[1]))
        for i in range(0,numImages):
            ims[i] = self.getImage()
        return ims.std(axis=0)

    def createMask(self, stddevs, threshold):
        w = np.where(stddevs < threshold)
        stddevs[w] *= 0.1
        return stddevs

    def getImage_doCorrs_inThread(self, num=0):

        if self.centroid==True:
            fit_func = correlateFrames.cent
        else:
            fit_func = correlateFrames.fitCorrGaussian_v2
        corr_func = correlateFrames.corr
        
        xdim_fft, ydim_fft = self.fft_target.shape
        
        current_image = self.getImage()

        xdim, ydim = current_image.shape

        if xdim>xdim_fft or ydim>ydim_fft:
            current_image = current_image[0:xdim_fft,0:ydim_fft]

        current_image = current_image - current_image.mean()
        fft_current_im = fft2(current_image)/(xdim_fft*ydim_fft)

        ''' Since not doing multiprocessing, don't need this
        result_below = Array('d', range(3))
        result_target = Array('d', range(3))
        result_above = Array('d', range(3))
        '''
        result_below = np.arange(3.)
        result_target = np.arange(3.)
        result_above = np.arange(3.)

        p_below = threading.Thread(target = self.doCorrAndFitInThread,
                          args = (fit_func, corr_func, self.fft_below, fft_current_im,
                                  self.gaussFitRegion, result_below))
        p_target = threading.Thread(target = self.doCorrAndFitInThread,
                           args = (fit_func, corr_func, self.fft_target, fft_current_im,
                                   self.gaussFitRegion, result_target))
        p_above = threading.Thread(target = self.doCorrAndFitInThread,
                          args = (fit_func, corr_func, self.fft_above, fft_current_im,
                                  self.gaussFitRegion, result_above))

        procs = [p_below, p_target, p_above]
        for p in procs:
            p.start()
            p.join()

        fit_target = result_target[:]
        fit_below = result_below[:]
        fit_above = result_above[:]

        #if self.save_images:
        #    self.saveImage(scipy.real(corr_target),"_corrtarg_"+str(num).zfill(4))
        #    self.saveImage(current_image, "_current_im_"+str(num).zfill(4))

        return fit_target[gs_ht], fit_below[gs_ht]-fit_above[gs_ht], fit_target[xval], fit_target[yval]

    def move_calibration(self, xyz, move, repeats, reflectionBased=False):

        gs_ht = 0
        if xyz == 'z':
            current_z = self.piezo.getPosition(3)-self.zoffset
            #Get zoffset:
            zoffset = 0
            for i in range(0,5):
                self.piezo.moveTo(3, current_z)
                time.sleep(0.05)
                zoffset = zoffset - (current_z - self.piezo.getPosition(3))
            self.zoffset = zoffset/5.0
            
        else:
            if not self.use_marz:
                current_x = self.piezo.getPosition(1)
                current_y = self.piezo.getPosition(2)
                print "Current x,y: ", (current_x, current_y)

        if self.num_refs == 3:
            signals = scipy.zeros((2*repeats,4))
        elif self.num_refs == 5:
            signals = scipy.zeros((2*repeats,5))
        if reflectionBased:
            signals = scipy.zeros((2*repeats))
        
        for i in range(0, repeats*2):
            if xyz == 'x':
                if self.use_marz:
                    self.xystage.goRelative(move,0)
                else:
                    self.piezo.moveTo(1, current_x+move, waitForConvergence=False)
            elif xyz == 'y':
                if self.use_marz:
                    self.xystage.goRelative(0, move)
                else:
                    self.piezo.moveTo(2, current_y+move, waitForConvergence=False)
            elif xyz == 'z':
                self.piezo.moveTo(3, current_z + move,waitForConvergence=False)

            time.sleep(self.targetAcquireDelay)

            if xyz=='z':
                current_z = self.piezo.getPosition(3)-self.zoffset

            if reflectionBased:

                xcent,modelx,sumx = self.getImage_reflectionBased()
                signals[i] = xcent

            else:

                if self.num_refs == 3:
                    targ, defo, xdrift, ydrift = self.getImage_doCorrs(num=i)
                    signals[i,0] = targ
                    signals[i,1] = defo
                    signals[i,2] = xdrift
                    signals[i,3] = ydrift
                elif self.num_refs == 5:
                    targ, defo, defo2, xdrift, ydrift = self.getImage_doCorrs(num=i)
                    signals[i,0] = targ
                    signals[i,1] = defo
                    signals[i,2] = defo2
                    signals[i,3] = xdrift
                    signals[i,4] = ydrift

            move *= -1

        return signals

            


    def monitor_limitedFrames_fast(self):
        if self.centroid==True:
            fit_func = correlateFrames.cent
        else:
            fit_func = correlateFrames.fitCorrGaussian_v2
        corr_func = correlateFrames.corr
        
        xdim_fft, ydim_fft = self.fft_target.shape
        i = 0

        print "Dry Run? ", self.dry_run
        
        while (i<self.frames) and (not self.quit):
            self.active=True
            time.sleep(self.monitorDelay)
            self.times.append(time.clock())
            if self.fr_change>0 and i>(self.unactive_until + self.initialFrames-1):
                ind = int(scipy.floor(i/self.fr_change))
                self.conversion = self.conversions[ind]
            
            current_image = self.getImage()
            print "Got image..."

            xdim, ydim = current_image.shape

            if xdim>xdim_fft or ydim>ydim_fft:
                current_image = current_image[0:xdim_fft,0:ydim_fft]
                
            current_image = current_image - current_image.mean()
            fft_current_im = fft2(current_image)/(xdim_fft*ydim_fft)
            
            corr_below = corr_func(self.fft_below, fft_current_im, one_is_ffted=True, two_is_ffted=True)
            corr_target = corr_func(self.fft_target, fft_current_im, one_is_ffted=True, two_is_ffted=True)
            corr_above = corr_func(self.fft_above, fft_current_im, one_is_ffted=True, two_is_ffted=True)
            print "Done correlations..."

            if self.num_refs==5:
                corr_below2 = corr_func(self.fft_below2, fft_current_im, one_is_ffted=True, two_is_ffted=True)
                corr_above2 = corr_func(self.fft_above2, fft_current_im, one_is_ffted=True, two_is_ffted=True)
            
            fit_below,cv = fit_func(abs(corr_below),self.gaussFitRegion)

            self.stage_zs.append(self.piezo.getPosition(3)-self.zoffset)
            
            fit_target,cv = fit_func(abs(corr_target),self.gaussFitRegion)
            fit_above,cv = fit_func(abs(corr_above),self.gaussFitRegion)
            print "Done fitting..."

            if self.num_refs==5:
                fit_below2,cv = fit_func(abs(corr_below2),self.gaussFitRegion)
                fit_above2,cv = fit_func(abs(corr_above2),self.gaussFitRegion)
                self.defocus_signal2.append(fit_below2[gs_ht]-fit_above2[gs_ht])

            self.target_signal.append(fit_target[gs_ht])
            self.defocus_signal.append(fit_below[gs_ht]-fit_above[gs_ht])

            self.xdrift.append(scipy.real(fit_target[xval]) - self.gaussFitRegion)
            self.ydrift.append(scipy.real(fit_target[yval]) - self.gaussFitRegion)

            #print "x, y:", (self.xdrift[-1], self.ydrift[-1])

            
            
            if not self.dry_run:
                if i>(self.unactive_until + self.initialFrames-1):
                    if i==self.unactive_until + self.initialFrames:
                        initial_ds = scipy.mean(self.defocus_signal[self.unactive_until:])
                        initial_z = scipy.mean(self.stage_zs)
                        initial_ts = scipy.mean(self.target_signal[self.unactive_until:])
                        print (initial_ds/initial_ts)
                        print initial_z
                    if i>(self.unactive_until+self.initialFrames):
                        if self.use_multiplane and self.newplane>0:
                            mp = int(i)/int(self.newplane)
                        else:
                            mp = 0
                        self.react_z(initial_ds, initial_ts, initial_z, toprint=False, multiplane=mp)
                        self.react_xy(toprint=False)
            print "Done moving..."

            i = i+1

        if self.laserShutOff:
            self.lasers.shutDown()
        self.safeSave('DefocusSignal', self.defocus_signal)
        self.safeSave('TargetSignal', self.target_signal)
        self.safeSave('Times', self.times)
        self.safeSave('xdrift', self.xdrift)
        self.safeSave('ydrift', self.ydrift)
        self.safeSave('zmovement', self.movedz)
        self.safeSave('zstage', self.stage_zs)
        self.safeSave('xmovement', self.movedx)
        self.safeSave('ymovement', self.movedy)
        self.active = False

    def getQuitState(self):
        return self.quit

    def monitor_limitedFrames_MP(self, i):
        if self.centroid==True:
            fit_func = correlateFrames.cent
        else:
            fit_func = correlateFrames.fitCorrGaussian_v2
        corr_func = correlateFrames.corr
        
        xdim_fft, ydim_fft = self.fft_target.shape
        #i = 0

        #print "Dry Run? ", self.dry_run
        
        if (i<self.frames) and (not self.quit):
            self.active=True
            time.sleep(self.monitorDelay)
            self.times.append(time.clock())
            if self.fr_change>0 and i>(self.unactive_until + self.initialFrames-1):
                ind = int(scipy.floor(i/self.fr_change))
                self.conversion = self.conversions[ind]

            self.stage_zs.append(self.piezo.getPosition(3)-self.zoffset)
            target_gs_ht, defocus_ht_diff, target_xval, target_yval = self.getImage_doCorrs_inThread()
            
            self.target_signal.append(target_gs_ht)
            self.defocus_signal.append(defocus_ht_diff)

            latestSig = defocus_ht_diff / target_gs_ht

            self.xdrift.append(scipy.real(target_xval) - self.gaussFitRegion)
            self.ydrift.append(scipy.real(target_yval) - self.gaussFitRegion)


            if True:
                if i>(self.unactive_until + self.initialFrames-1):
                    if i==self.unactive_until + self.initialFrames:
                        self.initial_ds = scipy.mean(self.defocus_signal[self.unactive_until:])
                        self.initial_z = scipy.mean(self.stage_zs)
                        self.initial_ts = scipy.mean(self.target_signal[self.unactive_until:])
                        self.sig0 = self.initial_ds/self.initial_ts
                        print self.initial_z
                    if i>(self.unactive_until+self.initialFrames):
                        if self.use_multiplane and self.newplane>0:
                            mp = int(i)/int(self.newplane)
                        else:
                            mp = 0
                        if not self.dry_run:
                            self.react_z(self.initial_ds, self.initial_ts, self.initial_z, toprint=False, multiplane=mp)
                            self.react_xy(rolling_av=self.rollingAvXY, toprint=False)
                    initialSigFound = True
                else:
                    initialSigFound = False

            if initialSigFound:
                return self.stage_zs[-1], self.xdrift[-1], self.ydrift[-1], self.sig0, latestSig
            else:
                return self.stage_zs[-1], self.xdrift[-1], self.ydrift[-1], -1, latestSig
        return 0,0,0,0,0

    def changeReflectionTarget(self, target, scale):
        self.initial_ts = np.float64(target)
        self.integrateStart = len(self.target_signal)
        self.scaleFactor = scale

    def changeMaybeTroubleThreshold(self, threshold):
        self.maybeTroubleThreshold = threshold

    def monitor_limitedFrames_Reflection(self, i, timesAndPositions=None):
        
        if (i<self.frames) and (not self.quit):
            self.active=True
            time.sleep(self.monitorDelay)
            self.times.append(time.clock())
            if self.fr_change>0 and i>(self.unactive_until + self.initialFrames-1):
                ind = int(scipy.floor(i/self.fr_change))
                self.conversion = self.conversions[ind]

            self.stage_zs.append(self.piezo.getPosition(3)-self.zoffset)
            xcent,modelx,sumx = self.getImage_reflectionBased()

            if xcent==0:
                return 0,0,0

            self.modelx = modelx
            self.sumx = sumx
            
            self.target_signal.append(xcent)
            self.times2.append(time.clock())
            
            if not self.dry_run:
                if i>(self.unactive_until + self.initialFrames-1):
                    if i==self.unactive_until + self.initialFrames:
                        self.initial_ts = scipy.mean(self.target_signal[self.unactive_until:])
                        self.initial_z = scipy.mean(self.stage_zs)
                        print self.initial_z
                    self.previousInitialTargets.append(self.initial_ts)
                    if i>(self.unactive_until+self.initialFrames):
                        self.react_reflection(self.initial_ts)
                    initialSigFound = True
                   
                    #print "initial_ts ", self.initial_ts
                else:
                    initialSigFound = False

        timenow = self.times[-1] - self.times[0]
        
        return xcent,modelx,self.stage_zs[-1], timenow

    def returnModels(self):
        return self.modelx,self.sumx

    def finishUpFocusLock(self, reflection=False):

        if self.laserShutOff:
            self.lasers.shutDown()

        focusLockDirectory = 'focuslock\\'
        if not os.path.exists(focusLockDirectory):
            os.makedirs(focusLockDirectory)

        flDir = os.getcwd()+'\\focuslock\\'
        
        
        self.safeSave(flDir+'TargetSignal', self.target_signal)
        self.safeSave(flDir+'Times_', self.times)
        self.safeSave(flDir+'zmovement', self.movedz)
        self.safeSave(flDir+'zstage', self.stage_zs)
        if not reflection:
            self.safeSave(flDir+'DefocusSignal', self.defocus_signal)
            self.safeSave(flDir+'xmovement', self.movedx)
            self.safeSave(flDir+'ymovement', self.movedy)
            self.safeSave(flDir+'xdrift', self.xdrift)
            self.safeSave(flDir+'ydrift', self.ydrift)
        else:
            self.safeSave(flDir+'TimesReacted', self.times2)
        self.active = False

    def safeSave(self, filenm, data):
        if os.path.exists(filenm+'.npy') or os.path.exists(filenm+str(0).zfill(3)+'.npy'):
            num = len(glob.glob(filenm+'*'))
            filenm = filenm+str(num).zfill(3)
        else:
            filenm = filenm+str(0).zfill(3)
        np.save(filenm, data)

    def moveZ(self, zpos):
        self.stage.moveTo(3, zpos)
        self.acquiredTargets = False

    def moveZRelative(self, zchange):
        zpos = self.stage.getPosition(3)
        self.stage.moveTo(3, zpos+zchange)
        self.acquiredTargets = False

    def react_reflection(self, initial, rolling_av=True):
        maybeTrouble = False
        if rolling_av:
            weights = scipy.exp((-1.*(scipy.arange(self.rolling,0,-1.)/self.rolling)**2)/2.)
            weighted_ts_mean = scipy.average(self.target_signal[(-1*self.rolling):], weights=weights)
            diff = weighted_ts_mean - initial
        else:
            diff = self.target_signal[-1] - initial
            
        self.diff_signal.append(diff)

        integrated_diff = scipy.sum(self.diff_signal)

        scaleFactor=1

        if len(self.previousInitialTargets)>100:
            if self.previousInitialTargets[-100]==initial:
                scaleFactor = 1
            else:
                scaleFactor = self.scaleFactor

        if len(self.stage_zs)>10:
            diffLast10 = np.array(self.stage_zs[-10:])-np.array(self.stage_zs[-11:-1])
            zchanged = diffLast10.sum()
            zmoved = abs(diffLast10).sum()
            ratio = np.array(zmoved)/np.array(zchanged)
            if ratio>self.maybeTroubleThreshold:
                maybeTrouble = True
            else:
                maybeTrouble = False


        if abs(diff)>self.threshold_diff:
            z_change = ((diff*self.conversion) + (integrated_diff*self.integrated_conv)) * scaleFactor
            if not self.dry_run_z and not maybeTrouble:
                last_z = self.stage_zs[-1]
                self.piezo.moveTo(3, scipy.mean(self.stage_zs[(-1*self.lastzs):]) + z_change, waitForConvergence=False)
                self.movedz.append(z_change)
            else:
                self.movedz.append(0)
        else:
            self.movedz.append(0)
                

    def react_z(self, initial_ds, initial_ts, init_z, rolling_av=True, toprint=True, multiplane=0):
        if multiplane>0:
            newsig = multiplane*self.multiplane_step
        else:
            newsig = 0
        if rolling_av:
            sig_to_add = newsig + self.addedSig
            weights = scipy.exp((-1.*(scipy.arange(self.rolling,0,-1.)/self.rolling)**2)/2.)
            weighted_ds_mean = scipy.average(self.defocus_signal[(-1*self.rolling):], weights=weights)
            weighted_ts_mean = scipy.average(self.target_signal[(-1*self.rolling):], weights=weights)
            diff = (weighted_ds_mean/weighted_ts_mean) - ((initial_ds/initial_ts)+sig_to_add)
        else:
            diff = (self.defocus_signal[-1]/self.target_signal[-1]) - (initial_ds/initial_ts)
        inst_diff = (self.defocus_signal[-1]/self.target_signal[-1]) - ((initial_ds/initial_ts)+newsig)
        last_slope_1 = self.defocus_signal[-1] - self.defocus_signal[-2]
        last_slope_2 = self.defocus_signal[-2] - self.defocus_signal[-3]
        integrated_diff = scipy.sum((np.array(self.defocus_signal)/np.array(self.target_signal)) - (initial_ds/initial_ts))

        conversion = self.conversion
        
        #if abs(diff) > abs(last_slope_1):
        #    conversion *= 0.9
        
        if abs(diff) > self.threshold_diff:  #2012/08/16: changed to 'diff' from 'inst_diff'
            z_change = (diff * conversion) + (integrated_diff * self.integrated_conv) + (last_slope_1*self.derivative_conv)
            if (not self.movedLastTime[-1]) or (not self.move_every_other):
                if self.dry_run_z:
                    print "Dry run. Would have moved",z_change
                else:
                    if toprint:
                        print "Z Movement: ", z_change
                    if self.move_every_other:
                        last_z = scipy.mean(self.stage_zs[-2:-1])
                    else:
                        last_z = self.stage_zs[-1]
                    self.piezo.moveTo(3, scipy.mean(self.stage_zs[(-1*self.lastzs):]) + z_change, waitForConvergence=False)
                    self.movedz.append(z_change)
                self.movedLastTime.append(True)
            else:
                self.movedLastTime.append(False)
                self.movedz.append(0)
        else:
            self.movedLastTime.append(False)
            self.movedz.append(0)

    def react_xy(self, rolling_av=False, toprint=True):

        if rolling_av:
            weights = scipy.exp((-1.*(scipy.arange(self.rolling,0,-1.)/self.rolling)**2)/2.)
            xd = scipy.average(self.xdrift[(-1*self.rolling):],weights=weights)
            yd = scipy.average(self.ydrift[(-1*self.rolling):],weights=weights)
        else:
            xd = self.xdrift[-1]
            yd = self.ydrift[-1]

        if len(self.xdrift)>1:
            last_slope_1_x = self.xdrift[-1] - self.xdrift[-2]
            last_slope_1_y = self.ydrift[-1] - self.ydrift[-2]

        integrated_diff_x = scipy.sum(self.xdrift)
        integrated_diff_y = scipy.sum(self.ydrift)

        move_x =  xd * self.micronperpixel_x
        move_y =  -1*yd * self.micronperpixel_y

        if not self.use_marz:
            last_x = self.piezo.getPosition(1)
            last_y = self.piezo.getPosition(2)

        
        if (not self.movedLastTime[-1]) or (not self.move_every_other):
            if (abs(xd) > self.xythreshold_pixels) and (not self.no_xy):
                if self.use_marz:
                    if toprint:
                        print "Moving x:", move_x
                    self.xystage.goRelative(move_x,0)
                    self.movedx.append(move_x)
                else:
                    if toprint:
                        print "Moving x:", move_x
                    self.piezo.moveTo(1, last_x+move_x, waitForConvergence=False)
                    self.movedx.append(move_x)
            else:
                self.movedx.append(0)
            if (abs(yd) > self.xythreshold_pixels) and (not self.no_xy):
                if self.use_marz:
                    if toprint:
                        print "Moving y:", move_y
                    self.xystage.goRelative(0,move_y)
                    self.movedy.append(move_y)
                else:
                    if toprint:
                        print "Moving y:", move_y
                    self.piezo.moveTo(2, last_y+move_y, waitForConvergence=False)
                    self.movedx.append(move_y)
            else:
                self.movedy.append(0)

'''
class FocusLockThread(threading.Thread):
    def __init__(self, focusLockParams, focusLock):
        threading.Thread.__init__(self)

        self.centroid = focusLockParams['centroid']
        self.frames = focusLockParams['frames']
        self.fft_target = focusLockParams['fft_target']
        self.unactive_until = focusLockParams['unactive_until']
        self.initialFrames = focusLockParams['initialFrames']
        self.monitorDelay = focusLockParams['monitorDelay']
        self.conversion = focusLockParams['conversion']
        self.dry_run = focusLockParams['dry_run']
        self.fft_below = focusLockParams['fft_below']
        self.fft_above = focusLockParams['fft_above']
        self.num_refs = focusLockParams['num_refs']

        self.active=False
        
        self.times = []
        self.stage_zs = []

        self.stoprequest = threading.Event()
        

    def run(self):
        if self.centroid==True:
            fit_func = correlateFrames.cent
        else:
            fit_func = correlateFrames.fitCorrGaussian_v2
        corr_func = correlateFrames.corr
        
        xdim_fft, ydim_fft = self.fft_target.shape
        i = 0

        print "Dry Run? ", self.dry_run
        
        while not self.stoprequest.isSet():
            self.active=True
            time.sleep(self.monitorDelay)
            self.times.append(time.clock())

            
            current_image = self.getImage()

            xdim, ydim = current_image.shape

            if xdim>xdim_fft or ydim>ydim_fft:
                current_image = current_image[0:xdim_fft,0:ydim_fft]
                
            current_image = current_image - current_image.mean()
            fft_current_im = fft2(current_image)/(xdim_fft*ydim_fft)
            
            corr_below = corr_func(self.fft_below, fft_current_im, one_is_ffted=True, two_is_ffted=True)
            corr_target = corr_func(self.fft_target, fft_current_im, one_is_ffted=True, two_is_ffted=True)
            corr_above = corr_func(self.fft_above, fft_current_im, one_is_ffted=True, two_is_ffted=True)

            if self.num_refs==5:
                corr_below2 = corr_func(self.fft_below2, fft_current_im, one_is_ffted=True, two_is_ffted=True)
                corr_above2 = corr_func(self.fft_above2, fft_current_im, one_is_ffted=True, two_is_ffted=True)
            
            fit_below,cv = fit_func(abs(corr_below),self.gaussFitRegion)

            self.stage_zs.append(self.piezo.getPosition(3)-self.zoffset)
            
            fit_target,cv = fit_func(abs(corr_target),self.gaussFitRegion)
            fit_above,cv = fit_func(abs(corr_above),self.gaussFitRegion)

            if self.num_refs==5:
                fit_below2,cv = fit_func(abs(corr_below2),self.gaussFitRegion)
                fit_above2,cv = fit_func(abs(corr_above2),self.gaussFitRegion)
                self.defocus_signal2.append(fit_below2[gs_ht]-fit_above2[gs_ht])

            self.target_signal.append(fit_target[gs_ht])
            self.defocus_signal.append(fit_below[gs_ht]-fit_above[gs_ht])

            self.xdrift.append(scipy.real(fit_target[xval]) - self.gaussFitRegion)
            self.ydrift.append(scipy.real(fit_target[yval]) - self.gaussFitRegion)

            #print "x, y:", (self.xdrift[-1], self.ydrift[-1])

            
            
            if not self.dry_run:
                if i>(self.unactive_until + self.initialFrames-1):
                    if i==self.unactive_until + self.initialFrames:
                        initial_ds = scipy.mean(self.defocus_signal[self.unactive_until:])
                        initial_z = scipy.mean(self.stage_zs)
                        initial_ts = scipy.mean(self.target_signal[self.unactive_until:])
                        print (initial_ds/initial_ts)
                        print initial_z
                    if i>(self.unactive_until+self.initialFrames):
                        if self.use_multiplane and self.newplane>0:
                            mp = int(i)/int(self.newplane)
                        else:
                            mp = 0
                        self.react_z(initial_ds, initial_ts, initial_z, toprint=False, multiplane=mp)
                        self.react_xy(toprint=False)

            i = i+1
'''
        
        
