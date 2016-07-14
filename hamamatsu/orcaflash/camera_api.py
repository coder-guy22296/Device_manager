#!/usr/bin/python


from ctypes import *
from ctypes import wintypes
from ctypes.wintypes import WORD, DWORD, BYTE, LONG, SHORT, DOUBLE
import numpy as np
import sys
import Queue
import threading
import time

# Hamamatsu constants.
DCAMCAP_EVENT_FRAMEREADY = int("0x0002", 0)

DCAMERR_NOERROR = 1  # I made this one up. It seems to be the "good" result.

DCAMPROP_ATTR_HASVALUETEXT = int("0x10000000", 0)
DCAMPROP_ATTR_READABLE = int("0x00010000", 0)
DCAMPROP_ATTR_WRITABLE = int("0x00020000", 0)

DCAMPROP_OPTION_NEAREST = int("0x80000000", 0)
DCAMPROP_OPTION_NEXT = int("0x01000000", 0)
DCAMPROP_OPTION_SUPPORT = int("0x00000000", 0)

DCAMPROP_TYPE_MODE = int("0x00000001", 0)
DCAMPROP_TYPE_LONG = int("0x00000002", 0)
DCAMPROP_TYPE_REAL = int("0x00000003", 0)
DCAMPROP_TYPE_MASK = int("0x0000000F", 0)

DCAMWAIT_TIMEOUT_INFINITE = int("0x80000000", 0)

DCAM_CAPTUREMODE_SNAP = 0
DCAM_CAPTUREMODE_SEQUENCE = 1

DCAM_DEFAULT_ARG = 0

DCAM_IDPROP_EXPOSURETIME = int("0x001F0110", 0)

DCAM_IDSTR_MODEL = int("0x04000104", 0)

class DCAM_SIZE(Structure):
    _fields_ = [("cx", LONG),
                ("cy", LONG)]

class DCAM_PARAM_PROPERTYATTR(Structure):
    _fields_ = [("cbSize", LONG),        #Size of this structure
               ("iProp", LONG),         # DCAM ID PROPERTY
               ("option", LONG),        # DAM PROP OPTION
               ("iReserved1", LONG),    #must be 0

               ("attribute", LONG),     #DCAMPROPATTRIBUTE
               ("iGroup", LONG),        #0 (DCAMIDGROUP)
               ("iUnit", LONG),         # DCAMPROPUNIT
               ("attribute2", LONG),    #DCAMPROPATTRIBUTE2
               ("valuemin", DOUBLE),    # minimum value
               ("valuemax", DOUBLE),    # maximum value
               ("valuestep", DOUBLE),   # minimum stepping value
               ("valuedefault", DOUBLE),#default
               ("nMaxChannel", LONG),   #max channel if supports
               ("iReserved3", LONG),    #0
               ("nMaxView", LONG),      #max view if supports
               ("iProp_NumberOfElement", LONG),
               ("iProp_ArrayBase", LONG),
               ("iPropStep_Element", LONG)]


# From Hazen's HAL4000
class HCamData():

    ## __init__
    #
    # Create a data object of the appropriate size.
    #
    # @param size The size of the data object in bytes.
    #
    def __init__(self, size):
        self.np_array = np.ascontiguousarray(np.empty(size/2, dtype=np.uint16))
        self.size = size

    ## __getitem__
    #
    # @param slice The slice of the item to get.
    #
    def __getitem__(self, slice):
        return self.np_array[slice]

    ## copyData
    #
    # Uses the C memmove function to copy data from an address in memory
    # into memory allocated for the numpy array of this object.
    #
    # @param address The memory address of the data to copy.
    #
    def copyData(self, address):
        ctypes.memmove(self.np_array.ctypes.data, address, self.size)

    ## getData
    #
    # @return A numpy array that contains the camera data.
    #
    def getData(self):
        return self.np_array

    ## getDataPtr
    #
    # @return The physical address in memory of the data.
    #
    def getDataPtr(self):
        return self.np_array.ctypes.data

class API:
    '''
    A ctypes based Python wrapper of the Software Development Kit for Cooke pco camera.

    '''

    def __init__(self):
        self.dcam_path = "D:\\Programs\\DCAM-SDK (1112)\\bin\\win32\\"
        self.dcam = windll.LoadLibrary(self.dcam_path+"dcamapi")

        res1 = c_void_p(0)
        count = c_int(0)
        f_init = self.dcam.__getattr__('dcam_init')
        err = f_init(res1, byref(count), wintypes.LPCSTR(0))

        self.f_getmodelinfo = self.dcam.__getattr__('dcam_getmodelinfo')
        self.f_open = self.dcam.__getattr__('dcam_open')
        self.f_close = self.dcam.__getattr__('dcam_close')
        self.f_uninit = self.dcam.__getattr__('dcam_uninit')
        self.f_getstring = self.dcam.__getattr__('dcam_getstring')
        self.f_getcapability = self.dcam.__getattr__('dcam_getcapability')
        self.f_getdatatype = self.dcam.__getattr__('dcam_getdatatype')
        self.f_getbitstype = self.dcam.__getattr__('dcam_getbitstype')
        self.f_setdatatype = self.dcam.__getattr__('dcam_setdatatype')
        self.f_setbitstype = self.dcam.__getattr__('dcam_setbitstype')
        self.f_getdatasize = self.dcam.__getattr__('dcam_getdatasize')
        self.f_getbitssize = self.dcam.__getattr__('dcam_getbitssize')
        self.f_queryupdate = self.dcam.__getattr__('dcam_queryupdate')
        self.f_getbinning = self.dcam.__getattr__('dcam_getbinning')
        self.f_getexposuretime = self.dcam.__getattr__('dcam_getexposuretime')
        self.f_gettriggermode = self.dcam.__getattr__('dcam_gettriggermode')
        self.f_gettriggerpolarity = self.dcam.__getattr__('dcam_gettriggerpolarity')
        self.f_setbinning = self.dcam.__getattr__('dcam_setbinning')
        self.f_setexposuretime = self.dcam.__getattr__('dcam_setexposuretime')
        self.f_settriggermode = self.dcam.__getattr__('dcam_settriggermode')
        self.f_settriggerpolarity = self.dcam.__getattr__('dcam_settriggerpolarity')
        self.f_precapture = self.dcam.__getattr__('dcam_precapture')
        self.f_getdatarange = self.dcam.__getattr__('dcam_getdatarange')
        self.f_getdataframebytes = self.dcam.__getattr__('dcam_getdataframebytes')
        self.f_allocframe = self.dcam.__getattr__('dcam_allocframe')
        self.f_getframecount = self.dcam.__getattr__('dcam_getframecount')
        self.f_capture = self.dcam.__getattr__('dcam_capture')
        self.f_firetrigger = self.dcam.__getattr__('dcam_firetrigger')
        self.f_idle = self.dcam.__getattr__('dcam_idle')
        self.f_wait = self.dcam.__getattr__('dcam_wait')
        self.f_getstatus = self.dcam.__getattr__('dcam_getstatus')
        self.f_gettransferinfo = self.dcam.__getattr__('dcam_gettransferinfo')
        self.f_freeframe = self.dcam.__getattr__('dcam_freeframe')
        self.f_attachbuffer = self.dcam.__getattr__('dcam_attachbuffer')
        self.f_releasebuffer = self.dcam.__getattr__('dcam_releasebuffer')
        self.f_lockdata = self.dcam.__getattr__('dcam_lockdata')
        self.f_lockbits = self.dcam.__getattr__('dcam_lockbits')
        self.f_unlockdata = self.dcam.__getattr__('dcam_unlockdata')
        self.f_unlockbits = self.dcam.__getattr__('dcam_unlockbits')
        self.f_setbitsinputlutrange = self.dcam.__getattr__('dcam_setbitsinputlutrange')
        self.f_setbitsoutputlutrange = self.dcam.__getattr__('dcam_setbitsoutputlutrange')
        self.f_extended = self.dcam.__getattr__('dcam_extended')
        self.f_getlasterror = self.dcam.__getattr__('dcam_getlasterror')
        self.f_getnextpropertyid = self.dcam.__getattr__('dcam_getnextpropertyid')
        self.f_getpropertyname = self.dcam.__getattr__('dcam_getpropertyname')
        self.f_getpropertyvalue = self.dcam.__getattr__('dcam_getpropertyvalue')
        self.f_setpropertyvalue = self.dcam.__getattr__('dcam_setpropertyvalue')

        '''
        Propety IDs:
        bits per channel, 4325680
        trigger source, 1048848
        trigger mode, 1049104
        trigger active, 1048864
        trigger polarity, 1049120
        trigger connector, 1049136
        trigger times, 1049152
        trigger delay, 1049184
        exposture time, 2031888
        defect correct mode, 4653072
        binning, 4198672
        subarray hpos, 4202768
        subarray hsize, 4202784
        subarray vpos, 4202800
        subarray vsize, 4202816
        subarray mode, 4202832  (1==Not used, 2==used)?
        timing readout time, 4206608
        timing cyclic trigger period, 4206624
        timing min trigger blanking, 4206640
        timing min trigger interval, 4206672
        timing exposure, 4206688
        timing invalid exposure period, 4206704
        internal frame rate, 4208656
        internal frame interval, 4208672
        image width, 4325904
        image height, 4325920
        image rowbytes, 4325936
        image framebytes, 4325952
        image top offset bytes, 4325968
        image pixel format, 4326000
        buffer rowbytes, 4326192
        buffer framebytes, 4326208
        buffer top offset bytes, 4326224
        buffer pixel type, 4326240
        number of output trigger connector, 1835024
        output trigger polarity[0], 1835296
        output trigger active[0], 1835312
        output trigger delay[0], 1835328
        output trigger period[0], 1835344
        output trigger kind[0], 1835360
        '''
        
        
        
        # Properties from HAL-4000
        self.frame_bytes = 0
        self.frame_x = 0
        self.frame_y = 0
        self.buffer_index = 0
        self.last_frame_number = 0
        self.max_backlog = 0
        self.number_image_buffers = 0
        
        
        #self.f_ = self.dcam.__getattr__('dcam_')
        

    def openCamera(self):
        self.hcam = wintypes.HANDLE()
        res = wintypes.LPCSTR(0)
        return self.f_open(byref(self.hcam), c_int(0), res)

    def closeCamera(self):
        return self.f_close(self.hcam)

    def getString(self, ID):
        '''
        Model (product name): 0x04000104
        CameraID: 0x04000102
        Camera Version: 0x04000105
        Driver Version: 0x04000106
        Module Version: 0x04000107
        DCAM-API Version: 0x04000108
        '''
        size_string_buffer = 64
        modelinfo = create_string_buffer(size_string_buffer)
        err = self.f_getstring(self.hcam, c_int(ID),
                               byref(modelinfo), DWORD(size_string_buffer))
        return modelinfo.value

    def uninit(self):
        res1 = c_void_p()
        res2 = c_char()
        return self.f_uninit(res1, byref(res2))

    def getDataType(self):
        dataType = c_int()
        err = self.f_getdatatype(self.hcam, byref(dataType))
        return dataType.value

    def getDataSize(self):
        size = DCAM_SIZE()
        err = self.f_getdatasize(self.hcam, byref(size))
        return size.cx, size.cy

    def getStatus(self):
        status = DWORD()
        err = self.f_getstatus(self.hcam, byref(status))
        return status.value

    def getTriggerMode(self):
        pMode = c_int()
        err = self.f_gettriggermode(self.hcam, byref(pMode))
        return pMode.value

    def setTriggerMode(self, value):
        pMode = c_int(value)
        err = self.f_settriggermode(self.hcam, pMode)
        return err

    def getInternalFrameRate(self):
        return self.getPropertyValue(4208656)

    def getExposureTime(self):
        expTime = c_double()
        err = self.f_getexposuretime(self.hcam, byref(expTime))
        return expTime.value

    def setExposureTime(self, expTime):
        return self.f_setexposuretime(self.hcam, DOUBLE(expTime))
        

    def getTransferInfo(self):
        pNewestFrame = c_int()
        pFrameCount = c_int()
        err = self.f_gettransferinfo(self.hcam, byref(pNewestFrame),
                                     byref(pFrameCount))
        return pNewestFrame.value, pFrameCount.value

    def preCapture(self, capturemode):
        ######## Block of code from HAL-4000
        self.buffer_index = -1
        self.last_frame_number = 0

        # Set sub array mode.
        #self.setSubArrayMode()

        # Get frame properties.
        #self.frame_x = self.getPropertyValue(4325904)[0]
        #self.frame_y = self.getPropertyValue(4325920)[0]
        #self.frame_bytes = self.getPropertyValue(4325952)[0]
        ########
        
        #0: Snap; 1: Sequence
        return self.f_precapture(self.hcam, c_int(capturemode))

    def allocFrame(self, framecount):
        return self.f_allocframe(self.hcam, c_int(framecount))

    def getFrameCount(self):
        count = c_int()
        err = self.f_getframecount(self.hcam, byref(count))
        return count.value

    def capture(self):
        return self.f_capture(self.hcam)

    def freeFrame(self):
        return self.f_freeframe(self.hcam)

    def getDataFrameBytes(self):
        size = DWORD()
        err = self.f_getdataframebytes(self.hcam, byref(size))
        return size.value

    def lockData(self, frame, size):
        #Use frame=-1 for last captured frame
        image_temp = np.zeros((size), dtype=np.uint16)
        pTop = c_void_p()
        pRowbytes = c_int()
        err = self.f_lockdata(self.hcam, byref(pTop), byref(pRowbytes), c_int(frame))
        if err!=0:
            memmove(image_temp.ctypes.data, pTop, size*2)
            image = np.frombuffer(image_temp, dtype=np.uint16)
        else:
            image = None
        return image

    def lockBits(self, frame, size):
        image_temp = np.zeros((size), dtype=np.uint8)
        pTop = c_void_p()
        pRowbytes = c_int()
        err = self.f_lockbits(self.hcam, byref(pTop), byref(pRowbytes), c_int(frame))
        memmove(image_temp.ctypes.data, pTop, size)
        image = np.frombuffer(image_temp, dtype=np.uint8)
        return image

    def unlockData(self):
        return self.f_unlockdata(self.hcam)

    def unlockBits(self):
        return self.f_unlockbits(self.hcam)

    def attachBuffer(self, frameSize, numFrames):
        one_frame = np.zeros((frameSize),dtype=np.uint16)
        frames = np.zeros((frameSize*numFrames), dtype=np.uint16)
        dSize = DWORD(frameSize*numFrames*2)
        err = self.f_attachbuffer(self.hcam,
                                  addressof(c_void_p(frames.ctypes.data)),
                                  dSize)
        return (err, frames)

    def attachBuffer2(self, frameSize, numFrames):
        frames = create_string_buffer(frameSize*numFrames*2)
        dSize = DWORD(4*numFrames)
        err = self.f_attachbuffer(self.hcam, byref(frames), dSize)
        return (err, frames)

    def attachBuffer3(self, frameSize, numFrames):
        oneframe = (c_uint16*frameSize)
        frames = oneframe(numFrames)
        #frames = c_void_p(frameSize*numFrames*2)
        dSize = DWORD(numFrames*4)
        err = self.f_attachbuffer(self.hcam, byref(frames), dSize)
        return (err, frames)

    def attachBuffer4(self, frameSize, numFrames):
        self.image_temp = np.zeros((numFrames*frameSize), dtype=np.uint16)
        pTop = c_void_p()
        dSize = DWORD(numFrames*4)
        #memset(byref(pTop), 0, 2*frameSize*numFrames)
        err = self.f_attachbuffer(self.hcam, byref(self.image_temp.ctypes.data_as(c_void_p)), dSize)
        #memmove(image_temp.ctypes.data, pTop, frameSize*numFrames*2)
        return (err,self.image_temp.reshape(numFrames,frameSize))

    def attachBuffer5(self, frameSize, numFrames):
        image_temp = np.zeros((frameSize*numFrames), dtype=np.uint16)
        pTop = (c_void_p * numFrames)()
        memset(byref(pTop), 0, 2*frameSize*numFrames)
        dSize = DWORD(frameSize*numFrames*2)
        err = self.f_attachbuffer(self.hcam, byref(pTop), dSize)
        memmove(image_temp.ctypes.data, pTop, frameSize*numFrames*2)
        return (err, image_temp)
                                  

    def releaseBuffer(self):
        return self.f_releasebuffer(self.hcam)
        
    def idle(self):
        return self.f_idle(self.hcam)

    def wait(self, waitCode, waitTime):
        '''
        waitCode:
        1 = frame start
        2 = frame end
        4 = cycle end
        8 = exposure end
        16 = capture end
        '''
        pCode = DWORD(waitCode)
        dwTime = DWORD(waitTime)
        hdCamSig = wintypes.HANDLE(0)
        return self.f_wait(self.hcam, byref(pCode), dwTime, hdCamSig)

    def getLastError(self):
        size_string_buffer = 64
        errBuff = create_string_buffer(size_string_buffer)
        err = self.f_getlasterror(self.hcam,
                                  byref(errBuff),
                                  DWORD(size_string_buffer))
        return errBuff.value

    def getNextPropertyID(self, iPropNum):
        iProp = c_int(iPropNum)
        err = self.f_getnextpropertyid(self.hcam, byref(iProp), c_int(0))
        return iProp.value

    def getPropertyName(self, iProp):
        size_string_buffer = 64
        propname = create_string_buffer(size_string_buffer)
        err = self.f_getpropertyname(self.hcam, c_int(iProp),
                                     byref(propname), c_int(size_string_buffer))
        return propname.value

    def getPropertyValue(self, iProp):
        pValue = DOUBLE()
        err = self.f_getpropertyvalue(self.hcam, c_int(iProp),
                                      byref(pValue))
        return pValue.value

    def setPropertyValue(self, iProp, value):
        return self.f_setpropertyvalue(self.hcam, c_int(iProp), DOUBLE(value))

    # HAL-4000 Code
    ## newFrames
    #
    # Return a list of the ids of all the new frames since the last check.
    #
    # This will block waiting for at least one new frame.
    #
    # @return [id of the first frame, .. , id of the last frame]
    #
    def newFrames(self):

        # Wait for a new frame.
        dwait = c_int(DCAMCAP_EVENT_FRAMEREADY)
        self.wait(dwait, c_int(DCAMWAIT_TIMEOUT_INFINITE))

        # Check how many new frames there are.
        b_index = c_int32(0)
        f_count = c_int32(0)          
        checkStatus(dcam.dcam_gettransferinfo(self.camera_handle,
                                              ctypes.byref(b_index),
                                              ctypes.byref(f_count)),
                    "dcam_gettransferinfo")
        cur_buffer_index,cur_frame_number = self.getTransferInfo()

        # Check that we have not acquired more frames than we can store in our buffer.
        # Keep track of the maximum backlog.
        backlog = cur_frame_number - self.last_frame_number
        if (backlog > self.number_image_buffers):
            print "warning: hamamatsu camera frame buffer overrun detected!"
        if (backlog > self.max_backlog):
            self.max_backlog = backlog
        self.last_frame_number = cur_frame_number


        # Create a list of the new frames.
        new_frames = []
        if (cur_buffer_index < self.buffer_index):
            for i in range(self.buffer_index + 1, self.number_image_buffers):
                new_frames.append(i)
            for i in range(cur_buffer_index + 1):
                new_frames.append(i)
        else:
            for i in range(self.buffer_index, cur_buffer_index):
                new_frames.append(i+1)
        self.buffer_index = cur_buffer_index

        #if self.debug:
        #    print new_frames

        return new_frames
        
        
    # HAL-4000 Code
    ## getFrames
    #
    # Gets all of the available frames.
    #
    # This will block waiting for new frames even if 
    # there new frames available when it is called.
    #
    # @return [frames, [frame x size, frame y size]]
    #
    def getFrames(self):
        frames = []
        for n in self.newFrames():

            # Lock the frame in the camera buffer & get address.
            data_address = c_void_p(0)
            row_bytes = c_int32(0)
            checkStatus(dcam.dcam_lockdata(self.camera_handle,
                                           ctypes.byref(data_address),
                                           ctypes.byref(row_bytes),
                                           ctypes.c_int32(n)),
                        "dcam_lockdata")

            # Create storage for the frame & copy into this storage.
            hc_data = HCamData(self.frame_bytes)
            hc_data.copyData(data_address)

            # Unlock the frame.
            #
            # According to the documentation, this would be done automatically
            # on the next call to lockdata, but we do this anyway.
            checkStatus(dcam.dcam_unlockdata(self.camera_handle),
                        "dcam_unlockdata")

            frames.append(hc_data)

        return [frames, [self.frame_x, self.frame_y]]
    


if __name__ == '__main__':
    cam_api = API()
    cam_api.openCamera()
    print "Camera: ", cam_api.getString(0x04000104)
    print "Version: ", cam_api.getString(0x04000105)
    print "DCAM-API Version: ", cam_api.getString(0x04000108)
    cam_api.setPropertyValue(4202768, 768)
    cam_api.setPropertyValue(4202784, 256)
    cam_api.setPropertyValue(4202800, 768)
    cam_api.setPropertyValue(4202816, 256)
    cam_api.setPropertyValue(4202832, 2)
    cam_api.setExposureTime(0.0005)
    data_type = cam_api.getDataType()
    xsize,ysize = cam_api.getDataSize()
    print "DataType: ", data_type
    print "x,y size: ", xsize, ysize
    print "Trigger mode: ", cam_api.getTriggerMode()
    print "Exposure time: ", cam_api.getExposureTime()
    cam_api.preCapture(1)
    numFrames = 200
    cam_api.allocFrame(numFrames)
    cam_api.allocFrame(numFrames)
    ims = np.zeros((numFrames,xsize*ysize),dtype=np.uint16)
    print "Frame count for allocFrame: ", cam_api.getFrameCount()
    bytesize = cam_api.getDataFrameBytes()
    print "ByteSize: ", bytesize
    cam_api.capture()
    t1 = time.clock()
    print "Transfer info: ", cam_api.getTransferInfo()
    cam_api.wait(4, 0x80000000)
    print "Transfer info: ", cam_api.getTransferInfo()
    t2 = time.clock()
    for i in range(0, numFrames):
        ims[i] = cam_api.lockData(i, xsize*ysize)
    cam_api.unlockData()
    print "Transfer info: ", cam_api.getTransferInfo()
    t3 = time.clock()
    cam_api.freeFrame()
    #cam_api.preCapture(0)
    
    '''
    err, im2 = cam_api.attachBuffer4(xsize*ysize,8)
    print "AttachBuffer error: ", cam_api.getLastError()
    print "Frame count for attachFrame: ", cam_api.getFrameCount()
    stat = cam_api.getStatus()
    print "Status: ", stat
    
    cam_api.capture()
    stat = cam_api.getStatus()
    print "Status: ", stat
    cam_api.wait(4, 0x80000000)

    time.sleep(1)

    
    #cam_api.wait(16, 0x80000000)
    '''
    #cam_api.wait(16, 0x80000000)
    #time.sleep(2.0)
    #print "Frame count for attachFrame: ", cam_api.getFrameCount()
    #stat = cam_api.getStatus()
    #print "Status: ", stat
    #errRelease = cam_api.releaseBuffer()
    #print "ReleaseBuffer error: ", cam_api.getLastError()
    print "Frame count for attachFrame: ", cam_api.getFrameCount()
    print "Transfer info: ", cam_api.getTransferInfo()
    
    cam_api.idle()

    #errRelease = cam_api.releaseBuffer()
    #print "ReleaseBuffer error: ", cam_api.getLastError()
    
    propID = 0
    for i in range(0,50):
        propID = cam_api.getNextPropertyID(propID)
        propName = cam_api.getPropertyName(propID)
        propVal = cam_api.getPropertyValue(propID)
        print "ID: ", propID
        print "Name: ", propName
        print "Value: ", propVal
    
    #cam_api.freeFrame()
    cam_api.closeCamera()
    cam_api.uninit()
