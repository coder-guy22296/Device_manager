''' Parameters '''
stabilization_time = 1
conventional_frames = 20
bleaching_time = 1
STORM_frames = 1000

# Make sure that the PFS is calibrated:
if not control.scope.getFocalPlaneCalibration():
    print 'Script canceled: Calibrate Nikon PFS first!'


else:
    # The PFS seems to be calibrated, so let's go:


    ''' General setup '''

    # Turn off lasers:
    control.hal4000.toggleChannel(0, False)   # 0 = 642
    control.hal4000.toggleChannel(4, False)   # 4 = 405
    # Set laser powers:
    control.hal4000.setPower(0, 0.5)
    control.hal4000.setPower(4, 0.0)
    # Import a Python library to support timings:
    import time
    # A list of hotspots:
    hotSpots = control.hotSpots.getHotSpots(optimize=True,
        priority='scope.PFS')
    # The number of hotSpots:
    nHotSpots = len(hotSpots)
    # Hal4000's working directory:
    workingDir = control.hal4000.getWorkingDirectory()

    # Loop through all hotspots:
    for i in range(nHotSpots):

        # Check if user clicked 'Stop'. If yes, break this loop:
        if control.scripts.stop:
            break

        # Set up the movie base name as in 'movie_0001':
        movieNameBase = 'movie_{0:04d}'.format(i+1)
        # Move to hotspot i:
        hotSpots[i].moveTo()

        # Wait for stabilization:
        time.sleep(stabilization_time)

        if control.scripts.stop:
            break

        ''' Conventional fluorescence image '''

        # Choose Hal4000's settings file:
        control.hal4000.parameters(0)
        control.hal4000.setRunShutters(False)
        control.motorizer.setOD(4)
        control.hal4000.toggleChannel(0, True)
        # Record movie:
        control.hal4000.movie(movieNameBase + '_642', conventional_frames)
        control.hal4000.toggleChannel(0, False)
        
        if control.scripts.stop:
            break

        ''' Out-of-focus bleaching '''

        control.motorizer.setOD(0)
        # Open TIRF ('--' = -0.1, '-' = -0.02):
        control.motorizer.stepTIRF('--')
        control.hal4000.toggleChannel(0, True)
        # Wait 1 second:
        time.sleep(bleaching_time)
        control.hal4000.toggleChannel(0, False)
        control.motorizer.stepTIRF('++')

        if control.scripts.stop:
            break

        ''' STORM image '''

        control.hal4000.parameters(1)
        control.hal4000.setRunShutters(True)
        control.hal4000.movie(movieNameBase, STORM_frames)
        # Get focal plane position:
        z = control.scope.getFocalPlanePosition()
        # Append focal plane position to movie's .inf file:
        with open(workingDir + movieNameBase + '.inf', 'a') as infFile:
            infFile.write('Focal plane position = {0}'.format(z))
