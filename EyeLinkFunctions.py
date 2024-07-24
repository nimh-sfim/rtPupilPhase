import pylink
import os
from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy
from string import ascii_letters, digits
from psychopy import core
import sys
import platform
import time
import config

from PsychoPy_funcs import instructions_screens, terminate_task

def validate_edf_fname(edf_fname):
    # Note: The script below is provided by SR Research, Inc.
    edf_fname = edf_fname.rstrip().split('.')[0]
    time_str = time.strftime("_%Y_%m_%d_%H_%M", time.localtime())

    allowed_char = ascii_letters + digits + '_'
    if not all([c in allowed_char for c in edf_fname]):
        state=False
        message="ERROR: *** Invalid EDF filename"
    elif len(edf_fname) > 8: 
        state=False
        message="ERROR: *** EDF filename should not exceed 8 characters"
    else: 
        state=True
        message="Valid filename"

    return edf_fname, state, message

def setup_eyelink(win, dummy_mode, edf_fname): 

    # Note: The script below is provided by SR Research, Inc.
     # Step 1: Connect to the EyeLink Host PC

    # The Host IP address, by default, is "100.1.1.1".
    # the "el_tracker" objected created here can be accessed through the Pylink
    # Set the Host PC address to "None" (without quotes) to run the script
    # in "Dummy Mode"
    if dummy_mode:
        el_tracker = pylink.EyeLink(None)
    else:
        try:
            el_tracker = pylink.EyeLink("100.1.1.1")
        
        except RuntimeError as error:
            print('ERROR:', error)
            core.quit()
            sys.exit()
        
    # Step 2: Open an EDF data file on the Host PC

    # Define edf filename
    edf_file = edf_fname + ".EDF"

    try:
        el_tracker.openDataFile(edf_file)

    except RuntimeError as err:    
        print('ERROR:', err)
        
        # Close the link if we have one open
        if el_tracker.isConnected():
            el_tracker.close()
        core.quit()
        sys.exit()

    # Step 3: Configure the tracker

    # Put the tracker in offline mode before we change tracking parameters
    el_tracker.setOfflineMode()

    # Get the software version:  1-EyeLink I, 2-EyeLink II, 3/4-EyeLink 1000,
    # 5-EyeLink 1000 Plus, 6-Portable DUO
    if dummy_mode:
        eyelink_ver = 0  # set version to 0, in case running in Dummy mode
        
    else:
        eyelink_ver = config.eyelink_ver
        
    if not dummy_mode:
        vstr = el_tracker.getTrackerVersionString()
        eyelink_ver = int(vstr.split()[-1].split('.')[0])
        
        # Print out some version info in the shell
        print('Running experiment on %s, version %d' % (vstr, eyelink_ver))

    # File and link data control
    # What eye events to save in the EDF file, include everything by default
    file_event_flags = 'LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT'

    # What eye events to make available over the link, include everything by default
    link_event_flags = 'LEFT,RIGHT,FIXATION,SACCADE,BLINK,BUTTON,FIXUPDATE,INPUT'

    # What sample data to save in the EDF data file and to make available
    # over the link, include the 'HTARGET' flag to save head target sticker
    # data for supported eye trackers
    if eyelink_ver > 3:
        file_sample_flags = 'LEFT,RIGHT,GAZE,HREF,RAW,AREA,HTARGET,GAZERES,BUTTON,STATUS,INPUT'
        link_sample_flags = 'LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT'

    else:
        file_sample_flags = 'LEFT,RIGHT,GAZE,HREF,RAW,AREA,GAZERES,BUTTON,STATUS,INPUT'
        link_sample_flags = 'LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT'

    if not dummy_mode:
        
        el_tracker.sendCommand("file_event_filter = %s" % file_event_flags)
        el_tracker.sendCommand("file_sample_data = %s" % file_sample_flags)
        el_tracker.sendCommand("link_event_filter = %s" % link_event_flags)
        el_tracker.sendCommand("link_sample_data = %s" % link_sample_flags)
        # Set EyeLink sample rate
        if eyelink_ver > 2 :
            el_tracker.sendCommand("sample_rate 1000")
        
        # Choose a calibration type, H3, HV3, HV5, HV13 (HV = horizontal/vertical),
        el_tracker.sendCommand("calibration_type = HV9")

        # Set a gamepad button to accept calibration/drift check target
        # You need a supported gamepad/button box that is connected to the Host PC
        el_tracker.sendCommand("button_function 5 'accept_target_fixation'")

        # Shrink the spread of the calibration/validation targets
        # if the default outermost targets are not all visible in the bore.
        # The default <x, y display proportion> is 0.88, 0.83 (88% of the display
        # horizontally and 83% vertically)
        el_tracker.sendCommand('calibration_area_proportion 0.88 0.83')
        el_tracker.sendCommand('validation_area_proportion 0.88 0.83')

    # Get the native screen resolution used by PsychoPy
    scn_width, scn_height = win.size

    # Resolution fix for Mac retina displays
    if 'Darwin' in platform.system():
        if config.use_retina:
            scn_width = int(scn_width/2.0)
            scn_height = int(scn_height/2.0)
            
    # Optional: online drift correction.
    # See the EyeLink 1000 / EyeLink 1000 Plus User Manual

    # Online drift correction to mouse-click position:
    # el_tracker.sendCommand('driftcorrect_cr_disable = OFF')
    # el_tracker.sendCommand('normal_click_dcorr = ON')

    if not dummy_mode:
        
        # Online drift correction to a fixed location, e.g., screen center
        el_tracker.sendCommand('driftcorrect_cr_disable = OFF')
        el_tracker.sendCommand('online_dcorr_refposn %d,%d' % (int(scn_width/2.0),
                                                                int(scn_height/2.0)))
        el_tracker.sendCommand('online_dcorr_button = ON')
        el_tracker.sendCommand('normal_click_dcorr = OFF')

        # Pass the display pixel coordinates (left, top, right, bottom) to the tracker
        # see the EyeLink Installation Guide, "Customizing Screen Settings"
        el_coords = "screen_pixel_coords = 0 0 %d %d" % (scn_width - 1, scn_height - 1)
        el_tracker.sendCommand(el_coords)

        # Write a DISPLAY_COORDS message to the EDF file
        # Data Viewer needs this piece of info for proper visualization, see Data
        # Viewer User Manual, "Protocol for EyeLink Data to Viewer Integration"
        dv_coords = "DISPLAY_COORDS  0 0 %d %d" % (scn_width - 1, scn_height - 1)
        el_tracker.sendMessage(dv_coords)

    # Configure a graphics environment (genv) for tracker calibration
    genv = EyeLinkCoreGraphicsPsychoPy(el_tracker, win)
    print(genv)  # print out the version number of the CoreGraphics library

    # Set background and foreground colors for the calibration target
    # in PsychoPy, (-1, -1, -1)=black, (1, 1, 1)=white, (0, 0, 0)=mid-gray
    foreground_color = (-1, -1, -1)
    background_color = win.color # Use the same background color as the entire study
    genv.setCalibrationColors(foreground_color, background_color)

    # Set up the calibration target

    # Use a picture as the calibration target
    genv.setTargetType('circle')
    genv.setTargetSize(24)
    #genv.setPictureTarget(os.path.join('images', 'fixTarget.bmp')) #CALIBRATION TARGET IMAGE

    # Configure the size of the calibration target (in pixels)
    # this option applies only to "circle" and "spiral" targets
    # genv.setTargetSize(24)

    # Beeps to play during calibration, validation and drift correction
    # parameters: target, good, error
    #     target -- sound to play when target moves
    #     good -- sound to play on successful operation
    #     error -- sound to play on failure or interruption
    # Each parameter could be ''--default sound, 'off'--no sound, or a wav file
    genv.setCalibrationSounds('off', 'off', 'off')

    # Resolution fix for macOS retina display issues
    if config.use_retina:
        genv.fixMacRetinaDisplay()

    #if not dummy_mode:
    # Request Pylink to use the PsychoPy window we opened above for calibration
    pylink.openGraphicsEx(genv)

def calibrate_eyelink(win, dummy_mode): 
    # Note: The script below is provided by SR Research, Inc.
    el_tracker = pylink.getEYELINK()
    if not dummy_mode: 
         # Gaze calibration
        task_msg = 'Press O to calibrate tracker'
        instructions_screens(win, task_msg)
        try:
            el_tracker.doTrackerSetup()
            
        # Error
        except RuntimeError as err:
            print('ERROR:', err)
            el_tracker.exitCalibration()
        
        el_tracker.setOfflineMode()  
        # Try start recording
        try:
            el_tracker.startRecording(1, 1, 1, 1)
        
        # Error
        except RuntimeError as error:
            print("ERROR:", error)
            terminate_task()
        eye_used = el_tracker.eyeAvailable()
        
        if eye_used == 1:
            el_tracker.sendMessage("EYE_USED 1 RIGHT")
        elif eye_used == 0 or eye_used == 2:
            el_tracker.sendMessage("EYE_USED 0 LEFT")
            eye_used = 0
        else:
            print("Error in getting the eye information!")
    pylink.pumpDelay(100)
