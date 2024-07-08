# ********************************************
# *** REAL TIME PUPILLOMETRY FIXATION TASK ***
# ********************************************

# This script can be run in PsychoPy to present a fixation task embedded with 
# the rtPupilPhase real time monitoring of pupil size fluctuations method.

# Written by: Sharif I. Kronemer and Tori Gobo
# Last Modified: 1/2/2024

# Version #
task_version = 'v2'

# ************************
# *** IMPORT LIBRARIES ***
# ************************

# Load Python and Psychopy Functions
from psychopy import visual, gui, data, core, event, logging
import numpy as np
import time
import pylink
import os
import platform
import sys
from scipy.signal import find_peaks
import argparse

# EyeLink Functions
#from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy
from string import ascii_letters, digits

from PsychoPy_funcs import *
from StimulusDecider import StimulusDecider

# *********************
# *** MAIN FUNCTION ***
# *********************

def main(task_name, behavioral_folder, eyelink_folder, block_length, max_num_blocks, baseline_duration_ms,
         max_search_window_duration_ms, num_random_events, IEI_duration_sec, ms_per_sample, 
         pupil_sample_duration_ms, peak_pupil_quantile, trough_pupil_quantile, 
         dilation_quantile, constriction_quantile, peak_threshold, trough_threshold,
         constriction_threshold, dilation_threshold):
    
    # Define globals
    global input_x
    global input_fraction
    global input_k_
    global input_lambda_


    # ***************************
    # ******* USER INPUTS *******
    # ***************************

    # Setup the subject info screen
    info = {'Session #': 1, 'Subject ID': 'Test', 'EyeLink': ['y','n'], 'EyeLink EDF': 'test.edf', '(1) Skip task instructions':['n', 'y']}

    dlg = gui.DlgFromDict(info, title = 'Real Time Pupillometry Fixation Experiment')

    # Find experiment date
    info['date'] = data.getDateStr()

    # Filename = Subject ID entered above
    filename = info['Subject ID']

    # EyeLink EDF filename
    tmp_str = info['EyeLink EDF']

    # ********************
    # *** SETUP WINDOW ***
    # ********************

    # Screen resolution in pixels
    resolution = [1920, 1080]
        
    # Setup the display window
    win = visual.Window(size = resolution, color = [0,0,0], monitor = 'testMonitor', fullscr = True, units ='cm')

    # *****************************************
    # *** MANAGE DATA FOLDERS AND FILENAMES ***
    # *****************************************

    # Behavioral data file
    behavioral_folder = 'Behavioral_Data'

    # Eyelink data file
    eyelink_folder = 'EyeLink_Data'

    # Define task name
    task_name = 'Fixation Task'

    # Define EyeLink folder - directory of task + folder name
    eyelinkFolder = (eyelink_folder + os.path.sep)

    # Check for the behavioral data directory, otherwise make it
    if not os.path.isdir(behavioral_folder):
            os.makedirs(behavioral_folder)  # If this fails (e.g. permissions) we will get error

    # Check for the EyeLink data directory, otherwise make it
    if not os.path.isdir(eyelink_folder):
            os.makedirs(eyelink_folder) # If this fails (e.g. permissions) we will get error

    # Show only critical log messages in the PsychoPy console
    logFile = logging.LogFile(behavioral_folder + os.path.sep + filename + '_Session_'+str(info['Session #'])+'_Real_Time_Pupillometry_Fixation_'+info['date']+'_'+task_version+'.log', level=logging.EXP)

    # Specify code directory 
    _thisDir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(_thisDir)

    # Note: The script below is provided by SR Research, Inc.

    # Strip trailing characters, ignore the ".edf" extension
    edf_fname = tmp_str.rstrip().split('.')[0]

    # Check if the EyeLink filename is valid (length <= 8 & no special char)
    allowed_char = ascii_letters + digits + '_'

    # If too many characters in EyeLink filename
    if not all([c in allowed_char for c in edf_fname]):
        
        print('ERROR: *** Invalid EDF filename')
        core.quit()  # abort experiment

    elif len(edf_fname) > 8:
        
        print('ERROR: *** EDF filename should not exceed 8 characters')
        core.quit()  # abort experiment

    # We download EDF data file from the EyeLink Host PC to the local hard
    # drive at the end of each testing session, here we rename the EDF to
    # include session start date/time
    time_str = time.strftime("_%Y_%m_%d_%H_%M", time.localtime())
    session_identifier = edf_fname + time_str

    # ********************
    # *** TASK STIMULI ***
    # ********************

    # Setup fixation cross
    fixation = visual.TextStim(win, text="+", color = 'black', pos = [0, 0], autoLog = False)

    # *************************
    # *** TASK INSTRUCTIONS ***
    # *************************

    # Create text screens to display later:
    instructions = visual.TextStim(win, text='', color='white', pos=[0, 2])  #This is an empty instructions screen to be filled with text below

    # *******************
    # *** TASK TIMERS ***
    # *******************

    block_timer = core.Clock() # Block timer
    general_timer = core.Clock() # Global timer
    pupil_phase_IEI_timer = core.Clock() # Inter-event interval timer
    random_IEI_timer = core.Clock() # Time for selecting random pupil times
    pupil_sample_timer = core.Clock() # Pupil sample timer

    # ************************
    # *** INITIATE EYELINK ***
    # ************************

    # Note: The script below is provided by SR Research, Inc.

    # EyeLink Dummy mode? - Set to False if testing with actual system
    if info['EyeLink'] == 'y':
        dummy_mode = False
        
    elif info['EyeLink'] == 'n':
        dummy_mode = True

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
        eyelink_ver = 5
        
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
        
    el_tracker.sendCommand("file_event_filter = %s" % file_event_flags)
    el_tracker.sendCommand("file_sample_data = %s" % file_sample_flags)
    el_tracker.sendCommand("link_event_filter = %s" % link_event_flags)
    el_tracker.sendCommand("link_sample_data = %s" % link_sample_flags)

    # Set EyeLink sample rate
    if eyelink_ver > 2 and not dummy_mode:
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

    # Set this variable to True if you use the built-in retina screen as your 
    # primary display device on macOS. If have an external monitor, set this 
    # variable True if you choose to "Optimize for Built-in Retina Display" 
    # in the Displays preference settings.
    use_retina = True

    # Resolution fix for Mac retina displays
    if 'Darwin' in platform.system():
        if use_retina:
            scn_width = int(scn_width/2.0)
            scn_height = int(scn_height/2.0)
            
    # Optional: online drift correction.
    # See the EyeLink 1000 / EyeLink 1000 Plus User Manual

    # Online drift correction to mouse-click position:
    # el_tracker.sendCommand('driftcorrect_cr_disable = OFF')
    # el_tracker.sendCommand('normal_click_dcorr = ON')

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
    if use_retina:
        genv.fixMacRetinaDisplay()

    # Request Pylink to use the PsychoPy window we opened above for calibration
    pylink.openGraphicsEx(genv)

    # Calibration task constants
    # Set random seed
    rng = np.random.default_rng()
    
    
    # Setup classes
    sd = StimulusDecider(task_name, ms_per_sample, block_length, baseline_duration_ms, 
                        max_search_window_duration_ms, pupil_sample_duration_ms, 
                        num_random_events, IEI_duration_sec, peak_pupil_quantile,
                        trough_pupil_quantile, dilation_quantile, constriction_quantile,
                        peak_threshold, trough_threshold, dilation_threshold, constriction_threshold,
                        online=True) # sets up stimulus decider object
    
    # *********************
    # *** Setup EyeLink ***
    # *********************
    
    # Note: The script below is provided by SR Research, Inc.
    
    # If running EyeLink
    if not dummy_mode:
        
        # Gaze calibration
        task_msg = 'Press O to calibrate tracker'
        instructions_screens(task_msg)
    
    # If running EyeLink
    if not dummy_mode:
        
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
    
    # Right eye 
    if eye_used == 1:
        el_tracker.sendMessage("EYE_USED 1 RIGHT")
    
    # Left eye
    elif eye_used == 0 or eye_used == 2:
        el_tracker.sendMessage("EYE_USED 0 LEFT")
        eye_used = 0
    
    else:
        print("Error in getting the eye information!")
    
    pylink.pumpDelay(100)
    
    # *******************************
    # *** Beginning of Experiment ***
    # *******************************
    
    # Define instruction text
    instructions_screens("Experiment is setup!\n\nLet's get started!")
    
    # Continue with task instructions
    if 'n' == info['(1) Skip task instructions']:

        # Instructions screen
        general_instruction_screens()
    
    # ************************
    # *** Main Task Phases ***
    # ************************
 
    # Setup block counter
    block_counter = 1
    
    # Log
    logging.log(level=logging.EXP,msg='Starting Main Task Phase')
    el_tracker.sendMessage("Starting Main Task Phase")
    
    # Loop over task blocks
    while block_counter <= max_num_blocks:
        
        # Count blocks
        block_instruction = 'Starting '+ task_name + ' Block #' + str(block_counter) + "\n\nAre you ready?"

        # Display instructions
        instructions_screens(block_instruction)
        
        # Block Trigger
        block_trigger()
        
        # Initialize variable
        decision_arr = []
        
        # Log
        logging.log(level=logging.EXP,msg='Starting ' + task_name + ' Block ' + str(block_counter))
        el_tracker.sendMessage('Starting ' + task_name + ' Block ' + str(block_counter))
        
        # Reset halfway logical 
        halfway_screen = False
        
        # Reset pupil phase IEI timer
        pupil_phase_IEI_timer.reset()
    
        # Reset block timer
        block_timer.reset()

        # Wait block duration
        while block_timer.getTime() < block_length:
            
            # Display a halfway completion screen
            if halfway_screen == False and block_timer.getTime() > block_length/2:

                # Reset
                halfway_screen = True

                # Setup progress screen
                progress_screen = visual.TextStim(win, text="You completed 50% of this block!", color='black')

                # Show progress screen/ turn off fixation
                fixation.setAutoDraw(False)
                progress_screen.setAutoDraw(True)

                # Reset timer
                general_timer.reset()
                
                # Present stimulus for target duration
                while general_timer.getTime() < 2:
                    win.update()

                # Turn on fixation
                progress_screen.setAutoDraw(False)
                fixation.setAutoDraw(True)
                win.update()
                
                # Reset timer
                general_timer.reset()
                
                # Timeout period before searching for pupil phase event
                while general_timer.getTime() < 2:
                    win.update()

            # Get all pressed keys
            # NOTE: allKeys will reset itself when getKeys is called
            allKeys = event.getKeys(['p','escape'])
    
            # If a key was pressed
            if allKeys != None:
    
                # Loop over key presses
                for thisKey in allKeys:
    
                    # End experiment
                    if np.in1d(thisKey,['escape','p']):
                        end_experiment()
            
            # Turn on fixation
            fixation.setAutoDraw(True)
            win.update()
        
            # Build search window
            sd.build_search_window()

            # Look for pupil phase events
            decision_arr.append(sd.detect_events_online())
        
        # Log
        logging.log(level=logging.EXP,msg='All pupil_sample Duration Times: ' + str(sd.get_pupil_sample_duration_time()))
        logging.log(level=logging.EXP,msg='All Search Window Detected Pupil Phase Events: ' + str(decision_arr))
        logging.log(level=logging.EXP,msg='Finished ' + task_name + ' Block ' + str(block_counter))
        el_tracker.sendMessage('Finished ' + task_name + ' Block ' + str(block_counter))
        
        # Turn off fixation point
        fixation.setAutoDraw(False)
        win.update()
        
        # Show instruction screen
        instructions_screens('Finished ' + task_name + ' Block #' + str(block_counter))
            
        # Add to block counter
        block_counter = block_counter + 1
        
    # *******************************
    # ***** End of Experiment *******
    # *******************************
    
    # Show instruction screen
    instructions_screens("Exiting task! \n\nAn experimenter will communicate with you shortly.")
    end_experiment()
            
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'rtPupilPhase: Real-Time Pupillometry')
    parser.add_argument("behavioral_folder",help="Directory where behavioral data should be stored")
    parser.add_argument("eyelink_folder", help="Directory where EyeLink data should be stored")
    parser.add_argument("task_name", help="Name of task")
    parser.add_argument("--max_num_blocks", help="Number of task blocks. Default: 10 blocks", 
                        type=int, default=10)
    parser.add_argument("--block_length", help="Duration of a single block, in seconds. Default: 600s", 
                        type=int, default=600)
    parser.add_argument("--baseline_duration_ms", help="Duration of baseline window in milliseconds. Default: 5000ms", 
                        type=int, default=5000)
    parser.add_argument("--max_search_window_duration_ms", help="Maximum duration of search window before resetting, in milliseconds. Default: 5000ms",
                         type=int, default=5000)
    parser.add_argument("--num_random_events", help="Number of random events per block. Default: 20 events", 
                        type=int, default=20)
    parser.add_argument("--IEI_duration_sec", help="Inter-event interval - how long to wait between valid events in seconds. Default: 3s", 
                        type=int, default=3)
    parser.add_argument("--ms_per_sample", help="Length of a single real-time pupil sample in milliseconds. Default: 17ms", 
                        type=int, default=17)
    parser.add_argument("--pupil_sample_duration_ms", help="How long we should consider a pupil sample in milliseconds. Default: 100ms", 
                        type=int, default=100)
    parser.add_argument("--peak_pupil_quantile", help="Quantile value a peak must be bigger than to accept. Default: 0.75",
                        type=float, default=0.75)
    parser.add_argument("--trough_pupil_quantile", help="Quantile value a trough must be smaller than to accept. Default: 0.25",
                        type=float, default=0.25)
    parser.add_argument("--dilation_quantile", help="Quantile value a dilation must be bigger than to accept. Default: 0.99",
                        type=float, default=0.99)
    parser.add_argument("--constriction_quantile", help="Quantile value a constriction must be smaller than to accept. Default: 0.01",
                        type=float, default=0.01)
    parser.add_argument("--peak_threshold", help="Initial threshold value that a pupil sample must be greater than to accept a peak. Default: 0.",
                    type=float, default=0.)
    parser.add_argument("--trough_threshold", help="Initial threshold value that a pupil sample must be smaller than to accept a trough. Default: 0.",
                    type=float, default=0.)
    parser.add_argument("--dilation_threshold", help="Initial threshold value that the change between pupil samples must be greater than to accept a dilation. Default: 50.",
                    type=float, default=50.)
    parser.add_argument("--constriction_threshold", help="Initial threshold value that the change between pupil samples must be smaller than to accept a constriction. Default: -50.",
                    type=float, default=-50.)
    
    args = parser.parse_args()


    main(args.task_name, args.behavioral_folder, args.eyelink_folder, args.block_length, 
        args.max_num_blocks,
        args.baseline_duration_ms, args.max_search_window_duration_ms, args.num_random_events, 
        args.IEI_duration_sec, args.ms_per_sample, args.pupil_sample_duration_ms, 
        args.peak_pupil_quantile, args.trough_pupil_quantile, args.dilation_quantile, 
        args.constriction_quantile, args.peak_threshold, args.trough_threshold, 
        args.constriction_threshold, args.dilation_threshold)