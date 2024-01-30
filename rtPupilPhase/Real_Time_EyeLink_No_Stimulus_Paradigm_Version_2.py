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

# EyeLink Functions
from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy
from string import ascii_letters, digits

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
results_folder = 'EyeLink_Data'

# Define task name
task_name = 'Fixation Task'

# Define EyeLink folder - directory of task + folder name
eyelinkFolder = (results_folder + os.path.sep)

# Check for the behavioral data directory, otherwise make it
if not os.path.isdir(behavioral_folder):
        os.makedirs(behavioral_folder)  # If this fails (e.g. permissions) we will get error

# Check for the EyeLink data directory, otherwise make it
if not os.path.isdir(results_folder):
        os.makedirs(results_folder) # If this fails (e.g. permissions) we will get error

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

# ***********************
# *** TASK PARAMETERS ***
# ***********************

# Block duration 
block_duration_sec = 600 # in seconds

# Maximum number of blocks (before quitting task)
max_num_blocks = 10 

# *****************************************
# *** REAL TIME PUPILLOMETRY PARAMETERS ***
# *****************************************

# Dictionary:
# ms = milliseconds

# Number of milliseconds between EyeLink real time samples (Note: EyeLink has 
# a real time/live stream sampling rate of 60 Hz; offline sampling rate is up
# 2000 Hz)
ms_per_sample = 17

# Peak and trough pupil size quantiles
peak_pupil_quantile = 0.75
trough_pupil_quantile = 0.25

# Dilation and constriction pupil phase quantiles
dilation_quantile = 0.99
constriction_quantile = 0.01

# Pupil sample parameters
pupil_sample_duration_ms = 100 # in milliseconds
pupil_sample_duration_samples = pupil_sample_duration_ms // ms_per_sample

# Search window parameters
max_search_window_duration_ms = 5000 # in milliseconds
max_search_window_duration_samples = max_search_window_duration_ms // ms_per_sample

# Baseline window - used to reset pupil size and pupil size derivative thresholds
baseline_duration_ms = 5000 # in milliseconds
baseline_duration_samples = round(baseline_duration_ms/ms_per_sample) # in samples

# Interstimulus interval: IEI
IEI_duration_sec = 3 # in seconds
IEI_duration_samples = IEI_duration_sec // ms_per_sample # in samples

# Pupil phase-independent or "random" event parameters
num_random_event = 20 # number of random stimuli per task block
random_event_time_sec = block_duration_sec//num_random_event # in seconds

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
 
# ******************************
# *** Task General Functions ***
# ******************************

def clear_screen(win):
    """ Clear window """ 
    
    win.fillColor = genv.getBackgroundColor()
    win.flip()

def terminate_task():
    """ Terminate the task gracefully and retrieve the EDF data file """
    
    el_tracker = pylink.getEYELINK()
    
    if el_tracker.isConnected():
        error = el_tracker.isRecording()
    
        if error == pylink.TRIAL_OK:
            # Stop EyeLink recording
            abort_trial()
    
        el_tracker.setOfflineMode()
        el_tracker.sendCommand('clear_screen 0')
        pylink.msecDelay(500)
        el_tracker.closeDataFile()         
        el_tracker.sendMessage('End EyeLink Recording')
        el_tracker.close()
    
    win.close()
    core.quit()
    sys.exit()
    
def abort_trial():   
    """Ends recording abruptly"""
    
    el_tracker = pylink.getEYELINK()
    
    # Stop EyeLink recording
    if el_tracker.isRecording():
        pylink.pumpDelay(100)
        el_tracker.stopRecording()  
    
    clear_screen(win)
    bgcolor_RGB = (116, 116, 116)
    el_tracker.sendMessage('!V CLEAR %d %d %d' % bgcolor_RGB)
    el_tracker.sendMessage('TRIAL_RESULT %d' % pylink.TRIAL_ERROR)
    
    return pylink.TRIAL_ERROR

def end_experiment() -> None:
    """End experiment"""
    
    # Log
    logging.log(level=logging.EXP,msg='*** END EXPERIMENT ***')
    el_tracker.sendMessage('*** END EXPERIMENT ***')
    
    # Stop EyeLink recording
    pylink.pumpDelay(100)
    el_tracker.stopRecording()
    
    # Terminate task
    terminate_task()
    core.wait(2)
    win.close()
   
def quit_task() -> None:
    """Quit task based off of key presses"""
    
    # Record key presses
    allKeys = event.getKeys(['p','escape'])
    
    # Check keypresses
    if allKeys != None:
        
        # Loop over key presses
        for thisKey in allKeys:
            
            # If target key found, end experiment
            if thisKey in ['p', 'escape']:
                end_experiment()

def block_trigger():
    """Display block trigger screen"""
    
    # Log
    logging.log(level=logging.EXP,msg='Waiting for block trigger')
    el_tracker.sendMessage('Waiting for block trigger')
        
    # On-screen text
    onscreen_instructions = visual.TextStim(win, text='Waiting to start. Standby...', color = genv.getForegroundColor(), wrapWidth = scn_width/2) 
    onscreen_instructions.draw()
    win.flip()

    # Wait for trigger key press
    key = event.waitKeys(keyList=['5','t','escape', 'p'])

    # If pressed escape key, quit task
    if np.in1d(key,['escape','p']):
       end_experiment()
       
    # Log
    logging.log(level=logging.EXP,msg='Block trigger received')
    el_tracker.sendMessage('Block trigger received')

def instruction_continue():
    """Continue/proceed from instruction screen"""

    # Wait for key press
    key = event.waitKeys(keyList=['space', 'escape', 'p'])

    # End experiment
    if np.in1d(key,['escape','p']):
        end_experiment()
    
def instructions_screens(instruction:str): 
    """Function presents all the instructions needed for the task"""
    
    # Setup text
    task_instructions = visual.TextStim(win, instruction, color = genv.getForegroundColor(), wrapWidth = scn_width/2, units = 'cm')
    
    # Clear window
    clear_screen(win)
    
    # Draw text screen
    task_instructions.draw()
    win.flip()   
    
    # Continue or quit from instructions
    instruction_continue()

def general_instruction_screens():
    """Present the task instructions"""
    
    # Log
    logging.log(level=logging.EXP,msg='General task instructions')

    # Define instructions text
    instructions = "Please fixate on the fixation point at all times."
    
    # Setup instructions
    task_instructions = visual.TextStim(win, text='', color='black', pos=[0, 5], units='cm')

    # Draw the text
    task_instructions.setText(instructions)
    task_instructions.draw()
    fixation.draw()
    win.flip()

    # Continue or quit from instructions
    instruction_continue()
    
# **********************************************
# *** Define rtPupilPhase Class Functions ***
# **********************************************

class rtPupilPhase():
    """Object to store data and detect pupil phase events in real time."""
    
    def __init__(
        self,task_name
    ):
        self._prior_search_window = []
        self._search_window_model_fits = []
        self._new_sample = None
        self._old_sample = None
        self._search_window = []
        self._baseline_window = []
        self._search_window_sample_times = []
        self._pupil_sample_duration_time = []
        self._peak_count = 0
        self._trough_count = 0
        self._dilation_count = 0
        self._constriction_count = 0
        self._accepted_pupil_event_bool = False
        self._peak_threshold_var = 0 
        self._trough_threshold_var = 0 
        self._dilation_threshold = 50
        self._constriction_threshold = -50
        self._idx_event = 0

    def accepted_pupil_event(self):
        """
        Log an accepted pupil event (i.e., event detected beyond the inter-event interval
        Note: This function can be used for building closed-loop paradigms where
        detect pupil phase events trigger a task event.

        PARAMETERS
            self
            Interacting with the globals:
            logging
            el_tracker
            pupil_phase_IEI_time 
        OUTPUTS
            Changes self._accepted_pupil_event_bool to True
            Resets pupil_phase_IEI_time
            Logs to messages
        """

        # Log
        logging.log(level=logging.EXP,msg='Accepted Pupil Event')
        el_tracker.sendMessage('Accepted Pupil Event')
        
        # Reset pupil phase IEI timer
        pupil_phase_IEI_timer.reset()
        
        # Set boolean to True
        self._accepted_pupil_event_bool = True
    
    def detect_events(self) -> int:
        """Check for possible pupil phase event/reset search window."""
        
        # If search window is empty
        if len(self._search_window) == 0:
            self._accepted_pupil_event_bool = False
            return 0
        
        # If search window is too long - reset it
        if len(self._search_window) > max_search_window_duration_samples:
            self._accepted_pupil_event_bool = False
            self._prior_search_window = self._search_window
            self._search_window = []
            self._search_window_sample_times = []
            self._search_window_model_fits = []
            return 0
        
        # Reset search window if a NaN sample is detected (e.g., blink event)
        if any(np.isnan(self._search_window)):
            self._accepted_pupil_event_bool = False
            self._prior_search_window = self._search_window
            self._search_window = []
            self._search_window_sample_times = []
            self._search_window_model_fits = []
            return 0
        
        # Reset search window if a NaN sample was found in a prior search window
        if any(np.isnan(self._prior_search_window)):
            self._accepted_pupil_event_bool = False
            self._prior_search_window = self._search_window
            self._search_window = []
            self._search_window_sample_times = []
            self._search_window_model_fits = []
            return 0
            
        # Demean search window        
        demeaned_search_window = demean_search_window(self._search_window)
        
        # Reset pupil phase pupil size and pupil size derivative thresholds - once the minimum baseline duration has been acquired
        if len(self._baseline_window) > baseline_duration_samples:
            
            # Input baseline window data to threshold updating function
            self.set_pupil_phase_thresholds(self._baseline_window)
        
        # Log random event if the minimum random IEI is exceed      
        if random_IEI_timer.getTime() > random_event_time_sec:
        
            # If IEI has been exceed (Note: The IEI timer gets reset in the accepted_pupil_event; 
            # no reset will happen if no stimulus is shown)
            if pupil_phase_IEI_timer.getTime() > IEI_duration_sec:
                
                # Define idx extrema
                self._idx_event = 3
                
                # Log
                logging.log(level=logging.EXP,msg='Random Event')
                el_tracker.sendMessage('Random Event') 
                
                # Accepted event
                self.accepted_pupil_event()

                # Reset timers
                random_IEI_timer.reset()
                general_timer.reset()
                
                return 0
            
            return 0
        
        # Find peaks (local maxima), troughs (local minima), dilation (dilation), and constriction (constriction) events
        # Index extrema dictionary: peak = 1; dilation = 2; trough = -1; constriction = -2; random = 3
        self._idx_event = self.find_pupil_phase_event(demeaned_search_window)
        
        # If a pupil phase event was found
        if self._idx_event!=0:
            
            # Confirm pupil phase IEI exceeded
            if pupil_phase_IEI_timer.getTime() > IEI_duration_sec:
                
                # Accepted event
                self.accepted_pupil_event()
                
                # Reset search window
                self._prior_search_window = self._search_window
                self._search_window = []
                self._search_window_sample_times = []
                self._search_window_model_fits = []
                self._idx_event = 0
            
            # If IEI is not exceeded 
            elif pupil_phase_IEI_timer.getTime() < IEI_duration_sec:
                
                # Not an accepted event
                self._accepted_pupil_event_bool = False
                
                # Log
                logging.log(level=logging.EXP,msg='Within IEI - Skipping this Pupil Event')
                el_tracker.sendMessage('Within IEI - Skipping this Pupil Event')
            
            # NOTE: search_window is not reset; will continue looking for a pupil phase events after adding another pupil_sample;
            # unless the search window gets too long
            
            # Not an accepted event
            self._accepted_pupil_event_bool = False
            return 0
        
        # Not an accepted event    
        self._accepted_pupil_event_bool = False
        return 0
    
    def get_pupil_sample_duration_time(self):
        """Used to log pupil sample duration array."""
        return self._pupil_sample_duration_time
    
    def set_pupil_phase_thresholds(self, baseline_window):
        """Reset pupil phase event thresholds according to current baseline window"""
        
        # Demean baseline window
        demeaned_baseline_window = list(baseline_window - np.nanmean(baseline_window))
        
        # Log
        logging.log(level=logging.EXP,msg='Demeaned Baseline Window: ' + str(demeaned_baseline_window))

        # Get rid of NaN samples
        np_demeaned_baseline_window = np.array(demeaned_baseline_window)
        demeaned_baseline_window = np_demeaned_baseline_window[~np.isnan(np_demeaned_baseline_window)]
        
        # Find peaks in baseline window
        retro_peaks_idx, _ = find_peaks(demeaned_baseline_window)
        
        # Store all peak pupil sizes
        retro_peaks = []
        for i in retro_peaks_idx:
            k = demeaned_baseline_window[i]
            retro_peaks.append(k)
            
        # Find troughs in retrospect
        # Note: Pupil data is inverted before searching for troughs
        retro_troughs_idx, _ = find_peaks(np.negative(demeaned_baseline_window))
        
        # Store all trough pupil sizes
        retro_troughs = []
        for i in retro_troughs_idx:
            k = demeaned_baseline_window[i]
            retro_troughs.append(k)

        # Find diff/gradient of pupil size in baseline window
        diff_demeaned_baseline_window = np.diff(demeaned_baseline_window)

        # Find peak and trough threshold quantiles
        retro_peaks_quantile = np.quantile(retro_peaks, peak_pupil_quantile)
        retro_troughs_quantile = np.quantile(retro_troughs, trough_pupil_quantile)
        
        # Find dilation and constriction threshold quantiles
        retro_dilation_quantile = np.quantile(diff_demeaned_baseline_window, dilation_quantile)
        retro_constriction_quantile = np.quantile(diff_demeaned_baseline_window, constriction_quantile)
        
        # Set new thresholds
        self._peak_threshold_var = retro_peaks_quantile
        self._trough_threshold_var = retro_troughs_quantile
        self._dilation_threshold = retro_dilation_quantile
        self._constriction_threshold = retro_constriction_quantile
        
        # Reset baseline window
        self._baseline_window = []
        
        # Log
        logging.log(level=logging.EXP,msg='Updated peak threshold: ' + str(self._peak_threshold_var))
        logging.log(level=logging.EXP,msg='Updated trough threshold: ' + str(self._trough_threshold_var))
        logging.log(level=logging.EXP,msg='Updated dilation threshold: ' + str(self._dilation_threshold))
        logging.log(level=logging.EXP,msg='Updated constriction threshold: ' + str(self._constriction_threshold))
        
    def find_pupil_phase_event(self, demeaned_search_window: np.array) -> int:
        """Find possible pupil phases"""
                
        # Fit polynomial
        sample_window_fit_coef = np.polyfit(list(range(len(demeaned_search_window))),demeaned_search_window, 2)
 
        # Find last sample fit value
        fit_value = np.polyval(sample_window_fit_coef,len(demeaned_search_window))

        # Store the last value of fitting curve
        self._search_window_model_fits = np.append(self._search_window_model_fits,fit_value)
        
        # Diff across fit values if there are more than one fit values
        if len(self._search_window_model_fits) > 1:
            
            # Find diff of fit values
            diff_fit = list(np.diff(self._search_window_model_fits))
            
            # *** Find trough event ***
            if diff_fit[-1] > 0 and demeaned_search_window[-1] < self._trough_threshold_var:
                
                # Log
                logging.log(level=logging.EXP,msg='Search Window Pupil: ' + str(demeaned_search_window))
                logging.log(level=logging.EXP,msg='Search Window Sample Time: ' + str(self._search_window_sample_times))
                logging.log(level=logging.EXP,msg='Search Window Fit Diff: ' + str(diff_fit))
                logging.log(level=logging.EXP,msg='Found Trough')
                el_tracker.sendMessage('Found Trough')
                
                # Add to event counter
                self._trough_count = self._trough_count + 1

                return -1

            # *** Find peak event ***
            elif diff_fit[-1] < 0 and demeaned_search_window[-1] > self._peak_threshold_var:
                
                # Log
                logging.log(level=logging.EXP,msg='Search Window Pupil: ' + str(demeaned_search_window))
                logging.log(level=logging.EXP,msg='Search Window Sample Time: ' + str(self._search_window_sample_times))
                logging.log(level=logging.EXP,msg='Search Window Fit Diff: ' + str(diff_fit))
                logging.log(level=logging.EXP,msg='Found Peak')
                el_tracker.sendMessage('Found Peak')
                
                # Add to event counter
                self._peak_count = self._peak_count + 1

                return 1
            
            # *** Find dilation event ***
            elif diff_fit[-1] > self._dilation_threshold:
                
                # Log
                logging.log(level=logging.EXP,msg='Search Window Pupil: ' + str(demeaned_search_window))
                logging.log(level=logging.EXP,msg='Search Window Sample Time: ' + str(self._search_window_sample_times))
                logging.log(level=logging.EXP,msg='Search Window Fit Diff: ' + str(diff_fit))
                logging.log(level=logging.EXP,msg='Found Dilation Event')
                el_tracker.sendMessage('Found Dilation Event')
                
                # Add to event counter
                self._dilation_count = self._dilation_count + 1

                return 2
            
            # *** Find constriction event ***
            elif diff_fit[-1] < self._constriction_threshold:
                
                # Log
                logging.log(level=logging.EXP,msg='Search Window Pupil: ' + str(demeaned_search_window))
                logging.log(level=logging.EXP,msg='Search Window Sample Time: ' + str(self._search_window_sample_times))
                logging.log(level=logging.EXP,msg='Search Window Fit Diff: ' + str(diff_fit))
                logging.log(level=logging.EXP,msg='Found Constriction Event')
                el_tracker.sendMessage('Found Constriction Event')
                
                # Add to event counter
                self._constriction_count = self._constriction_count + 1

                return -2
            
            # No extrema found
            else: 
                self._accepted_pupil_event_bool = False
                return 0
                
        # Not enough fit values to complete diff analysis
        else:
            self._accepted_pupil_event_bool = False
            return 0
            
    def build_search_window(self) -> float:
        """Add the pupil sample to search window"""
        
        # Initialize/reset pupil_sample variable
        pupil_sample = []
        pupil_sample_time = []
        
        # Setup pupil size/time variables
        p_size = float("nan")
        p_time = float("nan")
        
        # Reset pupil sample timer
        pupil_sample_timer.reset()
        
        # Keep adding pupil values until pupil_sample is long enough      
        while len(pupil_sample) < pupil_sample_duration_samples:
        
            # Update window
            win.update()
        
            # Get the latest EyeLink sample        
            self._new_sample = el_tracker.getNewestSample()
    
            # Check there is a new pupil value
            if self._new_sample is not None:
                
                # Check there is an old pupil value (i.e., not the first pupil sample)
                if self._old_sample is not None:
                    
                    # Check the new and old pupil values are not the same
                    if self._new_sample.getTime() != self._old_sample.getTime():
    
                        # Find the pupil size and time
                        p_size = self._new_sample.getRightEye().getPupilSize()
                        p_time = self._new_sample.getTime()
                        
                        # Replace pupil size of 0 with "nan"
                        if p_size == 0:
                            p_size = float("nan")
                        
                        # Add to pupil_sample
                        pupil_sample.append(p_size) # Add samples to pupil_sample
                        pupil_sample_time.append(p_time)
            
            # Quit task
            quit_task()
        
            # Replace old pupil value with new
            self._old_sample = self._new_sample
            
        # Add pupil sample and pupil sample times to search window and search window sample times 
        self._search_window.extend(pupil_sample)
        self._search_window_sample_times.extend(pupil_sample_time)
        
        # Store duration of pupil sample in time
        self._pupil_sample_duration_time.append(pupil_sample_timer.getTime())
        
        # Add pupil sample to baseline window
        self._baseline_window.extend(pupil_sample)

def demean_search_window(search_window: np.array) -> np.ndarray:
    """Calculate the mean of the pupil_sample and subtract from all data samples."""

    return list(search_window - np.mean(search_window))

# *********************
# *** MAIN FUNCTION ***
# *********************
def main():
    
    # Define globals
    global input_x
    global input_fraction
    global input_k_
    global input_lambda_
    global task_name
    
    # Setup classes
    sd = rtPupilPhase(task_name) # sets up stimulus decider object
    
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
        while block_timer.getTime() < block_duration_sec:
            
            # Display a halfway completion screen
            if halfway_screen == False and block_timer.getTime() > block_duration_sec/2:

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
            decision_arr.append(sd.detect_events())
        
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
    main()