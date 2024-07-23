import os

### Real-time parameters

## Initial Thresholds
# These are the initial thresholds for the beginning of the alghorithm. 
# They reflect the values used in Kroenemer et al., 2024. 
peak_threshold = 0.
trough_threshold = 0.
constriction_threshold = 50.
dilation_threshold = -50.

### Eyelink Information
recorded_eye = 1  # Note: EyeLink stores the pupil size data in a 2 x time/sample matrix. The first row = the left eye and the second row = the right eye.
ms_per_sample = 17 # assumes online sampling rate of 60Hz
downsample_value = 17 # Downsample value - to make offline recording match live-stream recording rate 
# CW: is ms_per_sample always the same as downsample rate? 
live_sampling_rate = 1000 # pupillometry offline sampling rate in Hz

### Display monitor information ###  
# This is information about the monitor that the task is run on.
# These may be necessary for calculating blinks, saccades and microsaccades
resolution = [1920, 1080]
use_retina = True     # Set this variable to True if you use the built-in retina screen as your 
    # primary display device on macOS. If have an external monitor, set this 
    # variable True if you choose to "Optimize for Built-in Retina Display" 
    # in the Displays preference settings.
monitor_height = 1024 # in pixels
monitor_width = 768 # in pixels
pixel_pitch = 0.254 # number retrieved from the monitor manufacturer's website - .254mm per pixel

### Psychopy information ### 
# Information needed for PsychoPy to create the experiment. 
bg_color = (116,116,116) # background color of the window
text_color = "black" # color of the text

### Filenames ### 
# These are folders specifying where to find the data and where to put the outputs of the simulations. 
data_fname = os.path.join("data", "human") # where to find previously collected data
results_fname = os.path.join("analysis", "subject_analysis","human") # where to store subject level analysis for the simulations

### Task Structure 
## These reflect the fixation task reflected in Kroenemer et al., 2024. 
num_blocks = 5 # number of blocks of task. 
block_duration_ms = 600000 # duration of block in milliseconds