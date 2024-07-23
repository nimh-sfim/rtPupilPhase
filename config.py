import os

## Real-time parameters

resolution = [1920, 1080]
use_retina = True     # Set this variable to True if you use the built-in retina screen as your 
    # primary display device on macOS. If have an external monitor, set this 
    # variable True if you choose to "Optimize for Built-in Retina Display" 
    # in the Displays preference settings.

bg_color = (116,116,116)
text_color = "black"
ms_per_sample = 17 # assumes online sampling rate of 60Hz
peak_threshold = 0.
trough_threshold = 0.
constriction_threshold = 50.
dilation_threshold = -50.

## Simulation Parameters 

### Eyelink Information
recorded_eye = 1  # Note: EyeLink stores the pupil size data in a 2 x time/sample matrix. The first row = the left eye and the second row = the right eye.
downsample_value = 17 # Downsample value - to make offline recording match live-stream recording rate
live_sampling_rate = 1000 # pupillometry offline sampling rate in Hz

### Filenames 
data_fname = os.path.join("data", "human")
results_fname = os.path.join("analysis", "subject_analysis","human")

### Task Structure 

num_blocks = 5
block_duration_ms = 600000 # duration of in milliseconds

# *** Display Parameters ***
# These may be necessary for calculating blinks, saccades and microsaccades

# Display monitor information
monitor_height = 1024 # in pixels
monitor_width = 768 # in pixels
pixel_pitch = 0.254 # number retrieved from the monitor manufacturer's website - .254mm per pixel