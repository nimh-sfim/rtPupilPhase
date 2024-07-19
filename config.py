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

recorded_eye = 1  
num_blocks = 5
downsample_value = 17 # Downsample value - to make offline recording match live-stream recording rate
live_sampling_rate = 1000 # pupillometry offline sampling rate in Hz

# *** Display Parameters ***
# These may be necessary for calculating blinks, saccades and microsaccades

# Display monitor information
monitor_height = 1024 # in pixels
monitor_width = 768 # in pixels
pixel_pitch = 0.254 # number retrieved from the monitor manufacturer's website - .254mm per pixel