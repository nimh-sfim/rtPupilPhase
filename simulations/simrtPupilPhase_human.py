import numpy as np 
import os 
import pickle

from utils import get_data, find_string_time, pull_pupil_sample, plot_mean_timecourses
from StimulusDecider import StimulusDecider
from EventCollector import EventCollector

# Dictionary:
# ms = milliseconds
# IEI = inter-event interval

# *** Subject and Recording Parameters ***

# Subject list
subject_list = ['046','048','073','074','078','079','080','081']
# The recorded eye (0 = left; 1 = right)
# Note: EyeLink stores the pupil size data in a 2 x time/sample matrix. The
# first row = the left eye and the second row = the right eye.
recorded_eye = 1

# Number of task blocks completed per participant
num_blocks = 5
block_duration_ms = 600000 # duration of in milliseconds

# Define pupillometry sampling rate
sampling_rate = 1000 # in Hz
ms_per_sample = 17 # EyeLink live recording rate is (60Hz or ~17 samples/s)
downsample_value = 17 # Downsample value - to make offline recording match live-stream recording rate

# *** Display Parameters ***
# These may be necessary for calculating blinks, saccades and microsaccades

# Display monitor information
monitor_height = 1024 # in pixels
monitor_width = 768 # in pixels
pixel_pitch = 0.254 # number retrieved from the monitor manufacturer's website - .254mm per pixel

# *** rtPupilPhase Parameters ***

# Pupil sample parameters
pupil_sample_duration_ms = 100 # in milliseconds
samples_in_pupil_sample = np.round(pupil_sample_duration_ms/ms_per_sample)

# Random event parameters
# Note: The number of random events specified and the block duration will
# determine the random IEI. Also, note that the number of random events
# selected will have implications on the number of pupil phase events.
num_random_events = 15
random_IEI = block_duration_ms/num_random_events # in milliseconds

# Baseline window duration for setting new pupil size and derivative thresholds
baseline_window_ms = 5000; # in milliseconds
samples_in_baseline_window = np.round(baseline_window_ms/ms_per_sample)

# Inter-event interval for pupil phase events
IEI_jitter_ms = 3000 # in milliseconds
samples_in_IEI = np.round(IEI_jitter_ms/ms_per_sample) # in samples

# Maximum length of search window
max_search_window_length_ms = 5000 # in milliseconds

# *** Other Parameters ***

# 50# of the epoch length to be extracted - use 50% so that can plot both before and after event
half_epoch_duration_ms = 2500 # in milliseconds

# No blinks or saccades interval
no_blinks_saccades = 500 # in milliseconds

# define paths 
root_path = os.getcwd()
data_base = os.path.join(root_path,"simulations", "data","human")
results_base = os.path.join(root_path, "simulations", "analysis", "subject_analysis","human")

for subjID in subject_list:

    print('running subject: '+subjID +"...")

    # make subject specific directories
    sub_dir = os.path.join(data_base,subjID, "EyeLink")
    results_dir = os.path.join(results_base, subjID)

    # get data
    # time, pupil size, x coord, y coord
    gaze_data, event_data = get_data(sub_dir)

    # if output directory doesn't exist, create it 
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir) 

    # create stimulus decider object 
    sd = StimulusDecider("fixation")

    # identify events across block - this matters because we care about the ends of blocks. 
    # event epochs are concatenated across blocks at end of loop 
    for block in range(1,num_blocks+1): 
        print("Starting block "+str(block))

        # initialize Event Collectors 
        random_events = EventCollector("random")
        trough_events = EventCollector("trough")
        peak_events = EventCollector("peak")
        constriction_events = EventCollector("constriction")
        dilation_events = EventCollector("dilation")
        accepted_pupil_event_times = []
        
        # pull the raw data for a block 

        start_time = find_string_time(time_array = event_data['sttime'][0], message_array = event_data['message'][0], 
                                      match_string="Starting Perception Rate Block "+str(block))
        end_time = find_string_time(time_array = event_data['sttime'][0], message_array = event_data['message'][0], 
                                      match_string="Finished Perception Rate Block "+str(block))
        
        block_data = gaze_data[:,np.where(gaze_data[0]==start_time)[0][0]:np.where(gaze_data[0]==end_time)[0][0]]
        
        # Downsampling matches th Eyelink real-time sampling in a live testing session
        # goes from 1000Hz to 60Hz 
        downsampled_block_data = block_data[:, 0::downsample_value]

        # ensures we only look at complete pupil samples 
        for pupil_sample_num in range(int(np.shape(downsampled_block_data)[1]/6)-1):
            # stage 1: get pupil sample 
            current_pupil_sample = pull_pupil_sample(downsampled_block_data, pupil_sample_num, samples_in_pupil_sample)

            # stage 2: create search window 
            sd.update_windows(current_pupil_sample)

            # if the baseline window is long enough (and not all nans), update thresholds so we account for drift
            if len(sd.get_baseline_window()) > samples_in_baseline_window:
                if not np.isnan(sd.get_baseline_window()).all(): 
                    sd.set_pupil_phase_thresholds()
                else: 
                    sd.reset_baseline_window()

            # validate search window - not too long, no blinks
            if not sd.validate_search_window(max_search_window_length_ms/ms_per_sample):
                sd.reset_search_window()

                continue
            
            # demean search window 
            demeaned_search_window = sd.get_search_window() - np.mean(sd.get_search_window())

            # Fit demeaned data in search window with polynomial function 
            sd.fit_polynomial(demeaned_search_window)

            # stage 4: find pupil events 

            # check if random event happened 
            current_time = sd.get_search_window_times()[-1]

            if len(random_events.get_times()) >= 1: 
                time_from_last_event = current_time - random_events.get_times()[-1]
            else: 
                # ensures that the first stimulus always triggers accepted random event
                time_from_last_event = random_IEI

            # if it's been enough time, log a random event 
            if time_from_last_event >= random_IEI:
                random_events.update_data(((pupil_sample_num+1)*samples_in_pupil_sample), current_time, 
                                        sd.get_search_window()[-1])
            
            # find pupil phase events 
            if len(sd.get_search_window()) > 1:
                # check if there was an event 
                found_event = sd.find_pupil_phase_event(pupil_sample_num=pupil_sample_num, current_time=current_time, 
                                                        peak_events=peak_events, trough_events=trough_events, 
                                                        constriction_events=constriction_events, dilation_events=dilation_events)
                
                # update internal tracking of events 
                sd.set_current_event_idx(found_event)

                # get times for all events (regardless of type) 
                event_times = peak_events.get_times() + trough_events.get_times() +  constriction_events.get_times() + dilation_events.get_times()
                event_times.sort() # make sure they're in order

                # if there is an event, check to see whether enough time has passed to consider it an accepted event
                # events get logged in respective EventCollector objects 
                if found_event != 0:
                    peak_events, trough_events, dilation_events, constriction_events = sd.validate_event_offline(event_times, accepted_pupil_event_times, IEI_jitter_ms, 
                                            peak_events, trough_events, dilation_events, constriction_events)
                    
        
        # Note that it may be recommended to perform additional blink and microsaccade detection using your method of choice here. 
        # In the MATLAB simulations, we interpolate over blinks and use a method to identify microsaccades based off of Engbert & Kliegl (2003)     
        # Here, we do not perform this data processing and simply plot the data. 

        # pull epoch data for each event and concatentate across block 
        if block == 1:
            accepted_constriction_epoch_data, all_constriction_epoch_data = constriction_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, ms_per_sample, False)
            accepted_dilation_epoch_data, all_dilation_epoch_data = dilation_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, ms_per_sample, False)
            accepted_peak_epoch_data, all_peak_epoch_data = peak_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, ms_per_sample, False)
            accepted_trough_epoch_data, all_trough_epoch_data = trough_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, ms_per_sample, False)
            all_random_epoch_data, _ = random_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, ms_per_sample, True) 

        else: 
            block_accepted_constriction_epoch_data, block_all_constriction_epoch_data = constriction_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, ms_per_sample, False)
            accepted_constriction_epoch_data = np.append(accepted_constriction_epoch_data, block_accepted_constriction_epoch_data, axis=1)
            all_constriction_epoch_data = np.append(all_constriction_epoch_data, block_all_constriction_epoch_data, axis=1)

            block_accepted_dilation_epoch_data, block_all_dilation_epoch_data = dilation_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, ms_per_sample, False)
            accepted_dilation_epoch_data = np.append(accepted_dilation_epoch_data, block_accepted_dilation_epoch_data, axis=1)
            all_dilation_epoch_data = np.append(all_dilation_epoch_data, block_all_dilation_epoch_data, axis=1)

            block_accepted_peak_epoch_data, block_all_peak_epoch_data = peak_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, ms_per_sample, False)
            accepted_peak_epoch_data = np.append(accepted_peak_epoch_data, block_accepted_peak_epoch_data, axis=1)
            all_peak_epoch_data = np.append(all_peak_epoch_data, block_all_peak_epoch_data, axis=1)

            block_accepted_trough_epoch_data, block_all_trough_epoch_data = trough_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, ms_per_sample, False)
            accepted_trough_epoch_data = np.append(accepted_trough_epoch_data, block_accepted_trough_epoch_data, axis=1)
            all_trough_epoch_data = np.append(all_trough_epoch_data, block_all_trough_epoch_data, axis=1)

            block_all_random_epoch_data, _ = random_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, ms_per_sample, True) 
            all_random_epoch_data = np.append(all_random_epoch_data, block_all_random_epoch_data, axis=1)

    # save subject level data across all blocks
    save_dict = {"accepted_peak_epochs": accepted_peak_epoch_data, 
                "accepted_trough_epochs": accepted_trough_epoch_data,
                "accepted_dilation_epochs": accepted_dilation_epoch_data, 
                "accepted_constriction_epochs": accepted_constriction_epoch_data, 
                "all_peak_epochs": all_peak_epoch_data, 
                "all_trough_epochs": all_trough_epoch_data,
                "all_dilation_epochs": all_dilation_epoch_data, 
                "all_constriction_epochs": all_constriction_epoch_data, 
                "random_epochs": all_random_epoch_data}        
    outfile = os.path.join(results_dir, "event_statistics.pckl")

    with open(outfile, 'wb') as pickle_file:
        pickle.dump(save_dict, pickle_file)

    # plot mean time courses for each subject - accepted and all events 
    plot_mean_timecourses(half_epoch_duration=half_epoch_duration_ms, title = "Accepted events - subject "+subjID, 
                     peak_epoch = accepted_peak_epoch_data, trough_epoch = accepted_trough_epoch_data, 
                     constriction_epoch = accepted_constriction_epoch_data, dilation_epoch = accepted_dilation_epoch_data, 
                     random_epoch = all_random_epoch_data, save_dir= os.path.join(results_dir, "accepted_trace.png"))
    
    plot_mean_timecourses(half_epoch_duration=half_epoch_duration_ms, title = "All events - subject "+subjID, 
                     peak_epoch = all_peak_epoch_data, trough_epoch = all_trough_epoch_data, 
                     constriction_epoch = all_constriction_epoch_data, dilation_epoch = all_dilation_epoch_data, 
                     random_epoch = all_random_epoch_data, save_dir= os.path.join(results_dir, "all_trace.png"))
