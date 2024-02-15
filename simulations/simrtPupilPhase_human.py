import numpy as np 
#from scipy.signal import find_peaks 
import os 
#import matplotlib.pyplot as plt
#from scipy.io import loadmat
import pickle

from utils import get_data, pull_pupil_sample, plot_mean_timecourses
from StimulusDecider import StimulusDecider
from EventCollector import EventCollector

# Dictionary:
# ms = milliseconds
# IEI = inter-event interval
# num = number
# EyeLink = the pupillometry system (SR Research, Inc.)
# idx = index

# *** Subject and Recording Parameters ***

# Subject list
#subject_list = ['046','048','073','074','078','079','080','081']
subject_list = ['048']

# The recorded eye (0 = left; 1 = right)
# Note: EyeLink stores the pupil size data in a 2 x time/sample matrix. The
# first row = the left eye and the second row = the right eye.
recorded_eye = 1

# Number of task blocks completed per participant
num_blocks = 5
block_duration_ms = 600000 # in milliseconds

# Define pupillometry sampling rate
sampling_rate = 1000 # in Hz
ms_per_sample = 17 # EyeLink live recording rate is (60Hz or ~17 samples/s)
downsample_value = 17 # Downsample value

# *** Display Parameters ***

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

# Warning flag (1 = yes; 0 = no)
# Note: There are multiple warnings scripted throughout the simulation that
# can help keep track of the simulation progress. However, for a cleaner
# terminal, you can supress these warnings. Scripted errors are not
# impacted by the warning flag.
warning_flag = 0

# 50# of the epoch length to be extracted
half_epoch_duration_ms = 2500 # in milliseconds

# No blinks or saccades interval
no_blinks_saccades = 500 # in milliseconds

# definepaths 
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

    if not os.path.isdir(results_dir):
        os.makedirs(results_dir) 

    sd = StimulusDecider("fixation")

    for block in range(1,num_blocks+1): 
        print("Starting block "+str(block))

        random_events = EventCollector("random")
        trough_events = EventCollector("trough")
        peak_events = EventCollector("peak")
        constriction_events = EventCollector("constriction")
        dilation_events = EventCollector("dilation")
        accepted_pupil_event_times = []

        start_time = event_data['sttime'][0][np.where(event_data['message'] == ["Starting Perception Rate Block "+str(block)])[1][0]][0][0]
        end_time = event_data['sttime'][0][np.where(event_data['message'] == ["Finished Perception Rate Block "+str(block)])[1][0]][0][0]

        block_data = gaze_data[:,np.where(gaze_data[0]==start_time)[0][0]:np.where(gaze_data[0]==end_time)[0][0]]
        
        # Downsampling matches th Eyelink real-time sampling in a live testing session
        # goes from 1000Hz to 60Hz 
        downsampled_block_data = block_data[:, 0::downsample_value]

        for pupil_sample_num in range(int(np.shape(downsampled_block_data)[1]/6)-1):
            # stage 1: get pupil sample 
            current_pupil_sample = pull_pupil_sample(downsampled_block_data, pupil_sample_num, samples_in_pupil_sample)

            # stage 2: create search window 
            sd.update_windows(current_pupil_sample)

            # if the baseline window is long enough (and not all nans), update thresholds
            if len(sd.get_baseline_window()) > samples_in_baseline_window:
                if not np.isnan(sd.get_baseline_window()).all(): 
                    sd.set_pupil_phase_thresholds(sd.get_baseline_window())
                else: 
                    sd.reset_baseline_window()

            # validate search window 
            if sd.validate_search_window(max_search_window_length_ms/ms_per_sample)==0:
                sd.reset_search_window()

                continue
            
            # demean search window 
            demeaned_search_window = sd.get_search_window() - np.mean(sd.get_search_window())

            sd.fit_polynomial(demeaned_search_window)

            # stage 4: find pupil events 

            # check if random event happened 
            current_time = sd.get_search_window_times()[-1]

            if len(random_events.get_times()) >= 1: 
                time_from_last_event = current_time - random_events.get_times()[-1]
            else: 
                # ensures that the first stimulus always triggers accepted random event
                time_from_last_event = random_IEI

            if time_from_last_event >= random_IEI:
                random_events.update_data(((pupil_sample_num+1)*samples_in_pupil_sample), current_time, 
                                        sd.get_search_window()[-1])
            
            # find pupil phase events 
            if len(sd.get_search_window()) > 1:
                found_event = sd.find_pupil_phase_event(pupil_sample_num=pupil_sample_num, current_time=current_time, peak_events=peak_events, 
                trough_events=trough_events, constriction_events=constriction_events, 
                dilation_events=dilation_events)
                sd.set_current_event_idx(found_event)

                event_times = peak_events.get_times() + trough_events.get_times() +  constriction_events.get_times() + dilation_events.get_times()
                event_times.sort()

                if found_event != 0:
                    peak_events, trough_events, dilation_events, constriction_events = sd.validate_event_offline(event_times, accepted_pupil_event_times, IEI_jitter_ms, 
                                            peak_events, trough_events, dilation_events, constriction_events)
                    
            
        # pull epoch data         
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

    # # save subject level data
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

    plot_mean_timecourses(half_epoch_duration=half_epoch_duration_ms, title = "Accepted events - subject "+subjID, 
                     peak_epoch = accepted_peak_epoch_data, trough_epoch = accepted_trough_epoch_data, 
                     constriction_epoch = accepted_constriction_epoch_data, dilation_epoch = accepted_dilation_epoch_data, 
                     random_epoch = all_random_epoch_data)
    
    plot_mean_timecourses(half_epoch_duration=half_epoch_duration_ms, title = "All events - subject "+subjID, 
                     peak_epoch = all_peak_epoch_data, trough_epoch = all_trough_epoch_data, 
                     constriction_epoch = all_constriction_epoch_data, dilation_epoch = all_dilation_epoch_data, 
                     random_epoch = all_random_epoch_data)


# TODO: refactor real-time code 
# TODO: documentation 


                            







             

        

        

    



