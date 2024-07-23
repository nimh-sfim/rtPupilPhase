import numpy as np 
import os 
import pickle
import argparse

from utils import get_data, find_string_time, pull_pupil_sample, plot_mean_timecourses
from StimulusDecider import StimulusDecider
from EventCollector import EventCollector

import config

def main(subject_list, pupil_sample_duration_ms, num_random_events, IEI_duration_ms, baseline_duration_ms,
         max_search_window_duration_ms, half_epoch_duration_ms, plot_timecourses):
    
    assert len(subject_list) > 0, "Please provide at least one subject to simulate"
    
    # number of eyelink samples in pupil sample
    samples_in_pupil_sample = np.round(pupil_sample_duration_ms/config.ms_per_sample) 
    random_IEI = config.block_duration_ms/num_random_events # duration of IEI in milliseconds
    # Inter-event interval for pupil phase events
    samples_in_IEI = np.round(IEI_duration_ms/config.ms_per_sample) # in samples

    # calculate number of samples in baseline window  
    samples_in_baseline_window = np.round(baseline_duration_ms/config.ms_per_sample) 

    # define paths 
    root_path = os.getcwd()
    data_base = os.path.join(root_path,config.data_fname)
    results_base = os.path.join(root_path,config.results_fname)

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
        for block in range(1,config.num_blocks+1): 
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
            downsampled_block_data = block_data[:, 0::config.downsample_value]

            # reset StimulusDecider for new block 
            sd.reset_baseline_window()
            sd.reset_search_window()

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
                if not sd.validate_search_window(max_search_window_duration_ms//config.ms_per_sample):
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
                        peak_events, trough_events, dilation_events, constriction_events = sd.validate_event_offline(event_times, accepted_pupil_event_times, IEI_duration_ms, 
                                                peak_events, trough_events, dilation_events, constriction_events)
                        
            
            # Note that it may be recommended to perform additional blink and microsaccade detection using your method of choice here. 
            # In the MATLAB simulations, we interpolate over blinks and use a method to identify microsaccades based off of Engbert & Kliegl (2003)     
            # Here, we do not perform this data processing and simply plot the data. 

            # pull epoch data for each event and concatentate across block 
            if block == 1:
                accepted_constriction_epoch_data, all_constriction_epoch_data = constriction_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, config.ms_per_sample)
                accepted_dilation_epoch_data, all_dilation_epoch_data = dilation_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, config.ms_per_sample)
                accepted_peak_epoch_data, all_peak_epoch_data = peak_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, config.ms_per_sample)
                accepted_trough_epoch_data, all_trough_epoch_data = trough_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, config.ms_per_sample)
                all_random_epoch_data, _ = random_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, config.ms_per_sample) 

            else: 
                block_accepted_constriction_epoch_data, block_all_constriction_epoch_data = constriction_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, config.ms_per_sample)
                accepted_constriction_epoch_data = np.append(accepted_constriction_epoch_data, block_accepted_constriction_epoch_data, axis=1)
                all_constriction_epoch_data = np.append(all_constriction_epoch_data, block_all_constriction_epoch_data, axis=1)

                block_accepted_dilation_epoch_data, block_all_dilation_epoch_data = dilation_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, config.ms_per_sample)
                accepted_dilation_epoch_data = np.append(accepted_dilation_epoch_data, block_accepted_dilation_epoch_data, axis=1)
                all_dilation_epoch_data = np.append(all_dilation_epoch_data, block_all_dilation_epoch_data, axis=1)

                block_accepted_peak_epoch_data, block_all_peak_epoch_data = peak_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, config.ms_per_sample)
                accepted_peak_epoch_data = np.append(accepted_peak_epoch_data, block_accepted_peak_epoch_data, axis=1)
                all_peak_epoch_data = np.append(all_peak_epoch_data, block_all_peak_epoch_data, axis=1)

                block_accepted_trough_epoch_data, block_all_trough_epoch_data = trough_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, config.ms_per_sample)
                accepted_trough_epoch_data = np.append(accepted_trough_epoch_data, block_accepted_trough_epoch_data, axis=1)
                all_trough_epoch_data = np.append(all_trough_epoch_data, block_all_trough_epoch_data, axis=1)

                block_all_random_epoch_data, _ = random_events.pull_valid_epochs(block_data[1,:], half_epoch_duration_ms, config.ms_per_sample) 
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

        if plot_timecourses:
            # plot mean time courses for each subject - accepted and all events 
            plot_mean_timecourses(half_epoch_duration=half_epoch_duration_ms, title = "Accepted events - subject "+subjID, 
                            peak_epoch = accepted_peak_epoch_data, trough_epoch = accepted_trough_epoch_data, 
                            constriction_epoch = accepted_constriction_epoch_data, dilation_epoch = accepted_dilation_epoch_data, 
                            random_epoch = all_random_epoch_data, save_dir= os.path.join(results_dir, "accepted_trace.png"))
            
            plot_mean_timecourses(half_epoch_duration=half_epoch_duration_ms, title = "All events - subject "+subjID, 
                            peak_epoch = all_peak_epoch_data, trough_epoch = all_trough_epoch_data, 
                            constriction_epoch = all_constriction_epoch_data, dilation_epoch = all_dilation_epoch_data, 
                            random_epoch = all_random_epoch_data, save_dir= os.path.join(results_dir, "all_trace.png"))


#subject_list = ['046','048','073','074','078','079','080','081'] # command line 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description= 'Simulated rtPupilPhase: simulating real-time pupillometry')
    parser.add_argument("subs", help="Subjects to simulate", 
                        nargs='+', default=[]) 
    parser.add_argument("--plot_timecourses", help="Whether or not to display mean timecourses.", 
                        action="store_true")
    parser.add_argument("--pupil_sample_duration_ms", help="Duration of pupil sample. Default: 100ms",
                        type = int, default=100)
    parser.add_argument("--num_random_events", help="Number of random events to occur per block. Note that this value will have implications on the number of pupil phase events. Default: 15",
                        type=int, default=15)
    parser.add_argument("--baseline_duration_ms", help="Duration of baseline window in milliseconds. Default: 5000ms",
                        type=int, default=5000)
    parser.add_argument("--max_search_window_duration_ms", help="Maximum duration of search window before resetting, in milliseconds. Default: 5000ms",
                         type=int, default=5000)
    parser.add_argument("--IEI_duration_ms", help="Inter-event interval - how long to wait between valid events in milliseconds. Default: 3000ms", 
                        type=int, default=3000)
    parser.add_argument("--half_epoch_duration_ms", help="Half epoch duration - how far before and after event to plot. Default: 2500ms",
                        type=int, default=2500)

    args = parser.parse_args()
    main(args.subs, args.pupil_sample_duration_ms, args.num_random_events, args.IEI_duration_ms, 
         args.baseline_duration_ms, args.max_search_window_duration_ms, args.half_epoch_duration_ms, 
         args.plot_timecourses)