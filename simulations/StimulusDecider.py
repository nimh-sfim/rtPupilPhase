import numpy as np
from scipy.signal import find_peaks
from psychopy import logging


class StimulusDecider():
    """Object to store data and detect pupil phase events in real time."""
    
    def __init__(
        self,task_name,peak_pupil_quantile=0.75, trough_pupil_quantile=0.25, 
        dilation_quantile=0.99, constriction_quantile=0.01, online=False
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
        self._peak_pupil_quantile = peak_pupil_quantile
        self._trough_pupil_quantile = trough_pupil_quantile 
        self._dilation_quantile = dilation_quantile
        self._constriction_quantile = constriction_quantile
        self._idx_event = 0
        self._online = online
        self._task_name = task_name

    def update_windows(self, sample, pupil_sample_timer=None): 
        sample[1, sample[1,:] == 0] = float("nan")
        self._search_window.extend(sample[1,:])
        self._search_window_sample_times.extend(sample[0,:])
        self._baseline_window.extend(sample[1,:])

        if pupil_sample_timer is not None:
            self._pupil_sample_duration_time.append(pupil_sample_timer.getTime())

    def set_pupil_phase_thresholds(self, baseline_window):
        """Reset pupil phase event thresholds according to current baseline window"""
        
        # Demean baseline window
        demeaned_baseline_window = list(baseline_window - np.nanmean(baseline_window))
        
        # Log
        if self._online:
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
        retro_peaks_quantile = np.quantile(retro_peaks, self._peak_pupil_quantile)
        retro_troughs_quantile = np.quantile(retro_troughs, self._trough_pupil_quantile)
        
        # Find dilation and constriction threshold quantiles
        retro_dilation_quantile = np.quantile(diff_demeaned_baseline_window, self._dilation_quantile)
        retro_constriction_quantile = np.quantile(diff_demeaned_baseline_window, self._constriction_quantile)
        
        
        # Set new thresholds
        self._peak_threshold_var = retro_peaks_quantile
        self._trough_threshold_var = retro_troughs_quantile
        self._dilation_threshold = retro_dilation_quantile
        self._constriction_threshold = retro_constriction_quantile
        
        # Reset baseline window
        self._baseline_window = []
        
        # Log
        if self._online:
            logging.log(level=logging.EXP,msg='Updated peak threshold: ' + str(self._peak_threshold_var))
            logging.log(level=logging.EXP,msg='Updated trough threshold: ' + str(self._trough_threshold_var))
            logging.log(level=logging.EXP,msg='Updated dilation threshold: ' + str(self._dilation_threshold))
            logging.log(level=logging.EXP,msg='Updated constriction threshold: ' + str(self._constriction_threshold))
        
    def get_baseline_window(self):
        return self._baseline_window
    
    def reset_baseline_window(self): 
        self._baseline_window = []

    def validate_search_window(self, max_search_window_duration_samples): 
        if len(self._search_window) == 0:
            self._accepted_pupil_event_bool = False
            return 0
        
        # If search window is too long - reset it
        if len(self._search_window) > max_search_window_duration_samples:
            self.reset_search_window()
            return 0
        
        # Reset search window if a NaN sample is detected (e.g., blink event)
        if any(np.isnan(self._search_window)):
            self.reset_search_window()
            return 0
        
        # Reset search window if a NaN sample was found in a prior search window
        if self._online:
            if any(np.isnan(self._prior_search_window)):
                self.reset_search_window()
                return 0
        
        return 1
    
    def get_search_window(self): 
        return self._search_window
    
    def get_search_window_times(self):
        return self._search_window_sample_times
        
    def reset_search_window(self):
        """
        
        """
        self._search_window = []
        self._search_window_sample_times = []
        self._search_window_model_fits = []
        if self._online: 
            self._prior_search_window = self._search_window
            self._accepted_pupil_event_bool = False
    
    def fit_polynomial(self, pupil_sample):
        sample_window_fit_coef = np.polyfit(list(range(1,len(pupil_sample)+1)),pupil_sample, 2)

        # Find last sample fit value
        fit_value = np.polyval(sample_window_fit_coef,len(pupil_sample))

        # Store the last value of fitting curve
        self._search_window_model_fits = np.append(self._search_window_model_fits,fit_value)
    
    def get_search_window_fit_vals(self): 
        return self._search_window_model_fits
    
    def log_found_event_live(self, el_tracker, kind, demeaned_search_window, diff_fit): 
        if kind == "trough": 
                logging.log(level=logging.EXP,msg='Search Window Pupil: ' + str(demeaned_search_window))
                logging.log(level=logging.EXP,msg='Search Window Sample Time: ' + str(self._search_window_sample_times))
                logging.log(level=logging.EXP,msg='Search Window Fit Diff: ' + str(diff_fit))
                logging.log(level=logging.EXP,msg='Found Trough')
                el_tracker.sendMessage('Found Trough')

                self._trough_count += 1
        
        elif kind == "peak":
            logging.log(level=logging.EXP,msg='Search Window Pupil: ' + str(demeaned_search_window))
            logging.log(level=logging.EXP,msg='Search Window Sample Time: ' + str(self._search_window_sample_times))
            logging.log(level=logging.EXP,msg='Search Window Fit Diff: ' + str(diff_fit))
            logging.log(level=logging.EXP,msg='Found Peak')
            el_tracker.sendMessage('Found Peak')
            
            self._peak_count += 1

        elif kind == "constriction":
            logging.log(level=logging.EXP,msg='Search Window Pupil: ' + str(demeaned_search_window))
            logging.log(level=logging.EXP,msg='Search Window Sample Time: ' + str(self._search_window_sample_times))
            logging.log(level=logging.EXP,msg='Search Window Fit Diff: ' + str(diff_fit))
            logging.log(level=logging.EXP,msg='Found Constriction Event')
            el_tracker.sendMessage('Found Constriction Event')

            self._constriction_count += 1
            
        elif kind == "dilation":             
            logging.log(level=logging.EXP,msg='Search Window Pupil: ' + str(demeaned_search_window))
            logging.log(level=logging.EXP,msg='Search Window Sample Time: ' + str(self._search_window_sample_times))
            logging.log(level=logging.EXP,msg='Search Window Fit Diff: ' + str(diff_fit))
            logging.log(level=logging.EXP,msg='Found Dilation Event')
            el_tracker.sendMessage('Found Dilation Event')
            
            self._dilation_count += 1

    def find_pupil_phase_event(self, pupil_sample_num=float('nan'), samples_in_pupil_sample=6, current_time = float("nan"),
                            peak_events = None, trough_events = None, 
                            constriction_events = None, dilation_events = None) -> int:
        """Find possible pupil phases"""

        demeaned_search_window =  list(self._search_window - np.nanmean(self._search_window))

        # Diff across fit values if there are more than one fit values
        if len(self._search_window_model_fits) > 1:
            
            # Find diff of fit values
            diff_fit = list(np.diff(self._search_window_model_fits))

            # *** Find trough event ***
            if diff_fit[-1] > 0 and demeaned_search_window[-1] < self._trough_threshold_var:
            
                if self._online:
                    self.log_found_event_live(kind = "trough")
                else: 
                    trough_events.update_data(((pupil_sample_num+1)*samples_in_pupil_sample), current_time, 
                                    self.get_search_window()[-1], diff_fit[-1])

                return -1

            # *** Find peak event ***
            elif diff_fit[-1] < 0 and demeaned_search_window[-1] > self._peak_threshold_var:
                
                if self._online:
                    self.log_found_event_live(kind = "peak")
                else: 
                    peak_events.update_data(((pupil_sample_num+1)*samples_in_pupil_sample), current_time, 
                                    self.get_search_window()[-1], diff_fit[-1])
                
                return 1
            
            # *** Find dilation event ***
            elif diff_fit[-1] > self._dilation_threshold:
                
                if self._online:
                    self.log_found_event_live(kind = "dilation")
                else: 
                    dilation_events.update_data(((pupil_sample_num+1)*samples_in_pupil_sample), current_time, 
                                    self.get_search_window()[-1], diff_fit[-1])
                    
                return 2
            
            # *** Find constriction event ***
            elif diff_fit[-1] < self._constriction_threshold:
                
                if self._online:
                    self.log_found_event_live(kind = "constriction")
                else: 
                    constriction_events.update_data(((pupil_sample_num+1)*samples_in_pupil_sample), current_time, 
                                    self.get_search_window()[-1], diff_fit[-1])

                return -2
            
            # No extrema found
            else: 
                if self._online: 
                    self._accepted_pupil_event_bool = False
                
                return 0
                
        # Not enough fit values to complete diff analysis
        else:
            if self._online:
                self._accepted_pupil_event_bool = False
            
            return 0
        
    def set_current_event_idx(self, found_event): 
        self._idx_event = found_event
    
    def validate_event_offline(self, all_event_times, accepted_pupil_event_times, IEI_jitter_ms, 
                   peak_events, trough_events, dilation_events, constriction_events):
        # get time from last events
        if len(all_event_times) > 1 and len(accepted_pupil_event_times) > 0:
            # calculate time from last event
            time_from_last_accepted_event = all_event_times[-1] - accepted_pupil_event_times[-1]
        else: 
            # guarantees first stimulus event triggers accepted event 
            time_from_last_accepted_event = IEI_jitter_ms
            
        # check if IEI is exceeded     
        if time_from_last_accepted_event >= IEI_jitter_ms:
            accepted_pupil_event_times.append(all_event_times[-1])
            if self._idx_event == 1: # peak
                peak_events.store_accepted_event()
            elif self._idx_event == -1: # trough
                trough_events.store_accepted_event() 
            elif self._idx_event == 2: # dilation
                dilation_events.store_accepted_event()
            elif self._idx_event == -2: # constriction
                constriction_events.store_accepted_event()
            
            self.reset_search_window()
            self._idx_event = 0
        
        return peak_events, trough_events, dilation_events, constriction_events

    def validate_event_online(self, IEI_duration_ms):

        if pupil_phase_IEI_event_timer.getTime() > IEI_duration_ms/1000: 
            self.accepted_pupil_event()
            self.reset_search_window()
            self.set_current_event_idx(0)
        else: 
            # Not an accepted event
            self._accepted_pupil_event_bool = False
                    
            # Log
            logging.log(level=logging.EXP,msg='Within IEI - Skipping this Pupil Event')
            el_tracker.sendMessage('Within IEI - Skipping this Pupil Event')
        
        # NOTE: search_window is not reset; will continue looking for a pupil phase events after adding another pupil_sample;
        # unless the search window gets too long
        
        # Not an accepted event
        self._accepted_pupil_event_bool = False              
