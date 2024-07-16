import numpy as np
from scipy.signal import find_peaks
from psychopy import logging, core
import pylink
import config
from PsychoPy_funcs import quit_task

class StimulusDecider():
    """Object to store data and detect pupil phase events in real time and in simulations 
    
    Attributes
    --------
    task_name : string
        name of task to be used with decider 
    online : boolean 
        whether object is used in real time data collection or in simulations
    baseline_duration_ms : int 
        duration of baseline window in milliseconds 
    max_search_window_duration_ms : int
        maximum length of search window before resetting search window (in milliseconds)
    pupil_sample_duration_ms : int
        duration of a pupil sample from the eye tracker in real time, in milliseconds 
    num_random_events : int
        how many random events per block 
    random_event_time_sec : int 
        how long between random events 
    IEI_duration_sec : int 
        inter-event interval (i.e. how long should we wait between pupil events), in seconds
    search_window : list
        set of pupil samples to be used to find pupil events 
    search_window_sample_times : list
        list of timestamps (ms) associated with search windows
    baseline_window : list
        set of pupil samples to be used to determine thresholds for finding pupil samples 
    prior_search_window : list
        previous search window; for online algorithm, must have no blinks for current pupil 
        sample to be valid 
    search_window_model_fits : list
        final values from model fit on pupil samples 
    peak_threshold_var : float
        absolute threshold for identifying a peak event - event must be larger than this value 
        to be identified as a peak. Taken from config.py file. 
    trough_threshold_var : float
        absolute threshold for identifying a trough event - event must be smaller than this value
        to be identified as a trough. Taken from config.py file. 
    dilation_threshold : float
        absolute threshold for identifying a dilation event - first derivative of search window 
        must be larger than this value to be identified as a dilation. Taken from config.py file.  
    constriction_threshold : float
        absolute threshold for identifying a constriction event - first derivative of search window
        must be smaller than this value to be identified as a constriction. Taken from config.py file. 
    peak_pupil_quantile : float
        quantile threshold for identifying peak values - pupil size must be above this quantile 
        of values in baseline window to be accepted as a peak 
    trough_pupil_quantile : float
        quantile threshold for identifying troughs - pupil size must be below this quantile of 
        values in baseline window to be accepted as a trough 
    dilation_quantile : float
        quantile threshold for identifying dilations - used to define absolute value threshold
    constriction_quantile : float
        quantile threshold for identifying constrictions - used to define absoulte value threhold
    new_sample : obj
        most recent EyeLink sample (collected in real time)
    old_sample : obj
        previous Eyelink sample (collected in real time)
    pupil_sample_duration_time : int 
        duration of pupil sample collected in real time 
    peak_count : int
        number of peaks identified in real time
    trough_count : int
        number of troughs indentified in real time
    dilation_count : int
        number of dilations identified in real time 
    constriction_count : int 
        number of constrictions identified in real time 
    idx_event : int
        flag for event that was identified. 
        1 = peak, -1 = trough, 2 = dilation, -2 = constriction, 0 = no event
    accepted_pupil_event_bool : boolean 
        internal variable for marking whether an event was accepted 

    Methods
    --------
    set_pupil_phase_thresholds():
        Update thresholds for determining pupil phase events according to current baseline window
    update_windows(sample, pupil_sample_timer=None): 
        Update search and baseline windows with values from pupil sample 
    validate_search_window(max_search_window_duration_samples):
        Determine whether search window meets criteria for finding pupil phase events  
    fit_polynomial(demeaned_pupil_sample): 
        Do polynomial fit on pupil sample and save final value
    log_found_event_live(kind, demeaned_search_window, diff_fit):
        Interface with eye-tracker to log an event in real-time
    find_pupil_phase_event(pupil_sample_num=float('nan'), samples_in_pupil_sample=6, 
        current_time = float("nan"), peak_events = None, trough_events = None, 
        constriction_events = None, dilation_events = None):
        Find possible pupil phase events and log them 
    validate_event_offline(all_event_times, accepted_pupil_event_times, IEI_jitter_ms, 
        peak_events, trough_events, dilation_events, constriction_events):
        Determine whether enough time has passed to accept an identified pupil event in simulation
    validate_event_online(IEI_duration_ms):
        Determine whether enough time has passed to accept an identified pupil event in real-time
    reset_baseline_window(): 
        Clear baseline window 
    reset_search_window():
        Clear search window and associated variables
    accepted_pupil_event(): 
        Log an accepted pupil event. Function can be used for building closed-loop paradigms where
        detected pupil events trigger some sort of other task event.
    detect_events_online():  
        Detect pupil events in real-time. Validate search window, update pupil phase thresholds (if
        necessary) and identify pupil phase event in search window. 
    build_search_window(): 
        Build baseline and search windows in real-time. 
    get_pupil_sample_duration_time(): 
        Used to log pupil sample duration array. 
    get_search_window(): 
        accesses pupil samples in currently stored search window 
    get_search_window_times():
        accesses times associated with currently stored pupil sample, as the EyeLink eye-tracker 
        takes samples stochastically. 
    get_baseline_window(): 
        accesses current baseline window 
    get_search_window_fit_vals(): 
        returns the fitted search window values. 
    
    Notes
    --------
    This class can be used both in real-time data collection and in simulations. As such, some 
    attributes will be not be used in all cases. 

    Default values for quantiles and thresholds reflect those used in Kroenemer et al., 2024. 
    """
    
    def __init__(
        self,task_name, block_duration_sec=600, 
        baseline_duration_ms=5000, max_search_window_duration_ms=5000,
        pupil_sample_duration_ms=100, num_random_event=20,
        IEI_duration_sec=3, peak_pupil_quantile=0.75, trough_pupil_quantile=0.25, 
        dilation_quantile=0.99, constriction_quantile=0.01, online=False, win=None
    ):
        """
        Sets up a StimulusDecider object with defaults that reflect options from Kronemer et al., 2024. 

        Parameters
        ----------
        task_name: str
            name of task to be used with decider 
        block_duration_sec : int 
            number of seconds in a block 
        baseline_duration_ms : int
            duration of baseline window in milliseconds 
        max_search_window_duration_ms : int
            maximum length of search window in ms
        pupil_sample_duration_ms : int 
            duration of real-time pupil sample from eye-tracker in milliseconds              
        num_random_events : int 
            number of random events to include per block 
        IEI_duration_sec : int
            duration of inter-event interval in seconds 
        peak_pupil_quantile : float 
            quantile threshold for identifying peak values - pupil size must be above this quantile 
            of values in baseline window to be accepted as a peak 
        trough_pupil_quantile : float 
            quantile threshold for identifying troughs - pupil size must be below this quantile of 
            values in baseline window to be accepted as a trough 
        dilation_quantile : float
            quantile threshold for identifying dilations - used to define absolute value threshold
        constriction_quantile : float
            quantile threshold for identifying constrictions - used to define absoulte value threhold
        online : boolean 
            whether object is used in real time data collection or in simulations    
        win : PsychoPy Window
            screen to be updated in real-time 
        """
    
        self._baseline_duration_ms = baseline_duration_ms
        self._max_search_window_duration_ms = max_search_window_duration_ms
        self._pupil_sample_duration_ms = pupil_sample_duration_ms
        self._num_random_event = num_random_event
        self._random_event_time_sec = block_duration_sec/num_random_event
        self._IEI_duration_sec = IEI_duration_sec

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
        self._peak_threshold_var = config.peak_threshold
        self._trough_threshold_var = config.trough_threshold
        self._dilation_threshold = config.dilation_threshold
        self._constriction_threshold = config.constriction_threshold
        self._peak_pupil_quantile = peak_pupil_quantile
        self._trough_pupil_quantile = trough_pupil_quantile 
        self._dilation_quantile = dilation_quantile
        self._constriction_quantile = constriction_quantile
        self._idx_event = 0
        self._online = online
        self._task_name = task_name
        self._win = win

        if self._online: 
            self._pupil_phase_IEI_timer = core.Clock() # Inter-event interval timer
            self._random_IEI_timer = core.Clock() # Time for selecting random pupil times
            self._pupil_sample_timer = core.Clock() # Pupil sample timer

    def update_windows(self, sample, duration = None): 
        """ 
        Update search and baseline windows with values from pupil sample. 

        Replaces 0s with NaNs and then adds pupil sample to search window and baseline window, in 
        addition to adding sample times to sample time collector. If live data collection, also 
        records duration of pupil sample, as live eye-tracker data is stochastic and duration of 
        samples is not consistent. 

        Parameters
        ----------
        sample :  numpy.ndarray
            pupil sample to be to be added to baseline and search windows. Dimenions 2 x # of 
            time points - first row = time, second row = pupil size
        duration : float
            duration of pupil sample, if collected in real time
    
        """
        sample[1, sample[1,:] == 0] = float("nan")
        self._search_window.extend(sample[1,:])
        self._search_window_sample_times.extend(sample[0,:])
        self._baseline_window.extend(sample[1,:])

        if duration is not None:
            self._pupil_sample_duration_time.append(duration)

    def set_pupil_phase_thresholds(self):
        """Update thresholds for determining pupil phase events according to current baseline window
        
        Demeans baseline window (stored internally in self._baseline_window) and identifies peaks,
        troughs, dilations and constrictions in baseline window. Peaks and troughs are identified 
        through find_peaks, dilations and constrictions are identified using the gradient of 
        subsequent pupil sizes. Thresholds are then updated according to quantiles provided 
        at StimulusDecider initialization. Once thresholds are updated, baseline window is reset. 

        If data is collected in realtime, baseline window and updated thresholds are saved to logfile. 
            
        """
        
        # Demean baseline window
        demeaned_baseline_window = list(self._baseline_window - np.nanmean(self._baseline_window))
        
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

    def validate_search_window(self, max_search_window_duration_samples): 
        """ 
        Determine whether search window meets criteria for finding pupil phase events 

        Search window will be reset if it is too long or if there are blinks detected in the 
        current or previous search windows

        Parameters
        ----------
        max_search_window_duration_samples :  int
            maximum number of pupil samples to keep in search window

        Returns
        ----------
        boolean 
            whether or not search window was considered valid 
        """
        # If search window is empty
        if len(self._search_window) == 0:
            self._accepted_pupil_event_bool = False
            return False
        
        # If search window is too long - reset it
        if len(self._search_window) > max_search_window_duration_samples:
            self.reset_search_window()
            return False
        
        # Reset search window if a NaN sample is detected (e.g., blink event)
        if any(np.isnan(self._search_window)):
            self.reset_search_window()
            return False
        
        # Reset search window if a NaN sample was found in a prior search window
        if self._online:
            if any(np.isnan(self._prior_search_window)):
                self.reset_search_window()
                return False
        
        return True
    
    def fit_polynomial(self, demeaned_pupil_sample):
        """
        Do polynomial fit on pupil sample and save final value

        Pupil sample fit with a second order polynomial and last fitted value from sample is saved 
        to attribute self._search_window_model_fits (list)

        Parameters
        ----------
        demeaned_pupil_sample : list 
            sample to fit polynomial on. Note: should already be demeaned
        
        Raises
        ----------
        AssertionError
            mean of pupil sample is not 0, suggesting that the pupil sample may not have
            been demeaned prior to running method 
        """
        assert np.isclose(np.mean(demeaned_pupil_sample),0.0), "Pupil sample may not be demeaned"
        sample_window_fit_coef = np.polyfit(list(range(1,len(demeaned_pupil_sample)+1)),demeaned_pupil_sample, 2)

        # Find last sample fit value
        fit_value = np.polyval(sample_window_fit_coef,len(demeaned_pupil_sample))

        # Store the last value of fitting curve
        self._search_window_model_fits = np.append(self._search_window_model_fits,fit_value)
    
    def log_found_event_live(self, kind, demeaned_search_window, diff_fit): 
        """
        Interface with eye-tracker to log an event in real-time

        Internally (i.e. in StimulusDecider object), update count of event types. 

        Parameters
        ----------
        kind : str
            kind of event to log 
        demeaned_search_window : np.ndarray 
            search window that is being fit 
        diff_fit : np.ndarray 
            gradient of search window being fit 
        """
        el_tracker = pylink.getEYELINK()
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
        """Find possible pupil phase events and log them 

        Function can be used in both simulations (in which case all event items should not be None)
        or in real-time data collection (in which case self._online = True and event objects should
        be None, as is the default). 

        If event is identified, update self._idx_event to the id of that event.
        
        Parameters
        ----------
        pupil_sample_num : float 
            index of pupil sample to identify events in 
        samples_in_pupil_sample : int 
            how many size samples in pupil sample 
        current_time : float
            time (ms) of pupil sample 
        peak_events : EventCollector object 
            object to log peak events and associated information 
        trough_events : EventCollector object 
            object to log trough events and associated information 
        constriction_events : EventCollector object 
            object to log constriction events and associated information 
        dilation_events : EventCollector object 
            object to log dilation events and associated information

        Returns 
        ----------
        int 
            integer reflecting kind of event detected. 1 = peak, -1 = trough, 
            2 = dilation, -2 = constriction, 0 = no event 
        
        """

        demeaned_search_window =  list(self._search_window - np.nanmean(self._search_window))

        if self._online: 
            self.fit_polynomial(demeaned_search_window)

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
    
    def validate_event_offline(self, all_event_times, accepted_pupil_event_times, IEI_jitter_ms, 
                   peak_events, trough_events, dilation_events, constriction_events):
        """
        Determine whether enough time has passed to accept an identified pupil event in simulation

        If event is accepted, update appropriate EventCollector object, record time and reset
        search window. 

        Parameters
        ----------
        all_event_times : list 
            list of all event times (regardless of type)
        accepted_pupil_event_times : list 
            list of accepted pupil event times
        IEI_jitter_ms : int
            amount of time (ms) that must have passed to be considered an accepted event 
        peak_events : EventCollector object 
            object to log peak events and associated information 
        trough_events : EventCollector object 
            object to log trough events and associated information 
        constriction_events : EventCollector object 
            object to log constriction events and associated information 
        dilation_events : EventCollector object 
            object to log dilation events and associated information

        Returns
        ----------
        peak_events : EventCollector object 
            object to log peak events with updated accepted event info
        trough_events : EventCollector object 
            object to log trough events with updated accepted event info
        constriction_events : EventCollector object 
            object to log constriction events with updated accepted event info
        dilation_events : EventCollector object 
            object to log dilation events with updated accepted event info
        """
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

    def validate_event_online(self, IEI_jitter_ms):
        """ 
        Determine whether enough time has passed to accept an identified pupil event in real-time

        Pulls time from global pupil_phase_IEI_event_timer. If enough time has passed, log event
        and reset search window. 

        Parameters
        ----------
        IEI_jitter_ms : int
            amount of time (ms) that must have passed to be considered an accepted event 
        
        """
        el_tracker = pylink.getEYELINK()
        if self._pupil_phase_IEI_event_timer.getTime() > IEI_jitter_ms/1000: 
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

    def reset_baseline_window(self): 
        """ Clear baseline window """
        self._baseline_window = []
        
    def reset_search_window(self):
        """
        Clear search window and associated variables 
        """
        self._search_window = []
        self._search_window_sample_times = []
        self._search_window_model_fits = []
        if self._online: 
            self._prior_search_window = self._search_window
            self._accepted_pupil_event_bool = False
    
    def set_current_event_idx(self, found_event): 
        self._idx_event = found_event

    def accepted_pupil_event(self):
        """
        Log an accepted pupil event (i.e., event detected beyond the inter-event interval
        Note: This function can be used for building closed-loop paradigms where
        detected pupil phase events trigger a task event. In this case, you might have another 
        PsychoPy function that does some sort of event that is called here. 

        Once an event is accepted, resets global pupil_phase_IEI_time, changes internal 
        accepted_pupil_event_bool to True and sends a logging event to PsychoPy logs. 

        """

        el_tracker = pylink.getEYELINK()

        # Log
        logging.log(level=logging.EXP,msg='Accepted Pupil Event')
        el_tracker.sendMessage('Accepted Pupil Event')
        
        # Reset pupil phase IEI timer
        self._pupil_phase_IEI_timer.reset()
        
        # Set boolean to True
        self._accepted_pupil_event_bool = True

    def detect_events_online(self): 
        """
        Detect pupil events in real-time. Search window is validated, then pupil phase threshold values are updated. 
        Once threshold values are updated, determine whether enough time has passed for a random event. If not, 
        pupil phase events are looked for. If event is identified, reset timers. 

        Returns 
        ----------
        int 
            integer reflecting kind of event detected. 1 = peak, -1 = trough, 
            2 = dilation, -2 = constriction, 3 = random event, 0 = no event 

        """
        el_tracker = pylink.getEYELINK()

        valid_window = self.validate_search_window(self._max_search_window_duration_ms//config.ms_per_sample)
        if not valid_window:
            return 0
        demeaned_search_window = list(self._search_window - np.mean(self._search_window))

        # Reset pupil phase pupil size and pupil size derivative thresholds - once the minimum baseline duration has been acquired
        if len(self._baseline_window) > round(self._baseline_duration_ms//config.ms_per_sample):
            self.set_pupil_phase_thresholds(self._baseline_window)
        
        # Log random event if the minimum random IEI is exceed      
        if self._random_IEI_timer.getTime() > self._random_event_time_sec:
        
            # If IEI has been exceed (Note: The IEI timer gets reset in the accepted_pupil_event; 
            # no reset will happen if no stimulus is shown)
            if self._pupil_phase_IEI_timer.getTime() > self._IEI_duration_sec:
                
                # Define idx extrema
                self._idx_event = 3
                
                # Log
                logging.log(level=logging.EXP,msg='Random Event')
                el_tracker.sendMessage('Random Event') 
                
                # Accepted event
                self.accepted_pupil_event()

                # Reset timers
                self._random_IEI_timer.reset()
                #self._general_timer.reset() # CW: why general timer here? 
                
                return self._idx_event
            
            return 0
        
        # Find peaks (local maxima), troughs (local minima), dilation (dilation), and constriction (constriction) events
        # Index extrema dictionary: peak = 1; dilation = 2; trough = -1; constriction = -2; random = 3
        self._idx_event = self.find_pupil_phase_event(demeaned_search_window)
        
        # If a pupil phase event was found
        if self._idx_event!=0:
            
            # Confirm pupil phase IEI exceeded
            if self._pupil_phase_IEI_timer.getTime() > self._IEI_duration_sec:
                
                # Accepted event
                self.accepted_pupil_event()
                idx_to_return = self._idx_event
                # Reset search window 
                self.reset_search_window()
                self._idx_event = 0
                return idx_to_return
            
            # If IEI is not exceeded 
            elif self._pupil_phase_IEI_timer.getTime() < self._IEI_duration_sec:
                
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

    def build_search_window(self) -> float:
        """Add the pupil sample to baseline and search window. 
        
        """
        el_tracker = pylink.getEYELINK()

        # Initialize/reset pupil_sample variable
        pupil_sample = []
        pupil_sample_time = []
        
        # Setup pupil size/time variables
        p_size = float("nan")
        p_time = float("nan")
        
        # Reset pupil sample timer
        
        self._pupil_sample_timer.reset()
        
        # Keep adding pupil values until pupil_sample is long enough      
        while len(pupil_sample) < self._pupil_sample_duration_ms//config.ms_per_sample:
        
            # Update window
            if self._win is not None:
                self._win.update()
        
            # Get the latest EyeLink sample        
            self._new_sample = el_tracker.getNewestSample()
            if self._new_sample is not None:
                  if self._old_sample is not None:

            # Check there is a new pupil value and that it is not the first pupil sample 
           # if self._new_sample is not None and self._old_sample is not None:
                    
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
            quit_task(self._win)
        
            # Replace old pupil value with new
            self._old_sample = self._new_sample
            
        # Add pupil sample and pupil sample times to search window and search window sample times 
        self._search_window.extend(pupil_sample)
        self._search_window_sample_times.extend(pupil_sample_time)
        
        # Store duration of pupil sample in time
        self._pupil_sample_duration_time.append(self._pupil_sample_timer.getTime())
        
        # Add pupil sample to baseline window
        self._baseline_window.extend(pupil_sample)

    def get_pupil_sample_duration_time(self):
        """Used to log pupil sample duration array."""
        return self._pupil_sample_duration_time

    def get_search_window(self): 
        return self._search_window

    def get_search_window_times(self):
        return self._search_window_sample_times

    def get_baseline_window(self):
        return self._baseline_window

    def get_search_window_fit_vals(self): 
        return self._search_window_model_fits
