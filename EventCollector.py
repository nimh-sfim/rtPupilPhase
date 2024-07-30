import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

def get_data(sub_dir, recorded_eye=1): 
    """
    Find and load eye-tracking data

    Parameters
    ----------
    sub_dir : str
         directory where .mat file is located  
    recorded_eye : int
         column of data to pull data from. For Eyelink data, 0 = left eye, 1 = right eye 

    Returns
    -------
    all_data : numpy.ndarray 
        array of eye-tracking data of dimensions (4, number of timepoints). 
        Rows reflect time (ms), pupil size (pixels), gaze X location, gaze Y location
    event_data : numpy.ndarray)
        array of messages from EyeLink events, with dimension (2, number of events)
        Rows reflect the message sent to the EyeLink and the time (ms) associated with message. 

    Raises
    ------
    AssertionError 
        More than one pupil data file found in subject directory 
    Assertion Error
        Fewer than 10 unique values in pupillometry data, suggesting incorrect eye used 
    """
    edf_files =  [x for x in os.listdir(sub_dir) if x.endswith(".mat")]
    
    assert len(edf_files) ==1, "More than one pupil data file found!"

    filename = os.path.join(sub_dir, edf_files[0])
    data = loadmat(filename)
    edf_data = data['edf_data']

    time_data = np.squeeze(np.transpose(edf_data[0,0]['FSAMPLE'][0,0]['time']))
    pupil_data = edf_data[0,0]['FSAMPLE'][0,0]['pa'][recorded_eye,:]
    gaze_Y_data = edf_data[0,0]['FSAMPLE'][0,0]['gy'][recorded_eye,:]
    gaze_X_data = edf_data[0,0]['FSAMPLE'][0,0]['gx'][recorded_eye,:]
    all_data = np.array([time_data,pupil_data, gaze_X_data, gaze_Y_data])

    # Check the data is actually coming from the correct eye
    assert len(np.unique(pupil_data)) > 10, "Low pupil data value diversity - might be using the wrong eye!"

    all_data = np.array([time_data, pupil_data, gaze_X_data, gaze_Y_data])

    event_data = edf_data[0,0]["FEVENT"][['message', 'sttime']]

    return all_data, event_data

def find_string_time(time_array, message_array, match_string):
    """
    Find the time in ms of a given event, when there are empty arrays 

    Parameters
    ----------
    time_array : numpy.ndarray
        array of timestamps associated with messages
    message_array : numpy.ndarray
        array of messages to search against
    match_string : str
        the string to find

    Returns
    -------
    int 
        Time (ms) of event matching inputted string

    Raises
    ------
    ValueError
        Message does not exist in array 
    
    """
    item_idx = None
    for idx, item in enumerate(message_array):
        if item.size > 0:
            if item == match_string:
                item_idx = idx
    
    if item_idx is None:
        raise ValueError("String does not exist in array")
    else:   
        time = time_array[item_idx][0][0]
        return time

def pull_pupil_sample(data, pupil_sample_num, samples_in_pupil_sample): 
    """
    Pull a specific pupil sample from gaze data 
    
    Parameters
    ----------
    data : numpy.ndarray
        downsampled eyetracking data from block. First 2 rows should be time, pupil size data
    pupil_sample_num : int
        which pupil sample to pull
    samples_in_pupil_sample : int
        how many pupil size readings are in each pupil sample

    Returns
    -------
    numpy.ndarray 
        extracted pupil sample 
    """
    if pupil_sample_num ==0: 
        current_pupil_sample = data[0:2, 0:int(samples_in_pupil_sample)]
    else: 
        current_pupil_sample = data[0:2, int(pupil_sample_num*samples_in_pupil_sample):int((pupil_sample_num+1)*samples_in_pupil_sample)]
    return current_pupil_sample

def plot_mean_timecourses(half_epoch_duration, title = "", peak_epoch=None, 
                          trough_epoch=None, constriction_epoch=None, dilation_epoch=None, 
                          random_epoch=None, save_dir = None):
    """
    Plot the mean epoch timecourses. Will plot as many event types as is provided. 
    Must provide labeled arguments to ensure that the color legend is correct. 

    Parameters
    ----------
    half_epoch_duration : int
        duration (ms) of half an epoch to plot (i.e. how much time before and after event)
    title : str
        title of the plot
    peak_epoch : numpy.ndarray
        array of peak event epochs to average and plot. 
        Should be dimensions (# of events x 2*half_epoch_duration + 1)
    trough_epoch : numpy.ndarray
        array of trough event epochs to average and plot. 
        Should be dimensions (# of events x 2*half_epoch_duration + 1)
    constriction_epoch : numpy.ndarray
        array of constriction event epochs to average and plot. 
        Should be dimensions (# of events x 2*half_epoch_duration + 1)
    dilation_epoch : numpy.ndarray
        array of dilation event epochs to average and plot. 
        Should be dimensions (# of events x 2*half_epoch_duration + 1)
    random_epoch : numpy.ndarray
        array of random event epochs to average and plot. 
        Should be dimensions (# of events x 2*half_epoch_duration + 1)
    save_dir : str
        directory to save plot as a .png file

    Returns
    -------
    Plot saved as .png, if `save_dir` is provided

    Raises
    -------
    Assertion Error 
        epoch duration does not match the length of half epoch duration
    
    """
    
    time_vector = range(-half_epoch_duration, half_epoch_duration+1)
    plt.axhline(y=0, color='silver')
    plt.axvline(x=0, color='silver')
    plt.xlabel("Time (ms)")
    plt.ylabel("Pupil Size (pixels)")   
    if peak_epoch is not None: 
        mean_peak_epoch = np.mean(peak_epoch, axis=1)
        assert len(time_vector)==len(mean_peak_epoch), "Peak epoch duration does not equal data length"
        plt.plot(time_vector, mean_peak_epoch, color="crimson", label="Peak")
    if trough_epoch is not None: 
        mean_trough_epoch = np.mean(trough_epoch, axis=1)
        assert len(time_vector)==len(mean_trough_epoch), "Trough epoch duration does not equal data length"
        plt.plot(time_vector, mean_trough_epoch, color="navy", label="Trough")
    if constriction_epoch is not None:
        mean_constriction_epoch = np.mean(constriction_epoch, axis=1)
        assert len(time_vector)==len(mean_constriction_epoch), "Constriction epoch duration does not equal data length"
        plt.plot(time_vector, mean_constriction_epoch, color="deepskyblue", label="Constriction")
    if dilation_epoch is not None: 
        mean_dilation_epoch = np.mean(dilation_epoch, axis=1)
        assert len(time_vector)==len(mean_dilation_epoch), "Dilation epoch duration does not equal data length"
        plt.plot(time_vector, mean_dilation_epoch, color = "gold", label="Dilation")
    if random_epoch is not None: 
        mean_random_epoch = np.mean(random_epoch, axis=1)
        assert len(time_vector)==len(mean_random_epoch), "Random epoch duration does not equal data length"
        plt.plot(time_vector, mean_random_epoch, color="limegreen", label="Random")
    plt.legend()
    plt.title(title)

    if save_dir is not None: 
        plt.savefig(save_dir)
    plt.show()

class EventCollector(): 
    """
    Collects information about pupil phase events 

    Attributes
    --------
    idx : list 
        list of indexes in raw block data where event occured 
    times : list 
        list of times (ms) in raw block data where event occurred 
    pupil_size : list
        size of pupil when event occurred 
    diff_fit : list 
        list of gradient of search window 
    count : int
        tally of events that occur 
    type : string 
        what kind of event 
    accepted_event_idx : list 
        indexes of accepted events (i.e. events that occur after a pre-defined IEI)

    Methods
    --------
    update_data(self, idx, times, pupil_size, diff_fit=None):
        Update stored information when an event is found 
    store_accepted_event(self):
        Store index of an accepted event
    validate_epoch(self, pupil_data, time, half_epoch_duration): 
        Determine whether an event is acceptable for plotting 
    pull_single_epoch(self, pupil_data, time, half_epoch_duration):
        Pull a single epoch for plotting 
    pull_valid_epochs(self, pupil_data, half_epoch_duration, ms_per_sample, random_event=False):
        Pull all valid epochs for plotting 
    """

    def __init__(self, type): 
        """
        Sets up an EventCollector object to keep track of pupil phase events when running simulations
        
        Parameters
        ----------
        type : string 
            what kind of event the EventCollector is tracking 

        """
        self._idx = []
        self._times = []
        self._pupil_size = []
        self._diff_fit = []
        self._count = 0
        self._type = type
        self._accepted_event_idx = []

    def update_data(self, idx, time, pupil_size, diff_fit=None):
        """ 
        Update stored information when an event is found 

        Also updates the internal tally of the number of total events. Only meaningful to include
        diff_fit if not random event 

        Parameters
        ----------
        idx : int 
            index of event 
        time : int
            time (ms) of event 
        pupil_size : float 
            size of pupil 
        diff_fit : int 
            last value of the search window gradient 
        
        This function updates the following attributes: 
        ----------------------------------------------
        _idx : adds to the list of event indices 
        _times : adds to the list of times of events
        _pupil_size : adds to the list of pupil sizes
        _diff_fit : adds to the gradient fit 
        _count : increments tally of event

        """
        self._idx.append(idx-1) 
        self._times.append(time)
        self._pupil_size.append(pupil_size)
        if diff_fit is not None:
            self._diff_fit.append(diff_fit)
        self._count += 1

    def store_accepted_event(self): 
        """
        Store index of an accepted event, if identified as valid. 

        This function updates the following attributes: 
        ----------------------------------------------
        _accepted_event_idx : stores index of current event in list of accepted event indices

        """
        accepted_idx = self._idx[-1]
        self._accepted_event_idx.append(accepted_idx)
    
    def validate_epoch(self, pupil_data, time, half_epoch_duration): 
        """ 
        Determine whether an event is acceptable for plotting 

        Events considered valid if there is at least half_epoch_duration worth of time on 
        either side of event 

        Parameters
        ----------
        pupil_data : np.ndarray
            raw or preprocessed pupil size data 
        time : int 
            the time (in ms) of the event 
        half_epoch_duration : int 
            length of half epoch (ms) - how long to plot in before and after an event 
        
        Returns
        ----------
        boolean 
            True if event is valid for plotting
        """
        if time-half_epoch_duration >= 0 and time+half_epoch_duration <= len(pupil_data):
            return True 
        else: 
            return False
    
    def pull_single_epoch(self, pupil_data, time, half_epoch_duration):
        """
        Pull a single epoch for plotting 

        Parameters
        ----------
        pupil_data : np.ndarray 
            raw or preprocessed pupil size data 
        time : int
            time of the event 
        half_epoch_duration : int
            length of half epoch (ms) - how long to plot in before and after an event 
        
        Returns
        ----------
        demeaned_epoch_data : np.ndarray 
            demeaned epoch data, with dimensions (1, half_epoch_duration*2 + 1)

        """
        epoch_data = pupil_data[time-half_epoch_duration:time+half_epoch_duration+1] # add 1 to get equal number of events on each side (centered around event) 
        demeaned_epoch_data = epoch_data - np.nanmean(epoch_data)  
        return demeaned_epoch_data

    def pull_valid_epochs(self, pupil_data, half_epoch_duration, ms_per_sample):
        """
        Pull all valid epochs for plotting 

        pupil_data should have one timepoint per unit of half_epoch_duration (i.e. if half_epoch_duration
        is in ms, each value in pupil_data should be 1 ms apart).
        
        Parameters
        ----------        
        pupil_data : np.ndarray 
            raw or preprocessed pupil size data 
        half_epoch_duration : int
            length of half epoch (ms) - how long to plot in before and after an event 
        ms_per_sample : int 
            how many ms per pupil_sample 
        
        Returns
        ----------
        accepted_events : np.ndarray 
            array of the epochs surrounding each accepted event. Rows reflect unique events, columns 
            reflect time point 
        all-events : np.ndarray 
            array of epochs surrounding all identified events. Rows reflect unique events, columns 
            reflect time point 
                
        Notes 
        ---------- 
        If random event (i.e. if self._type == "random"), there are no accepted events, so 
        the returned accepted event array includes all events

        """ 

        all_events = [int(x*ms_per_sample) for x in self._idx]
        if self._type == "random":
           accepted_events = [int(x*ms_per_sample) for x in self._idx]
        else: 
            accepted_events = [int(x*ms_per_sample) for x in self._accepted_event_idx]

        all_epoch_data = None
        accepted_epoch_data = None

        for time in all_events:
            if self.validate_epoch(pupil_data, time, half_epoch_duration):
                demeaned_epoch_data = self.pull_single_epoch(pupil_data, time, half_epoch_duration) #already demeaned here 
                if all_epoch_data is None:
                    all_epoch_data = demeaned_epoch_data.reshape(((half_epoch_duration*2)+1,1))
                else: 
                    all_epoch_data = np.append(all_epoch_data, demeaned_epoch_data.reshape(((half_epoch_duration*2)+1,1)), axis=1)

                if time in accepted_events: 
                    if accepted_epoch_data is None: 
                        accepted_epoch_data = demeaned_epoch_data.reshape(((half_epoch_duration*2)+1,1))
                    else: 
                        accepted_epoch_data = np.append(accepted_epoch_data, demeaned_epoch_data.reshape(((half_epoch_duration*2)+1,1)), axis=1)


        return accepted_epoch_data, all_epoch_data
    
    def get_idx(self):
        return(self._idx)
    
    def get_times(self):
        return(self._times)
    
    def get_pupil_size(self):
        return(self._pupil_size)
    
    def get_diff_fit(self):
        return(self._diff_fit)
    
    def get_count(self):
        return(self._count)
    
    def get_accepted_event_idx(self): 
        return self._accepted_event_idx