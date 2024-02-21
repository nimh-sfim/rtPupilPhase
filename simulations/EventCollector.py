import numpy as np

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
            
        """
        self._idx.append(idx-1) 
        self._times.append(time)
        self._pupil_size.append(pupil_size)
        if diff_fit is not None:
            self._diff_fit.append(diff_fit)
        self._count += 1

    def store_accepted_event(self): 
        """
        Store index of an accepted event

        Method called if event is identified as valid. Pull the final value from the list of all
        event indexes and add it to the internal tracker self._accepted_event_idx
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
            demeaned epoch data, with dimensions (1 x half_epoch_duration*2 + 1)

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