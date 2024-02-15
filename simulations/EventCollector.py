import numpy as np

class EventCollector(): 
    def __init__(self, type): 
        self._idx = []
        self._times = []
        self._pupil_size = []
        self._diff_fit = []
        self._count = 0
        self._type = type
        self._accepted_event_idx = []

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

    def update_data(self, idx, times, pupil_size, diff_fit=None):
        self._idx.append(idx-1) 
        self._times.append(times)
        self._pupil_size.append(pupil_size)
        if diff_fit is not None:
            self._diff_fit.append(diff_fit)
        self._count += 1

    def store_accepted_event(self): 
        accepted_idx = self._idx[-1]
        self._accepted_event_idx.append(accepted_idx)
    
    def get_accepted_event_idx(self): 
        return self._accepted_event_idx
    
    def validate_epoch(self, pupil_data, time, half_epoch_duration): 
        if time-half_epoch_duration >= 0 and time+half_epoch_duration <= len(pupil_data):
            return True 
        else: 
            return False
    
    def pull_single_epoch(self, pupil_data, time, half_epoch_duration):
        epoch_data = pupil_data[time-half_epoch_duration:time+half_epoch_duration+1] # add 1 to get equal number of events on each side (centered around event) 
        demeaned_epoch_data = epoch_data - np.nanmean(epoch_data)  
        return demeaned_epoch_data

    def pull_valid_epochs(self, pupil_data, half_epoch_duration, ms_per_sample, random_event=False): 

        all_events = [int(x*ms_per_sample) for x in self._idx]
        if random_event:
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