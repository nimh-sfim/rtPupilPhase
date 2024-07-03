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
        duration (ms) of half an epoch to plot (i.e. how much time before and after evenbt)
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