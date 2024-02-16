import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

def get_data(sub_dir, num_blocks = 5, recorded_eye=1): 
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
    for idx, item in enumerate(message_array):
        if item.size > 0:
            if item == match_string:
                break

    time = time_array[idx][0][0]
    return time

def pull_pupil_sample(data, pupil_sample_num, samples_in_pupil_sample): 
    if pupil_sample_num ==0: 
        current_pupil_sample = data[0:2, 0:int(samples_in_pupil_sample)]
    else: 
        current_pupil_sample = data[0:2, int(pupil_sample_num*samples_in_pupil_sample):int((pupil_sample_num+1)*samples_in_pupil_sample)]
    return current_pupil_sample

def plot_mean_timecourses(half_epoch_duration, title = "", peak_epoch=None, 
                          trough_epoch=None, constriction_epoch=None, dilation_epoch=None, 
                          random_epoch=None, save_dir = None):
    
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