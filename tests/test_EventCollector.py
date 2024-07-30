import pytest
import numpy as np
from simulations.EventCollector import EventCollector

@pytest.mark.parametrize("input_test, expected", [([1,1,1], [1,1,1,0] ), ([1,1,1,1], [1,1,1,1])])
def test_update_data(input_test, expected):
    e = EventCollector("test")
    e.update_data(*input_test)
    assert len(e._idx)==expected[0]
    assert len(e._times)==expected[1]
    assert len(e._pupil_size)==expected[2]
    assert len(e._diff_fit)==expected[3]
    assert e._count == 1
    
def test_store_accepted_event(): 
    e=EventCollector("test")
    e.update_data(1,1,1)
    e.store_accepted_event()
    assert e._accepted_event_idx == [0] # check it adds when empty
    e.update_data(2,2,2)
    e.store_accepted_event()
    assert e._accepted_event_idx == [0,1] # check it adds when there's something in it 

@pytest.mark.parametrize("time, expected",[ (5,True), (2,False), (9,False)])
def test_validate_epoch( time, expected):
    pupil_data = [1,2,3,4,5,6,7,8,9,10]
    half_epoch_duration = 5
    e= EventCollector("test")
    valid = e.validate_epoch(pupil_data, time, half_epoch_duration)
    assert valid == expected

def test_pull_single_epoch():
    pupil_data = [1,2,3,4,5,6,7,8,9,10]
    time = 5
    half_epoch_dur = 2
    e = EventCollector("test")
    epoch_data = e.pull_single_epoch(pupil_data, time, half_epoch_dur) 
    assert len(epoch_data) == 2*half_epoch_dur+1 # pulling the correct number of values
    assert np.mean(epoch_data)==0 # demeaning is working correctly
    assert np.array_equal(epoch_data , [-2,-1,0,1,2]) # expected output

def test_pull_valid_epochs():
    pupil_data = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
    half_epoch_duration = 2 
    ms_per_sample = 1 
    e = EventCollector("test")
    e.update_data(4,1,1,1)
    e.store_accepted_event()
    e.update_data(5,1,1,1)
    e.update_data(10,1,1,1)
    e.store_accepted_event()

    # added 3 events, 2 were accepted 
    # epoch should be 5 - 2 * half epoch duration + 1 


    accepted_events, all_events = e.pull_valid_epochs(pupil_data, half_epoch_duration, ms_per_sample)
    assert np.shape(accepted_events) == (5,2)
    assert np.shape(all_events) == (5,3)

    # data should be demeaned 
    assert np.mean(accepted_events[:,0]) == 0

def test_pull_valid_epochs_random(): 
    pupil_data = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
    half_epoch_duration = 2 
    ms_per_sample = 1 
    e = EventCollector("random")
    e.update_data(4,1,1,1)
    e.update_data(5,1,1,1)
    e.update_data(10,1,1,1)
    # with random events, no "accepted" events so should be same size 

    # added 3 events
    # epoch should be 5 - 2 * half epoch duration + 1 

    accepted_events, all_events = e.pull_valid_epochs(pupil_data, half_epoch_duration, ms_per_sample)
    assert np.shape(accepted_events) == (5,3)
    assert np.shape(all_events) == (5,3)

    # data should be demeaned 
    assert np.mean(accepted_events[:,0]) == 0


def test_get_idx(): 
    e = EventCollector("test")
    e.update_data(1,1,1,1)
    idx = e.get_idx()
    assert idx == [0]

def test_get_times():     
    e = EventCollector("test")
    e.update_data(1,1,1,1)
    times = e.get_times()
    assert times == [1]

def test_get_diff_fit():     
    e = EventCollector("test")
    e.update_data(1,1,1,1)
    diff_fit = e.get_diff_fit()
    assert diff_fit == [1]

def test_get_count():     
    e = EventCollector("test")
    e.update_data(1,1,1,1)
    count = e.get_count()
    assert count==1

def test_get_accepted_idx():     
    e = EventCollector("test")
    e.update_data(1,1,1,1)
    e.store_accepted_event()
    accepted_idx = e.get_accepted_event_idx()
    assert accepted_idx == [0]