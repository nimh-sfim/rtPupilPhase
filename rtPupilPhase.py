# ********************************************
# *** REAL TIME PUPILLOMETRY FIXATION TASK ***
# ********************************************

# This script can be run in PsychoPy to present a fixation task embedded with 
# the rtPupilPhase real time monitoring of pupil size fluctuations method.

# Written by: Sharif I. Kronemer, Tori Gobo and Catherine Walsh
# Last Modified: 7/23/2024

# ************************
# *** IMPORT LIBRARIES ***
# ************************

import os
import time
import argparse
import numpy as np
from psychopy import visual, core, logging
from psychopy.gui import DlgFromDict
from psychopy.data import getDateStr
from psychopy.event import getKeys
import pylink

import rtPupil_config 
from EyeLinkFunctions import validate_edf_fname, setup_eyelink, calibrate_eyelink
from PsychoPyFunctions import set_up_directories, instructions_screens, general_instruction_screens, block_trigger, end_experiment
from StimulusDecider import StimulusDecider

def main(block_length, max_num_blocks, baseline_duration_ms,
         max_search_window_duration_ms, num_random_events, IEI_duration_sec, 
         pupil_sample_duration_ms, peak_pupil_quantile, trough_pupil_quantile, 
         dilation_quantile, constriction_quantile):

    ## User Inputs ## 

    # Setup the subject info screen
    info = {'Session #': 1, 'Subject ID': 'Test', 'EyeLink': ['y','n'], 'EyeLink EDF': 'test.edf', '(1) Skip task instructions':['n', 'y']}
    dlg = DlgFromDict(info, title = 'Real Time Pupillometry Fixation Experiment')
    # Find experiment date
    info['date'] = getDateStr()
    # Filename = Subject ID entered above
    sub_id = info['Subject ID']
    # EyeLink EDF filename
    tmp_str = info['EyeLink EDF']

    # set up directories - if they don't exist, make them, and change directory to directory where this script is located
    behav_fname = os.path.join(rtPupil_config.data_fname, sub_id,'Behavior')
    eyelink_fname = os.path.join(rtPupil_config.data_fname, sub_id,'EyeLink')
    set_up_directories(behav_fname, eyelink_fname)

    # Show only critical log messages in the PsychoPy console
    logFile = logging.LogFile(behav_fname + os.path.sep + sub_id + '_Session_'+str(info['Session #'])+'_'+info['date']+'.log', level=logging.EXP)

    # log input parameters
    param_log_message = "Input Parameters: block length: "+str(block_length)+", max_num_blocks: "+str(max_num_blocks)+\
    ", baseline_duration_ms: "+str(baseline_duration_ms)+", max_search_window_duration_ms: "+str(max_search_window_duration_ms)+\
    ", num_random_events: "+str(num_random_events)+", IEI_duration_sec: "+str(IEI_duration_sec)+", pupil_sample_duration_ms: "+\
    str(pupil_sample_duration_ms)+", peak_pupil_quantile: "+str(peak_pupil_quantile)+", trough_pupil_quantile: "+str(trough_pupil_quantile)+\
    ", dilation_quantile: "+str(dilation_quantile)+", constriction_quantile: "+str(constriction_quantile)
    logging.log(level=logging.EXP,msg=param_log_message)

    # validate edf file name (length <= 8 & no special char)
    edf_fname, edf_state, edf_message = validate_edf_fname(tmp_str)
    if not edf_state:
        print(edf_message)
        core.quit()

    # We download EDF data file from the EyeLink Host PC to the local hard
    # drive at the end of each testing session, here we rename the EDF to
    # include session start date/time
    time_str = time.strftime("_%Y_%m_%d_%H_%M", time.localtime())
    session_identifier = edf_fname + time_str

    ### PsychoPy Items ### 

    # Basic window to use later 
    win = visual.Window(size = rtPupil_config.resolution, color = [0,0,0], monitor = 'testMonitor', fullscr = True, units ='cm')
    # Setup fixation cross
    fixation = visual.TextStim(win, text="+", color = rtPupil_config.text_color, pos = [0, 0], autoLog = False)
    block_timer = core.Clock() # Block timer
    general_timer = core.Clock() # Global timer

    ### Initiate Eyelink

    # EyeLink Dummy mode? - Set to False if testing with actual system
    if info['EyeLink'] == 'y':
        dummy_mode = False
        
    elif info['EyeLink'] == 'n':
        dummy_mode = True
        logging.log(level=logging.EXP,msg='Experiment run in dummy mode - no EyeLink')

    setup_eyelink(win, dummy_mode, edf_fname)
    calibrate_eyelink(win, dummy_mode)
    el_tracker = pylink.getEYELINK()
    
    # Set up StimulusDecider object for rtPupilPhase algorithm
    sd = StimulusDecider(block_length, baseline_duration_ms, 
                        max_search_window_duration_ms, pupil_sample_duration_ms, 
                        num_random_events, IEI_duration_sec, peak_pupil_quantile,
                        trough_pupil_quantile, dilation_quantile, constriction_quantile,
                        online=True, win=win) 
    
    ### Beginning of PsychoPy Experiment ### 

    # Define instruction text
    instructions_screens(win, "Experiment is setup!\n\nLet's get started!")
    
    # Continue with task instructions
    if 'n' == info['(1) Skip task instructions']:

        # Instructions screen
        general_instruction_screens(win, fixation)
    
    # Setup block counter
    block_counter = 1
    
    # Log
    logging.log(level=logging.EXP,msg='Starting Main Task Phase')
    el_tracker.sendMessage("Starting Main Task Phase")
    
    # Loop over task blocks
    while block_counter <= max_num_blocks:
        
        # Start block 
        block_instruction = 'Starting Block #' + str(block_counter) + "\n\nAre you ready?"
        instructions_screens(win, block_instruction)
        block_trigger(win)
        
        # Initialize variable to track events 
        decision_arr = []
        
        # Log
        logging.log(level=logging.EXP,msg='Starting Block ' + str(block_counter))
        el_tracker.sendMessage('Starting Block ' + str(block_counter))
        
        # Reset halfway logical 
        halfway_screen = False
        
        # Reset block timer to start block
        block_timer.reset()

        # Wait block duration
        while block_timer.getTime() < block_length:
            # Display a halfway completion screen
            if halfway_screen == False and block_timer.getTime() > block_length/2:
                halfway_screen = True

                # Show progress screen 
                progress_screen = visual.TextStim(win, text="You completed 50% of this block!", color=rtPupil_config.text_color)
                fixation.setAutoDraw(False)
                progress_screen.setAutoDraw(True)

                # Reset timer
                general_timer.reset()
                
                # Present block update for target duration
                while general_timer.getTime() < 2:
                    win.update()

                # Turn fixation back on
                progress_screen.setAutoDraw(False)
                fixation.setAutoDraw(True)
                win.update()
                
                # Reset timer
                general_timer.reset()
                
                # Timeout period before searching for pupil phase event
                while general_timer.getTime() < 2:
                    win.update()

            # Get all pressed keys. NOTE: allKeys will reset itself when getKeys is called
            allKeys = getKeys(['p','escape'])
    
            # If a key was pressed, check for quit keys
            if allKeys != None:
                for thisKey in allKeys:
                    if np.in1d(thisKey,['escape','p']):
                        end_experiment(win)
            
            # Turn on fixation
            fixation.setAutoDraw(True)
            win.update()
            if not dummy_mode:
                # Build search window
                sd.build_search_window()
                # Look for pupil phase events
                decision_arr.append(sd.detect_events_online())
        # Log
        logging.log(level=logging.EXP,msg='All pupil_sample Duration Times: ' + str(sd.get_pupil_sample_duration_time()))
        logging.log(level=logging.EXP,msg='All Search Window Detected Pupil Phase Events: ' + str(decision_arr))
        logging.log(level=logging.EXP,msg='Finished Block ' + str(block_counter))
        el_tracker.sendMessage('Finished Block ' + str(block_counter))
        
        # End block
        fixation.setAutoDraw(False)
        win.update()
        instructions_screens(win, 'Finished Block #' + str(block_counter))
        block_counter = block_counter + 1
        
    ### End experiment ###
    instructions_screens(win, "Exiting task! \n\nAn experimenter will communicate with you shortly.")
    end_experiment(win)
            
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'rtPupilPhase: Real-Time Pupillometry', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--max_num_blocks", help="Number of task blocks.", 
                        type=int, default=10)
    parser.add_argument("--block_length", help="Duration of a single block, in seconds.", 
                        type=int, default=600)
    parser.add_argument("--baseline_duration_ms", help="Duration of baseline window in milliseconds.", 
                        type=int, default=5000)
    parser.add_argument("--max_search_window_duration_ms", help="Maximum duration of search window before resetting, in milliseconds.",
                         type=int, default=5000)
    parser.add_argument("--num_random_events", help="Number of random events per block.", 
                        type=int, default=20)
    parser.add_argument("--IEI_duration_ms", help="Inter-event interval - how long to wait between valid events in seconds.", 
                        type=int, default=3000)
    parser.add_argument("--pupil_sample_duration_ms", help="How long we should consider a pupil sample in milliseconds.", 
                        type=int, default=100)
    parser.add_argument("--peak_pupil_quantile", help="Quantile value a peak must be bigger than to accept.",
                        type=float, default=0.75)
    parser.add_argument("--trough_pupil_quantile", help="Quantile value a trough must be smaller than to accept.",
                        type=float, default=0.25)
    parser.add_argument("--dilation_quantile", help="Quantile value a dilation must be bigger than to accept.",
                        type=float, default=0.99)
    parser.add_argument("--constriction_quantile", help="Quantile value a constriction must be smaller than to accept.",
                        type=float, default=0.01)
    
    args = parser.parse_args()

    main(args.block_length, 
        args.max_num_blocks, args.baseline_duration_ms, args.max_search_window_duration_ms,
        args.num_random_events, args.IEI_duration_ms/1000, args.pupil_sample_duration_ms, 
        args.peak_pupil_quantile, args.trough_pupil_quantile, args.dilation_quantile, 
        args.constriction_quantile)
