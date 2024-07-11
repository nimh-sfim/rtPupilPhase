import os
import pylink
import sys
from psychopy import visual, logging, core, event

import numpy as np

def set_up_directories( behavioral_folder, eyelink_folder): 
    eyelinkFolder = (eyelink_folder + os.path.sep)

    # Check for the behavioral data directory, otherwise make it
    if not os.path.isdir(behavioral_folder):
            os.makedirs(behavioral_folder)  # If this fails (e.g. permissions) we will get error

    # Check for the EyeLink data directory, otherwise make it
    if not os.path.isdir(eyelink_folder):
            os.makedirs(eyelink_folder) # If this fails (e.g. permissions) we will get error
    # Specify code directory 
    _thisDir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(_thisDir)

def clear_screen():
    """ Clear window """ 
    
    win.fillColor = genv.getBackgroundColor()
    win.flip()

def terminate_task():
    """ Terminate the task gracefully and retrieve the EDF data file """
    
    el_tracker = pylink.getEYELINK()
    
    if el_tracker.isConnected():
        error = el_tracker.isRecording()
    
        if error == pylink.TRIAL_OK:
            # Stop EyeLink recording
            abort_trial()
    
        el_tracker.setOfflineMode()
        el_tracker.sendCommand('clear_screen 0')
        pylink.msecDelay(500)
        el_tracker.closeDataFile()         
        el_tracker.sendMessage('End EyeLink Recording')
        el_tracker.close()
    
    win.close()
    core.quit()
    sys.exit()
    
def abort_trial():   
    """Ends recording abruptly"""
    
    el_tracker = pylink.getEYELINK()
    
    # Stop EyeLink recording
    if el_tracker.isRecording():
        pylink.pumpDelay(100)
        el_tracker.stopRecording()  
    
    clear_screen(win)
    bgcolor_RGB = (116, 116, 116)
    el_tracker.sendMessage('!V CLEAR %d %d %d' % bgcolor_RGB)
    el_tracker.sendMessage('TRIAL_RESULT %d' % pylink.TRIAL_ERROR)
    
    return pylink.TRIAL_ERROR

def end_experiment() -> None:
    """End experiment"""
    el_tracker = pylink.getEYELINK()

    # Log
    logging.log(level=logging.EXP,msg='*** END EXPERIMENT ***')
    el_tracker.sendMessage('*** END EXPERIMENT ***')
    
    # Stop EyeLink recording
    pylink.pumpDelay(100)
    el_tracker.stopRecording()
    
    # Terminate task
    terminate_task()
    core.wait(2)
    win.close()
   
def quit_task(win) -> None:
    """Quit task based off of key presses"""
    
    # Record key presses
    allKeys = event.getKeys(['p','escape'])
    
    # Check keypresses
    if allKeys != None:
        
        # Loop over key presses
        for thisKey in allKeys:
            
            # If target key found, end experiment
            if thisKey in ['p', 'escape']:
                end_experiment(win)

def block_trigger(genv, win):
    """Display block trigger screen"""
    el_tracker = pylink.getEYELINK()

    # Log
    logging.log(level=logging.EXP,msg='Waiting for block trigger')
    el_tracker.sendMessage('Waiting for block trigger')
        
    # On-screen text
    onscreen_instructions = visual.TextStim(win, text='Waiting to start. Standby...', color = genv.getForegroundColor(), wrapWidth = scn_width/2) 
    onscreen_instructions.draw()
    win.flip()

    # Wait for trigger key press
    key = event.waitKeys(keyList=['5','t','escape', 'p'])

    # If pressed escape key, quit task
    if np.in1d(key,['escape','p']):
       end_experiment()
       
    # Log
    logging.log(level=logging.EXP,msg='Block trigger received')
    el_tracker.sendMessage('Block trigger received')

def instruction_continue():
    """Continue/proceed from instruction screen"""

    # Wait for key press
    key = event.waitKeys(keyList=['space', 'escape', 'p'])

    # End experiment
    if np.in1d(key,['escape','p']):
        end_experiment()
    
def instructions_screens(instruction:str): 
    """Function presents all the instructions needed for the task"""

    scn_width, scn_height = win.size
    # Setup text
    task_instructions = visual.TextStim(win, instruction, color = genv.getForegroundColor(), wrapWidth = scn_width/2, units = 'cm')
    
    # Clear window
    clear_screen(win)
    
    # Draw text screen
    task_instructions.draw()
    win.flip()   
    
    # Continue or quit from instructions
    instruction_continue()

def general_instruction_screens():
    """Present the task instructions"""
    
    # Log
    logging.log(level=logging.EXP,msg='General task instructions')

    # Define instructions text
    instructions = "Please fixate on the fixation point at all times."
    
    # Setup instructions
    task_instructions = visual.TextStim(win, text='', color='black', pos=[0, 5], units='cm')

    # Draw the text
    task_instructions.setText(instructions)
    task_instructions.draw()
    fixation.draw()
    win.flip()

    # Continue or quit from instructions
    instruction_continue()