# rtPupilPhase

This is a repository for code associated with [Kronemer et al., 2024](https://www.biorxiv.org/content/10.1101/2024.02.12.579393v1). This repo includes three components:

1. PsychoPy code implementing a fixation task with real-time pupillometry. This task consists of a black fixation cross on a grey background. In this task, participants are asked to maintain their gaze on the central fixation point at all times and keep a steady head position in the head/chin rest system. Participants completed five 10-minute fixation task blocks. Between each block participants were given break until they indicated they were prepared to begin the next fixation period.
1. Human pupillometry data collected from the fixation task described above, and monkey pupillometry data.
1. Simulation code (in Python and MATLAB) to test different real-time pupillometry parameters. Although the parameters used in rtPupilPhase are robust for human pupillometry data, if you are recording data from other species, or you want to determine the effect of changing the parameters, you may wish to simulate the effects of the changes prior to collecting data.

## Getting Started

### Setting up the Environment

#### Python

Because SR-Research requires you to download their Developer's Toolkit in order to interface their eye tracker with PsychoPy, the process of getting this code is a little more complicated than just installing some Python packages. To that end, we have provided installation instructions for getting set up to use this code.

1. Ensure that you have the Eyelink Developer's ToolKit installed from SR-Research. This can be downloaded once you have made a free account with SR-Research.
1. Ensure that you have downloaded [PsychoPy](https://www.psychopy.org/). This code was developed with PsychoPy version 2022.2.4, but has also been tested with version 2024.1.5.
1. Clone this repo to your desired location (for more details on cloning a repo, see GitHub's instructions [here](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)).
1. Create and activate a new conda environment by opening a Terminal, move to the directory that you cloned this repo to and running:

    ```bash
    conda env create --name rtPupil --file environment.yml
    conda activate rtPupil
    ```

1. Install the appropriate version of `pylink` by running the `install_pylink.py` script from the EyeLink Developers ToolKit (see [their instructions](https://www.sr-research.com/support/thread-48.html) for more details). Note that the version on PyPi (that you might install via `pip`) is 0.3.3; this code was developed with `pylink` version 2.1.762.0.
1. Move `error.wav`, `type.wav` and `qbeep.wav` from one of the PsychoPy Coder examples from SR-Research (for example, those stored in `~/Applications/Eyelink/SampleExperiments/Python/examples/Psychopy_examples/Builder/EyeLinkMRIdemo_Builder`) into the cloned directory. These files are required for EyeLink calibration. We don't actually use these files (in fact, we actively turn off sounds in the calibration), but the `pylink` calibration code from SR-Research will crash if they don't exist.

The same Python environment can be used for the real-time PsychoPy code and the Python simulation code.

#### MATLAB

If you are running the MATLAB simulation code, these scripts additional require the following MATLAB toolboxes:

- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Curve Fitting Toolbox

Additionally, the code requires a `naninterp` function, which can be downloaded from the [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/8225-naninterp) and should be placed in the `utils` directory or otherwise added to the MATLAB path in an alternative manner.

#### Data Structure

We have provided pupillometry data for the human and monkey simulations (both data sets are previously unpublished; see manuscript for details). Mouse pupillometry data are [freely available online](https://www.sciencedirect.com/science/article/pii/S2211124723005387) and should be downloaded and placed in the `simulations/data/mouse` directory for the code to function properly. We have provided an example of what the data structure should look like for the simulation code to run properly; the real-time pupillometry script automatically creates this structure if it does not already exist. If you would like to change this data structure, you can change it in the `rtPupil_config.py` file (see below for details).

```bash
data
├── human
│   ├── 046
│   │   ├── Behavior
│   │   └── EyeLink
├── monkey
│   ├── Monkey_1
│   │   └── Monkey_1_pupil_data_230829.mat
│   └── Monkey_2
│       └── Monkey_2_pupil_data_210924.mat
└── mouse
    └── Mouse_pupillometry_data_sessions.mat

```

The simulation script will create directories for each subject in the following data structure:

```bash
analysis
└── subject_analysis
    ├── human
    ├── monkey
    └── mouse
```

Each simulation script will save images of the mean time courses for each pupil event and the underlying epoch data for those peaks. For the Python simulation code, it will look something like this: ![All traces](/static/all_trace.png)

## Running the code

### Real-time experiment

#### Command Line Interface

Once you have the environment set up, the code can be easily run from a command-line interface (such as Terminal on a Mac). For example, the real-time pupillometry task can be run using the default parameters from Kronemer et al., 2024 (listed in Table 1 of the manuscript) with the command:

```bash
python3 rtPupilPhase.py
```

You can also change the parameters of the rtPupilPhase algorithm easily from the command line. For example, to change the maximum length of the search window:

```bash
python3 rtPupilPhase.py --max_search_window_duration_ms 6000
```

You can see all of the options that you can specify by using the `-h` command:

```bash
python3 rtPupilPhase.py -h
```

When you run the script, you will get a startup screen from PsychoPy that will ask about a few details:

1. Skip task instructions (n = No, y = Yes): if you select "no", the participant will not be reminded to stay fixated on the central cross.
1. Eyelink (n = No, y = Yes): whether you are running the code with an EyeLink eye-tracker currently active. If you are intending on running the code in dummy mode (i.e. without an EyeLink), you must select "no", otherwise the code will crash.
1. EyeLink EDF: This defines the prefix to the EyeLink EDF file that will be created. This filename must be 8 characters or less (before the `.edf` extension) and only contain basic numbers and letters.
1. Session #: Experimenter can define any number to be logged with the behavioral file along with the timestamp of running the task.
1. Subject ID: any character/numeric value that the user specifies to define the log file.

When running the task, most screens are advanced by pressing the space bar to continue. The only exception to this is the `Waiting to start` screen, which will only advance by pressing the `5` or `t` key.
You can quit the task at any time by pressing `p` or the `escape` key. If the task refuses to quit, you can force quit by using the keystroke combination `option + command + escape`.

#### PsychoPy GUI

For the easiest use and most flexibility, we recommend running the fixation experiment from the command line. However, we have also provided two options to implement the rtPupilPhase algorithm in the PsychoPy GUI.

First, you can run the real-time script (i.e. `rtPupilPhase.py`) directly from the PsychoPy GUI. If you prefer this method and wish to change the real-time parameters, you can directly change the inputs to the `main` function at the end of the `rtPupilPhase.py` script.

We have additionally provided a Builder implementation of the fixation task (`rtPupilPhase_builder.psyexp`). This basic code includes setup for the EyeLink eye-tracker and a basic fixation task set up with the default parameters from Kronemer et al., 2024. To change the `rtPupilPhase` algorithm parameters, please edit the code in the `Before Experiment` tab of the `StimulusDecider` code block in the `fixation_routine` routine. A basic loop for a block structure is provided; you can edit this in `BlockData.csv`. We hope that providing a Builder example will lower the barrier to entry for users who have less experience writing PsychoPy code from scratch so that rtPupilPhase can be used in new, exciting kinds of tasks.

### Simulations

Similarly, the simulated rtPupilPhase can be run from the command line with the command:

```bash
python3 simrtPupilPhase_human.py 046 048
```

There are no required parameters for the real-time pupillometry code; the simulation Python code requires a list of subject IDs to run through (see example above for syntax). Just as with the real-time code, you can specify the parameters of the real-time algorithm using command line options. To see all of the options, you can run

```bash
python3 simrtPupilPhase_human.py -h
```

Note that the simulation code assumes that you have already converted EDF files to `.mat` files. A tool acccomplish this is provided through the EyeLink Developers ToolKit - please see [this thread](https://www.sr-research.com/support/thread-54.html) in the SR-Research forum for more information.

The MATLAB scripts can be run out of the box in MATLAB.

### Known issues with running the code

- If you run the PsychoPy code intending it to be in dummy mode (i.e. no eye-tracker connected) and tell it there is an EyeTracker, the code will try to connect to the non-existant tracker, hang there and eventually crash.
- This code is written so it expects the EyeLink to be tracking the right eye. If the left eye is selected on the EyeLink, the real-time pupillometry will crash. If you need to collect data from the left eye, you can adjust the method used in the `build_search_window` function from the `StimulusDecider` class.
- The Python simulation code does not include steps for identifying and processing blinks/microsaccades. You may wish to apply additional cleaning to the data.

### Optional Parameters

Although the code is set up to use the default parameter values from the human data analysis in Kronemer et al., 2024, the scripts are set up to easily take in parameters as command line options, as detailed below.

| Parameter | Description | Default | Relevant Script |
| ---| ---| --- | --- |
| `max_num_blocks` | Maximum number of task blocks to run through | 10 | `rtPupilPhase.py` |
| `block_length` | Duration of block, in seconds | 600 | `rtPupilPhase.py` |
| `baseline_duration_ms` | Duration of baseline window in milliseconds | 5000 |  `rtPupilPhase.py` and `simrtPupilPhase_human.py`|
| `max_search_window_duration_ms` | Maximum duration of search window before resetting, in milliseconds| 5000 | `rtPupilPhase.py` and `simrtPupilPhase_human.py`|
| `num_random_events` | Number of random events per block. Note that in the simulations, this parameter (along with the block length) will impact the number of pupil phase events, as the code will force this number of events and will only accept events after the inter-event interval. | 20 | `rtPupilPhase.py` and `simrtPupilPhase_human.py` |
| `IEI_duration_ms` | Inter-event interval - how long to wait between valid events, in milliseconds| 3000 |  `rtPupilPhase.py` and `simrtPupilPhase_human.py` |
| `pupil_sample_duration_ms` | Pupil sample length in milliseconds | 100 | `rtPupilPhase.py` and `simrtPupilPhase_human.py`|
| `peak_pupil_quantile` | Quantile value a peak must be bigger than to accept | 0.75 |  `rtPupilPhase.py` and `simrtPupilPhase_human.py` |
| `trough_pupil_quantile` | Quantile value a trough must be smaller than to accept.| 0.25 | `rtPupilPhase.py` and `simrtPupilPhase_human.py` |
| `dilation_quantile` | Quantile value a dilation must be bigger than to accept.| 0.99 | `rtPupilPhase.py` and `simrtPupilPhase_human.py` |
| `constriction_quantile` | Quantile value a constriction must be smaller than to accept. | 0.01 | `rtPupilPhase.py` and `simrtPupilPhase_human.py` |
| `plot_timecourses` | Flag determining whether or not to plot mean timecourses in simulation. | Flag optional |  `simrtPupilPhase_human.py` |
| `half_epoch_duration_ms` |  How far before and after event to plot | 2500 | `simrtPupilPhase_human.py` |

### Additional Configuration

Additional configuration options are provided in the `rtPupil_config.py` file:

- Real-time parameters - Initial thresholds: these are the initial thresholds for pupil values from the real-time algorithm. They get updated whenever there is a new pupil sample.
- EyeLink information: This is information about the EyeLink eye tracker used. This code was developed for an EyeLink 1000 Plus, with 60Hz online sampling and 1000Hz offline sampling and tracking the right eye.
- Display monitor information: This is information about the monitor that the task is being run on, which may impact the PsychoPy task and some blink/microsaccade detection. Note that it also includes the `use_retina` variable, which should be set to `True` if you use the built-in retina screen as your primary display device on MacOS, or if you choose to "Optimize for Built-in Retina Display" for an external monitor.
- PsychoPy information: Settings for the color of the background and text of the PsychoPy experiment (default to grey and black, respectively)
- Filenames: This sets up the default data structure that the Python simulations expect to see the data in and where the real-time code saves the data to.
- Task structure (for simluations): number of blocks and block duration for the simulations.

### Customizing the PsychoPy script

The scripts provided are built around a fixation task: participants were completed 10 5-minute blocks of a black fixation cross on a grey screen. The PsychoPy script included in this repo will run this task. If you would like to alter this script to create a closed-loop paradigm such that an event happens upon detection of a pupil phase event, we recommend inserting code in the `accepted_pupil_event()` method from the `StimulusDecider` class.

## API

The main workhorse of the rtPupilPhase algorithm is contained in the `StimulusDecider` class - the `StimulusDecider` object is used in both the real-time PsychoPy code and the offline Python simulation code. The Python simulation code additionally requires the `EventCollector` module to keep track of simulated events and manipulate the offline data. We provide two additional modules (`EyeLinkFunctions.py` and `PsychoPyFunctions.py`) that support the real-time PsychoPy code.

### `StimulusDecider` Module

An object of class `StimulusDecider` stores data and detects pupil phase events in real time and in simulations.

Please see Kronemer et al., 2024 for more detailed descriptions of algorithm parameters.

#### Properties defined at initialization that never change

|Property | Type| Description|
| ---- | ----| ----|
| `online`| boolean | Whether the object is used in real-time data collection or in simulations |
| `baseline_duration_ms` | int | Duration of baseline window in milliseconds |
| `max_search_window_duration_ms` | int | Maximum length of search window (in milliseconds) before resetting |
| `pupil_sample_duration_ms` | int | Duration of a single pupil sample, in milliseconds |
| `num_random_events` | int | How many random events per block |
| `random_event_time_sec` | int | How long between random events during real-time, in seconds |
| `IEI_duration_sec` | int | Inter-event interval (i.e. how long should we wait between pupil events), in seconds |
| `peak_pupil_quantile` | float | Quantile threshold for identifying peak values - pupil size must be above this quantile of values in baseline window to be accepted as a peak  |
| `trough_pupil_quantile` | float | Quantile threshold for identifying trough values - pupil size must be below this quantile of values in baseline window to be accepted as a trough  |
| `dilation_quantile` | float | Quantile threshold for identifying dilations - used to define absolute value threshold  |
| `constriction_quantile` | float | Quantile threshold for identifying constrictions - used to define absoulte value threhold  |

#### rtPupil Algorithm Thresholds

|Property | Type| Description|
| ---- | ----| ----|
| `peak_threshold_var` | float | Absolute threshold for identifying a peak event - event must be larger than this value to be identified as a peak. Taken from `rtPupil_config.py` file.  |
| `trough_threshold_var` | float | Absolute threshold for identifying a trough event - event must be smaller than this value to be identified as a trough. Taken from `rtPupil_config.py` file.  |
| `constriction_threshold` | float | Absolute threshold for identifying a constriction event - first derivative of search window must be smaller than this value to be identified as a constriction. Taken from `rtPupil_config.py` file.  |
| `dilation_threshold` | float | Absolute threshold for identifying a dilation event - first derivative of search window must be larger than this value to be identified as a dilation. Taken from `rtPupil_config.py` file.   |

#### Search window attributes that are updated internally

|Property | Type| Description|
| ---- | ----| ----|
| `baseline_window` | list | Set of pupil samples to be used to determine thresholds for finding pupil samples |
| `search_window` | list | Set of pupil samples to be used to find pupil events |
| `search_window_sample_times` | list | List of timestamps (ms) associated with search window |
| `prior_search_window` | list | Previous search window. For online algorithm, must have no blinks for current pupil sample to be valid. |
| `search_window_model_fits` | list | Final values from the model fits on pupil samples |
| `pupil_sample_duration_time` | list | Durations of pupil samples collected in real time  |

#### Attributes to internally track events

|Property | Type| Description|
| ---- | ----| ----|
| `new_sample` | obj | Most recent EyeLink sample (collected in real time)  |
| `old_sample` | obj | Previous EyeLink sample (collected in real time)  |
| `peak_count` | int | Number of peaks identified in real-time  |
| `trough_count` | int | Number of troughs identified in real-time  |
| `dilation_count` | int | Number of dilations identified in real-time  |
| `constriction_count` | int | Number of constrictions identified in real-time  |
| `idx_event` | int | Flag for event that was identified. 1 = peak, -1 = trough, 2 = dilation, -2 = constriction, 0 = no event |
| `accepted_pupil_event_bool` | boolean | Internal variable for marking whether an event was accepted |
| `pupil_phase_IEI_timer` | PsychoPy core.Clock() object | Timer to keep track of inter-event interval in real-time |
| `random_IEI_timer` | PsychoPy core.Clock() object | Timer to keep track of random events in real-time |
| `pupil_sample_IEI_timer` | PsychoPy core.Clock() object | Timer to keep track of duration of pupil sample in real-time |
| `win` | PsychoPy Window | PsychoPy window used in real-time experiment |

#### Methods

The following methods are available for an instance of the `StimulusDecider` class.

| Call | Inputs | Outputs | Description |
| ----| ----| ----| --- |
| `init()` | _1._ `block_duration_sec`: Number of seconds in a block. Default: 600 _2._ `baseline_duration_ms`: Duration of baseline window in ms. Default: 5000 _3._ `max_search_window_duration_ms`: Maximum length of search window in ms. Default: 5000 _4._ `pupil_sample_duration_ms`: duration of real-time pupil sample from eye tracker in ms. Default: 100 _5._ `num_random_events`: number of random events to include per block. Default: 20 _6._ `IEI_duration_sec`: duration of inter-event interval in seconds. Default: 3 _7._ `peak_pupil_quantile`: quantile threshold for identifying peak values. Default: 0.75 _8._ `trough_pupil_quantile`: quantile threshold for identifying trough values. Default: 0.25 _9._ `dilation_quantile`: quantile threshold for identifying dilations. Default: 0.99 _10._ `constriction_quantile`: quantile threshold for identifying constrictions. Default: 0.01 _11._ `online`: whether object is used in real-time or simulations. Default: False _12._ `win`: screen for PsychoPy, if real-time. Default: None | | Initialize `StimulusDecider` object with defaults that reflect options from Kronemer et al., 2024 |
| `build_search_window()` | | | Build baseline and search window in real-time |
| `detect_events_online()` | | Integer reflecting kind of event detected. 1 = peak, -1 = trough, 2 = dilation, -2 = constriction, 3 = random event, 0 = no event | Detect pupil events in real-time. Validate search window, update pupil phase thresholds (if necessary) and identify pupil phase event in search window.  |
| `validate_search_window()` | _1._ `max_search_window_duration_samples`: maximum number of pupil samples in each search window | boolean: `True` if valid search window, `False` otherwise | Determine whether search window meets criteria for finding pupil phase events |
| `update_pupil_phase_thresholds()` | | | Update thresholds for determining pupil phase events according to current baseline window. Only updates internal variables. |
| `find_pupil_phase_event()` | _1._ `pupil_sample_num`: index of pupil sample. Default: `np.nan` _2._ `samples_in-pupil_sample`: how many size values in pupil sample. Default: 6. _3._ `current_time`: time (ms) of pupil sample. Default: np.nan. _4._ `peak_events`: EventCollector object (or `None`) to log peak events. Default: None. _5._ `trough_events`: EventCollector object (or `None`) to log trough events. Default: None. _6._ `constriction_events`: EventCollector object (or `None`) to log constriction events. Default: None. _7._ `dilation_events`: EventCollector object (or `None`) to log dilation events. Default: None | Integer reflecting kind of event detected. 1 = peak, -1 = trough, 2 = dilation, -2 = constriction, 0 = no event. | Find possible pupil phase events and log them |
| `fit_polynomial()` | _1._ `demeaned_pupil_sample`: sample to fit polynomial on. Should already be demeaned | | Do polynomial fit on pupil sample and save final value |
| `accepted_pupil_event()` | | | Log an accepted pupil event. This function can be used to build closed-loop paradigms where a detected pupil event triggers some sort of other task event. |
| `update_windows()` | _1._ `sample`: Pupil sample to be added _2._ `duration`: duration of pupil sample, if collected in real-time. Default: None | | Update search and baseline windows with values from pupil sample |
| `log_found_event_live()` | _1._ `kind`: kind of event to log. _2._ `demeaned_search_window`: search window that is being fit 3. `diff_fit`: gradient of search window being fit | | Interface with eye-tracker to log an event in real-time |
| `validate_event_offline()` | _1._ `all_event_times`: list of all event times. _2._ `accepted_pupil_event_times`: list of accepted pupil event times. _3._ `IEI_jitter_ms`: amount of time that must have passed for event to be accepted._4._ `peak_events`: EventCollector object to log peak events. _5._ `trough_events`: EventCollector object to log trough events. _6._ `constriction_events`: EventCollector object to log constriction events. _7._ `dilation_events`: EventCollector object to dilation events. | _1._ `peak_events`: EventCollector object with updated event info. _2._ `trough_events`: EventCollector object with updated event info. _3._ `constriction_events`: EventCollector object with updated event info. _4._ `dilation_events`: EventCollector object with updated event info. | Determine whether enough time has passed to accept an identified pupil event in simulation |
| `reset_baseline_window()` | | | Clear baseline window |
| `reset_search_window()` | | | Clear search window and related variables |
| `get_pupil_sample_duration_time()` | | Integer duration of pupil sample duration | Used to log pupil sample duration array |
| `get_search_window()` | | List containing current search window | Accesses pupil samples in currently stored search window |
| `get_search_window_times()` | | List containing current search window times | Accesses times associated wiht currently stored pupil sample, as the EyeLink takes samples stochastically |
| `get_baseline_window()` | | List containing baseline window | Accesses current baseline window |
| `get_search_window_fit_vals()` | | List containing current search window fit values. | Accesses fitted search window values |
| `set_current_event()` | _1._ `found_event`: integer coding found event | | Update internal event tracker |

### `EventCollector` Module

#### Module level functions

The below methods can be called directly. These functions mostly relate to formatting eye tracking data so it can be used with a `StimulusDecider` object.

| Call | Inputs | Outputs | Description |
| ----| ----| ----| --- |
| `get_data()`| _1._ `recorded_eye`: Which eye was used to track. Default: 1 (reflecting right eye). | _1._ `all_data`: array of eye tracking data, of dimensions (4, number of time points). Rows reflect time (ms), pupil size (pixels), gaze X location, gaze Y location. _2._ `event_data`: array of messages from EyeLink events, with dimensions (2, number of events). Rows reflect message sent to EyeLink and time (ms) associated with message | Find and load eye-tracking data |
| `find_string_time()` | _1._ `time_array`:  array of timestamps associated with messages. _2._ `message_array`: array of messages to search against. _3._ `match_string`: string to find. | _1._ Time of event matching input string. | Find the time (ms) of a given event. |
| `pull_pupil_sample()` | _1._ `data`: numpy.ndarray of downsampled eyetracking data from block. First 2 rows should be time and pupil size data. _2._ `pupil_sample_num`: integer reflecting sample to pull. _3._ `samples_in_pupil_sample`: how many pupil size readings are in each pupil_sample. | numpy.ndarray of extracted pupil sample | Pull a specific pupil sample from gaze data |
| `plot_mean_timecourses()` | _1._ `half_epoch_duration`: duration (ms) of half an epoch to plot (so both before and after event is plotted). _2._ `title`: title of plot. Default: "". _3._ `peak_epoch`: numpy.ndarray of peak epochs to average and plot. Default = None _4._ `trough_epoch`: numpy.ndarray of trough epochs to average and plot. Default: None. _5._  `constriction_epoch`: numpy.ndarray of constrcition epochs to average and plot. Default: None. _6._ `dilation_epoch`: numpy.ndarray of dilation epochs to average and plot. Default: None. _7._  `random_epoch`: numpy.ndarray of random epochs to average and plot. Default: None. _8._ `save_dir`: directory to save plot as a .png file. Default: None. | Plot saved as a .png file, if `save_dir` is provided | Plot the mean epoch timecourses. Will plot as many event types as is provided, but must provide them as labeled arguments to ensure the color legend is correct. |

#### Properties

The following are properties of an `EventCollector` object that will keep track of pupil phase events.

| Property| Type | Description|
| ---- | ---- | ---- |
| `idx` | list | list of indexes in raw block data where event occurred |
| `times` | list | list of times (ms) in raw block data where event occurred |
| `pupil_size` | list | list of size of pupil when event occurred |
| `diff_fit` | list | gradient of the search window |
| `count` | int | tally of events that occur |
| `type` | string | what kind of event |
| `accepted_event_idx` | list | Indexes of accepted events (i.e. events that occur after a pre-defined IEI) |

#### Methods

The following methods are available for an instance of class `EventCollector`.

| Call | Inputs | Outputs | Description |
| ----| ----| ----| --- |
| `init()` | _1._ `type`: what kind of event the EventCollector is tracking | | Set up an EventCollector object to keep track of pupil phase events when runnign simulations|
| `update_data()` | _1._ `idx`: integer index of event. _2._ `time`: integer time (ms) of event. _3._ `pupil_size` size of pupil. _4._ `diff_fit` last value of search window gradient. | | Update stored information when an event is found |
| `store_accepted_events()` | | | Store index of an accepted event|
| `validate_epoch()` | _1._ `pupil_data`: numpy.ndarray of raw or preprocessed pupil size data. _2._ `time` time (ms) of event. _3._ `half_epoch_duration`: length of half epoch (ms) | boolean: `True` if event is valid for plotting | Determine whether an event is acceptable for plotting. Events considered valid if there is at least half_epoch_duration worth of time on either side of event. |
| `pull_single_epoch()` | _1._ `pupil_data`: numpy.ndarray of raw or preprocessed pupil size data. _2._ `time`: time of event. _3._ `half_epoch_duration`: length of half epoch (ms) to plot before and after event. | `demeaned_epoch_data`: np.ndarray of demeaned epoch data, with dimensions (1, half_epoch_duration*2 + 1) | Pull a single epoch for plotting from offline data |
| `pull_valid_epochs()` | _1._ `pupil_data`: numpy.ndarray of raw or preprocessed pupil size data. _2._ `half_epoch_duration`: length of half epoch (ms) to plot before and after event. _3._ `ms_per_sample`: how many ms per pupil sample | _1._ `accepted_events`: numpy.ndarray of epochs surrounding each accepted events. Rows reflect unique events, columns reflect time point. _2._ `all_events`: numpy.ndarray of epochs surrounding all events. Rows reflect unique events, columns reflect time point. | Pull all valid epochs for plotting. Note that if random event, there are no accepted events so the accepted event array includes all events. |

### `EyeLinkFunctions` Module

These are functions for interfacing with an SR Research, Inc. eye tracker. Many of these functions are based on functions provided by SR Research, Inc. Note that some of these scripts depend on settings defined in `rtPupil_config.py`.

#### Module level functions

| Call | Inputs | Outputs | Description |
| ----| ----| ----| --- |
| `validate_edf_fname()` | _1._ `edf_fname`: file name to validate | _1._ `edf_fname`: valid filename _2._ `state`: boolean, `True` if valid filename. _3._ `message`: message indicating why filename is invalid. |  Determine whether provided EDF filename is valid. |
| `setup_eyelink()` | _1._ `win`: PsychoPy screen _2._ `dummy_mode`: whether eye tracker is to be used _3._ `edf_fname`: filename to save EDF data. | | Set up EyeLink eye-tracker |
| `calibrate_eyelink()` | _1._ `win`: PsychoPy screen _2._ `dummy_mode`: whether eye tracker is to be used | | Calibrate EyeLink eye tracker |

### `PsychoPyFunctions` Module

These are functions relating to the PsychoPy implementation of the real-time experiment.

#### Module level functions

| Call | Inputs | Outputs | Description |
| ----| ----| ----| --- |
| `set_up_directories()` | _1._ `behavioral_folder`: where to save behavioral log files. _2._ `eyelink_folder`: where to save EyeLink EDF files | | Set up directories to save real-time data. If directories do note exist, make them. |
| `clear_screen()` | _1._ `win`: PsychoPy screen for experiment | | Clear window |
| `terminate_task()` | _1._ `win`: PsychoPy screen for experiment | | Terminate the task gracefully and retrieve the EDF data file |
| `abort_trial()` | _1._ `win`: PsychoPy screen for experiment | | Ends recording abruptly |
| `end_experiment()` | _1._ `win`: PsychoPy screen for experiment | | End experiment, stop EyeLink recording |
| `quit_task()` | _1._ `win`: PsychoPy screen for experiment | | Quit task based off of key presses |
| `block_trigger()` | _1._ `win`: PsychoPy screen for experiment | | Display the block trigger screen |
| `instruction_continue()` | _1._ `win`: PsychoPy screen for experiment | | Continue from instruction screen |
| `instructions_screens()` | _1._ `win`: PsychoPy screen for experiment _2._ `instruction`: string to be presented | | Presents instructions needed for task |
| `general_instruction_screens()` | _1._ `win`: PsychoPy screen for experiment _2._ `fixation`: PsychoPy TextStim for fixation | | Present the general task instructions |
