# Simulated rtPupilPhase

The MATLAB scripts in this directory (developed in MATLAB version 2022b) is designed to simulate the real-time pupillometry in rtPupil. The goal of this code is to take pre-existing pupillometry data and simulate the detection of events, as though it were collected in real-time. This process allow a user to test the impact of different parameter choices or task designs prior to online data collection.

This directory includes scripts that can be used to simulate human, macaque and rat eye-tracking data. The scripts identify pupil phase events (peak, trough, constriction and dilation), in addition to blinks, saccades and microsaccades.

## Requirements

- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Curve Fitting Toolbox

Additionally, these scripts require helper functions, which are provided in the `utils` directory. Please ensure that these files are available, or that you change the paths to in the scripts as necessary.

## Data

We have provided pupillometry data for the human and monkey simulations. Mouse pupillometry data are [freely available online](https://www.sciencedirect.com/science/article/pii/S2211124723005387) and should be downloaded and placed in the `simulations/data/mouse` directory for the code to function properly.

Some parameters associated with these scripts are unique to the [experimental task and pupillometry acquisition](link) used in Kronemer et al., 2024. Therefore, updates to these parts of the script may be required when testing alternative data sets. Depending on your specific hardware (i.e. eye-tracker and display monitor), you may need to adjust values in the code to ensure correct calculations. Please ensure that you update the appropriate variables in the `Display Parameters` section of the scripts. Also please note that all data should be converted to `.mat` files prior to running these simulations - if data were collected using an Eyelink eye-tracker, the Eyelink Developers Toolkit provides an application to convert `.edf` files to `.mat` files.

## Real-time parameters

Although the details about the data across species may differ (see below for species specific code information), the underlying algorithm remains the same. By default, the parameters were developed using an Eyelink 1000 Plus. Of potential interest are:

- `sampling_rate`: overall sampling rate for eye-tracker
- `ms_per_sample`: online sampling rate for eye-tracker
- `downsample_value`: the Eyelink 1000 Plus ultimately samples at 1000Hz, but only has an online sampling rate of 60Hz. This value allows for the offline data to be downsampled to be effectively the same rate as though data were being streamed in real time.
- `pupil_sample_duration_ms`: how long to sample for while looking for events
- `num_random_events`: how many random events (i.e. events unrelated to a pupil event) should be given per block
- `baseline_ms`: duration of pupil baseline period for setting new thresholds
- `peak_threshold`, `trough_threshold`, `rising_threshold`, `falling_threshold`, `peak_pupil_quantile`, `trough_pupil_quantile`, `rising_quantile`, `falling_quantile`: thresholds for determining an event during a given sample window, although note that these will be updated automatically as new data is collected. Please see [Kronemer et al., 2024](link) for details on how these values are used.
- `IEI_jitter_ms`: how long to wait before recording sequential accepted events
- `max_search_window_length_ms` how long to search for events before discarding data
- `half_epoch_duration_ms`: how long to visualize before and after time-locked event

## Outputs

The following outputs will be saved for each subject:

- a text file containing the total number and number of accepted identified events (dilations, peaks, constrictions, troughs and random).
- figures plotting all epochs for each kind of pupil event, mean time course of each kind of pupil event, and the mean timecourse surrounding blinks.
