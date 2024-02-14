# rtPupilPhase

This code runs real-time pupillometry script in Pyschopy.

In brief, this code streams in the pupillometry data as it is read by the eye-tracker and determines the phase of the pupil. For details about the algorithm used to determine pupil phase, please see [Kronemer et al., 2024](https://biorxiv.org/cgi/content/short/2024.02.12.579393v1).

## Requirements

- Python 3.6+
- [PsychoPy](https://www.psychopy.org/): version 2022.2.4
- if on a Mac or a Linux system, Eyelink Developers Kit is also necessary to work with PsychoPy - see the SR Research support site for more information.
  - running this script in PsychoPy will fail if you have not included `EyeLinkCoreGraphicsPsychoPy.py`, `error.wav`, `type.wav` and `qbeep.wav` in the directory with this script. These files can be found in any of the PsychoPy Coder examples that are provided with the Eyelink Developers Kit.

This script additionally requires the following Python packages:

- `pylink` : for connecting with the eye tracker
- `scipy` : for functions used in algorithm to find peaks
- `numpy` : general numerical processing functions

## Hardware Requirements/Notes

By default, this code was developed to work with an Eyelink 1000 Plus (SR Research, Inc.). If you are using a different Eyelink eye-tracker, please be sure to change the `eyelink_ver` variable (line 248). Given that the EyeLink 1000 Plus has a maximum online recording limit of 60Hz, pupil samples are taken every 17ms. This value may change depending on your hardware - if so, be sure to change the `ms_per_sample` (line 150) to reflect your hardware to ensure that the real-time algorithm performs as expected.

Additionally, if the primary display device is a macOS device with a built-in retina screen, or if you have an external monitor with the "Optimize for Built-in Retina Display" preference setting, ensure that the `use_retina` variable (line 305) is set to `True`.

## Real-time Parameters

Recommended parameters for the real-time processing of pupillometry data can be found in lines 165-198. If you would like to change these parameters, we suggest using the [MATLAB simulation procedures](https://github.com/nimh-sfim/rtPupilPhase/tree/main/simulations) provided in this repository to test the impact of these changes prior to online data collection.

## Task Parameters

The script provided is designed to include a fixation task with 10-min blocks, and a maximum of 10 blocks. In the manuscript, participants completed 5 blocks. If you would like to change these parameters, please adjust the `block_duration_sec` and `max_num_blocks` parameters (lines 135-138).

Alternative versions of this task may include the presentation of stimulus when a specific pupil phase is detected. This functionality may be implemented in the `accepted_pupil_event` function (line 559).

## Additional Eye-tracker settings

The default calibration used for this code is HV9 (i.e., 9 calibration points). By default, calibration sounds are turned off. If you would like to change this for any reason, see lines 237-365 for these settings.

## Outputs

Behavioral data and EyeLink data will be saved in the directories specified in lines 65-68.
