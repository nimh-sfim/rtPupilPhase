# rtPupilPhase

This is a repository for code associated with [Kronemer et al., 2024](link).

This repository includes two main parts:

- Psychopy code associated with collecting real-time pupillometry data. This includes a fixation task where pupil events (peak, trough, dilation, constriction) are identified in real-time and recorded. This code could be extended to include the presentation of a stimulus upon the detection of a given pupil event, for example.
- Simulated real-time pupillometry. This includes a MATLAB script that reads in pupillometry data that has been collected offline and is processed as though it were collected in real-time.

## Real-time Pupillometry

The `rtPupilPhase` directory provides the code and instructions needed to run the real-time pupillometry experimental code.

The version provided is for a fixation task, where a fixation cross is presented on a solid grey screen. In this task, participants are asked to maintain their gaze on the central fixation point at all times and keep a steady head position in the head/chin rest system. Participants completed five 10-minute fixation task blocks. Between each block participants were given break until they indicated they were prepared to begin the next fixation period.

## Simulated real-time pupillometry

The `simulations` directory provides the instructions, code and example for running the MATLAB script for *simulated* real-time pupillometry for human, monkey and mouse data. Although the parameters used in rtPupilPhase are robust for human pupillometry data, if you are recording data from other species, or you want to determine the effect of changing the parameters, you may wish to simulate the effects of the changes prior to collecting data.
