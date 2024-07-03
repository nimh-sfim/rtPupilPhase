%% Simulated Real Time Pupil Phase - Human Data Set

% This script simulates running rtPupilPhase on human pupillometry data
% during a live recording session. The parameters are set to copy those 
% used in rtPupilPhase. 
% 
% Note: Some parameters associated with this script are
% unique to the experimental task and pupillometry acquisition used in
% Kronemer et al., 2024. Therefore, updates to these parts of the script
% may be required when testing alternative data sets. 

% Written by: Sharif I. Kronemer
% Last Modified: 12/27/2023

clear all

tic

%% Directories and Paths

% Root path
root_path = pwd;

% Add paths (Note: Paths are added that houses functions used for pupil
% preprocessing and extracting saccade and microsaccade occurrences)
addpath(fullfile(root_path,'utils'))
%addpath(genpath(fullfile(root_path,'Real_Time_Perception_Study/Analysis/Analysis_Code/EyeLink_Analysis')))

%% Parameters

% Dictionary:
% ms = milliseconds
% IEI = inter-event interval
% num = number
% EyeLink = the pupillometry system (SR Research, Inc.)
% idx = index

% *** Subject and Recording Parameters ***

% Subject list
subject_list = {'046','048','073','074','078','079','080','081'};

% The recorded eye (1 = left; 2 = right)
% Note: EyeLink stores the pupil size data in a 2 x time/sample matrix. The
% first row = the left eye and the second row = the right eye.
recorded_eye = 2;

% Number of task blocks completed per participant
num_blocks = 5;
block_duration_ms = 600000; % in milliseconds

% Define pupillometry sampling rate
sampling_rate = 1000; % in Hz
ms_per_sample = 17; % EyeLink live recording rate is (60Hz or ~17 samples/s)
downsample_value = 17; % Downsample value

% *** Display Parameters ***

% Display monitor information
monitor_height = 1024; % in pixels
monitor_width = 768; % in pixels
pixel_pitch = 0.254; % number retrieved from the monitor manufacturer's website - .254mm per pixel

% *** rtPupilPhase Parameters ***

% Pupil sample parameters
pupil_sample_duration_ms = 100; % in milliseconds
samples_in_pupil_sample = round(pupil_sample_duration_ms/ms_per_sample);

% Random event parameters
% Note: The number of random events specified and the block duration will
% determine the random IEI. Also, note that the number of random events
% selected will have implications on the number of pupil phase events.
num_random_events = 15;
random_IEI = block_duration_ms/num_random_events; % in milliseconds

% Baseline window duration for setting new pupil size and derivative thresholds
baseline_window_ms = 5000; % in milliseconds
samples_in_baseline_window = round(baseline_window_ms/ms_per_sample);

% Pupil event thresholds 
% Note: These are default values that initiate the script simulation but
% are subsequently updated approximately each baseline interval.
peak_threshold = 0;
trough_threshold = 0;
dilation_threshold = 50;
constriction_threshold = -50;

% Quantile thresholds for pupil size and derivative
peak_pupil_quantile = 0.75;
trough_pupil_quantile = 0.25;
dilation_quantile = 0.99;
constriction_quantile = 0.01;

% Inter-event interval for pupil phase events
IEI_jitter_ms = 3000; % in milliseconds
samples_in_IEI = round(IEI_jitter_ms/ms_per_sample); % in samples

% Maximum length of search window
max_search_window_length_ms = 5000; % in milliseconds

% *** Other Parameters ***

% Warning flag (1 = yes; 0 = no)
% Note: There are multiple warnings scripted throughout the simulation that
% can help keep track of the simulation progress. However, for a cleaner
% terminal, you can supress these warnings. Scripted errors are not
% impacted by the warning flag.
warning_flag = 0;

% 50% of the epoch length to be extracted
half_epoch_duration_ms = 2500; % in milliseconds

% No blinks or saccades interval
no_blinks_saccades = 500; % in milliseconds

%% Begin rtPupilPhase Simulation

% Loop over subjects
for human = 1:length(subject_list)

    % Define the subject ID
    subID = subject_list{human};

    disp(['Running subject ', subID,'...'])

    % Data directory
    data_dir = fullfile(root_path,'data', 'human');
    
    % Output directory 
    output_dir = fullfile(root_path,'analysis','subject_analysis', 'human',subID);
    
    % Make output directory
    if ~exist(output_dir)
    
       mkdir(output_dir)
    
    end
    
    %% Load Pupil Data
    
    % Note: EyeLink data is stored as a EDF file and then manually
    % converted to mat formatting using edfmex (see details on the SR
    % Research support page: https://www.sr-research.com/support/thread-54.html
    
    disp('Loading pupil data...')

    % Find mat files - data filename
    edf_files = dir(fullfile(data_dir,'*.mat'));
    
    % If only one mat file
    if size(edf_files, 1) == 1
    
        filename = edf_files.name;
    
    % If more than one mat file
    elseif size(edf_files, 1) > 1
    
        error('More than one pupil data file found!')
    
    end

    % Load data
    load(fullfile(data_dir,filename))

    %% Extract Pupil Data

    % Note: EyeLink data is stored in a MATLAB struct. 

    disp('Extract pupil and gaze data...')
    
    % Time data
    time_data = edf_data.FSAMPLE.time;

    % Pupil data
    pupil_data = edf_data.FSAMPLE.pa(recorded_eye,:);

    % Gaze X and Y data
    gaze_Y_data = edf_data.FSAMPLE.gy(recorded_eye,:);
    gaze_X_data = edf_data.FSAMPLE.gx(recorded_eye,:);
    
    % Check the data is actually coming from the correct eye
    if length(unique(pupil_data)) < 10
    
        error('Low pupil data value diversity - might be using the wrong eye!')
    
    end
    
    % Note: Time vs log indices: time index is relative to the EyeLink data
    % timecourse; log index is relative to the EyeLink FEVENT event log.
    
    % Find task block onset times
    % Note: Specific to fixation task used in Kronemer et al., 2024.
    
    % Initialize block start and end variables
    start_block_times = [];
    end_block_times = [];
    
    % Loop over blocks
    for block = 1:num_blocks
    
        % Find block start string
        start_block_times = [start_block_times; edf_data.FEVENT(find(strcmp({edf_data.FEVENT(:).message}',...
            ['Starting Perception Rate Block ',num2str(block)]))).sttime];
        
        % Find block end string
        end_block_times = [end_block_times; edf_data.FEVENT(find(strcmp({edf_data.FEVENT(:).message}',...
            ['Finished Perception Rate Block ',num2str(block)]))).sttime];
    
    end
    
    % Initialize epoch variables
    % Note: "all" epochs store the timecourses for all detected events,
    % while "accepted" epochs store the timecourses for detected events
    % that exceed that IEI. The results reported in Kronemer et al., 2024
    % are from the accepted epochs. 

    % Pupil
    all_peak_event_pupil_epochs = [];
    all_trough_event_pupil_epochs = [];
    all_constriction_event_pupil_epochs = [];
    all_dilation_event_pupil_epochs = [];
    all_random_event_pupil_epochs = [];
    
    accepted_peak_event_pupil_epochs = [];
    accepted_trough_event_pupil_epochs = [];
    accepted_constriction_event_pupil_epochs = [];
    accepted_dilation_event_pupil_epochs = [];
    accepted_random_event_pupil_epochs = [];
    
    % Blinks
    all_peak_event_blink_epochs = [];
    all_trough_event_blink_epochs = [];
    all_constriction_event_blink_epochs = [];
    all_dilation_event_blink_epochs = [];
    all_random_event_blink_epochs = [];

    accepted_peak_event_blink_epochs = [];
    accepted_trough_event_blink_epochs = [];
    accepted_constriction_event_blink_epochs = [];
    accepted_dilation_event_blink_epochs = [];
    accepted_random_event_blink_epochs = [];

    % Saccades
    all_peak_event_saccade_epochs = [];
    all_trough_event_saccade_epochs = [];
    all_constriction_event_saccade_epochs = [];
    all_dilation_event_saccade_epochs = [];
    all_random_event_saccade_epochs = [];

    accepted_peak_event_saccade_epochs = [];
    accepted_trough_event_saccade_epochs = [];
    accepted_constriction_event_saccade_epochs = [];
    accepted_dilation_event_saccade_epochs = [];
    accepted_random_event_saccade_epochs = [];

    % Microsaccades
    all_peak_event_microsaccade_epochs = [];
    all_trough_event_microsaccade_epochs = [];
    all_constriction_event_microsaccade_epochs = [];
    all_dilation_event_microsaccade_epochs = [];
    all_random_event_microsaccade_epochs = [];

    accepted_peak_event_microsaccade_epochs = [];
    accepted_trough_event_microsaccade_epochs = [];
    accepted_constriction_event_microsaccade_epochs = [];
    accepted_dilation_event_microsaccade_epochs = [];
    accepted_random_event_microsaccade_epochs = [];
    
    % Setup variables to count the number of events
    peak_count = 0;
    trough_count = 0;
    dilation_count = 0;
    constriction_count = 0;
    random_count = 0;
    
    % Loop over blocks
    for block = 1:num_blocks
    
        disp(['Running block #', num2str(block),' ...'])
    
        % Start reviewing pupil data within recording blocks
        block_pupil_data = pupil_data(find(start_block_times(block) == time_data):find(end_block_times(block) == time_data));
        block_time_data = time_data(find(start_block_times(block) == time_data):find(end_block_times(block) == time_data));
        block_gaze_Y_data = gaze_Y_data(find(start_block_times(block) == time_data):find(end_block_times(block) == time_data));
        block_gaze_X_data = gaze_X_data(find(start_block_times(block) == time_data):find(end_block_times(block) == time_data)); 

        % Downsample data
        % Note: Downsampling matches the EyeLink real time sampling
        % rate in a live testing session (1000Hz to 60Hz)
        
        % Downsample pupil and time data
        downsampled_pupil_data = block_pupil_data(1:downsample_value:end);
        downsampled_time_data = block_time_data(1:downsample_value:end);
        
        %% Run Simulation
        
        % Initialize variables

        % Search window
        search_window_pupil = [];
        search_window_time = [];

        % Baseline window
        baseline_window_pupil_data = [];
    
        % Event sample number/index
        all_peak_idx = [];
        all_trough_idx = [];
        all_dilation_idx = [];
        all_constriction_idx = [];
        all_random_idx = [];
        accepted_peak_idx = [];
        accepted_trough_idx = [];
        accepted_dilation_idx = [];
        accepted_constriction_idx = [];
    
        % Event pupil size
        all_peak_pupil = [];
        all_trough_pupil = [];
        all_dilation_pupil = [];
        all_constriction_pupil = [];
        all_random_pupil = [];
    
        % Event times
        all_peak_times = [];
        all_trough_times = [];
        all_dilation_times = [];
        all_constriction_times = [];
        all_random_times = [];
        all_pupil_event_times = [];
            
        % Model fit parameters
        search_sample_fit_vals = [];
        all_trough_diff_fit = [];
        all_peak_diff_fit = [];
        all_dilation_diff_fit = [];
        all_constriction_diff_fit = [];
        
        % Pupil event threshold
        peak_threshold_array = [];
        trough_threshold_array = [];
        dilation_threshold_array = [];
        constriction_threshold_array = [];
        
        % Loop over pupil sample windows
        for pupil_sample_num = 1:length(downsampled_pupil_data)/samples_in_pupil_sample - 1
        
            %% Stage 1: Fill pupil sample
        
            % If the first pupil sample
            if pupil_sample_num == 1
        
                % If first pupil sample, start from the first index of the pupil data vector
                current_pupil_sample = downsampled_pupil_data(1:pupil_sample_num*samples_in_pupil_sample);
                current_pupil_sample_time = downsampled_time_data(1:pupil_sample_num*samples_in_pupil_sample);
                
            % If 2nd or later pupil sample
            else
        
                % Count from +1 sample data from the previous pupil sample to last sample of the current pupil sample
                current_pupil_sample = downsampled_pupil_data(((pupil_sample_num-1)*samples_in_pupil_sample)+1:...
                    pupil_sample_num*samples_in_pupil_sample);
                current_pupil_sample_time = downsampled_time_data(((pupil_sample_num-1)*samples_in_pupil_sample)+1:...
                    pupil_sample_num*samples_in_pupil_sample);
                
            end
           
            % Replace 0s in pupil sample with NaN ("not a number")
            current_pupil_sample(current_pupil_sample == 0) = nan;
            
            %% Stage 2: Create Search Window
        
            % Append pupil sample to search window
            search_window_pupil = [search_window_pupil, current_pupil_sample];
            search_window_time = [search_window_time, current_pupil_sample_time];
        
            % Append pupil sample to baseline window
            baseline_window_pupil_data = [baseline_window_pupil_data, current_pupil_sample];

            % Find the pupil phase event thresholds
            
            % If the number of samples in the baseline window exceeds the minimum
            % samples in the basline window, attempt to update the
            % thresholds
            if length(baseline_window_pupil_data) > samples_in_baseline_window
        
                % Check that not all values in baseline window are NaN
                if sum(isnan(baseline_window_pupil_data)) ~= length(baseline_window_pupil_data)
        
                    % Demean data
                    demean_baseline_data = baseline_window_pupil_data - mean(baseline_window_pupil_data,"omitnan");
            
                    % Get rid of NaNs
                    demean_baseline_data = demean_baseline_data(~isnan(demean_baseline_data));
            
                    % Find derivative of demeaned data
                    total_data_detrend_diff = diff(demean_baseline_data);
    
                    % Find peaks
                    [peaks,peaks_locs,w,peaks_prom] = findpeaks(demean_baseline_data);
                    
                    % Find troughs
                    % Note: The data is first inverted. 
                    [troughs,troughs_locs,w,troughs_prom] = findpeaks(-demean_baseline_data);
                    
                    % Flip trough pupil size files to restore original data
                    troughs = -troughs;
    
                    % Find new pupil size and derivative quantile thresholds
                    peak_threshold = quantile(peaks,peak_pupil_quantile);
                    trough_threshold = quantile(troughs,trough_pupil_quantile);
                    dilation_threshold = quantile(total_data_detrend_diff,dilation_quantile);
                    constriction_threshold = quantile(total_data_detrend_diff,constriction_quantile);
            
                    % Add updated threshold to array of threshold values
                    peak_threshold_array = [peak_threshold_array; peak_threshold];
                    trough_threshold_array = [trough_threshold_array; trough_threshold];
                    dilation_threshold_array = [dilation_threshold_array; dilation_threshold];
                    constriction_threshold_array = [constriction_threshold_array; constriction_threshold];
            
                    % Reset baseline window
                    baseline_window_pupil_data = [];
        
                % If all values in the baseline window is NaN
                else
        
                    % Reset baseline window
                    baseline_window_pupil_data = [];
        
                end
        
            end
            
            % Reset the search window if NaN is found (i.e. blink event)
            if any(isnan(search_window_pupil))
        
                % Warning
                if warning_flag == 1
                
                    warning('Blink/artifactual event detected - skipping search window!')
        
                end

                % Reset variables
                search_window_pupil = [];
                search_window_time = [];
                diff_fit_vals = [];
                search_sample_fit_vals = [];
        
                % Skip to the next search window
                continue
        
            % Reset if search window is too long    
            elseif length(search_window_pupil)*ms_per_sample > max_search_window_length_ms
        
                % Warning
                if warning_flag == 1

                    warning('Resetting search window because its too long!')
        
                end

                % Reset variables
                search_window_pupil = [];
                search_window_time = [];
                diff_fit_vals = [];
                search_sample_fit_vals = [];
        
                % Skip to the next search window
                continue
                
            end
        
            %% Stage 3 - Model search window pupil data with polynomial fit
        
            % Setup fitting sample vector
            search_window_sample_vector = 1:length(search_window_pupil);
        
            % Demean search window
            demean_search_window_pupil = search_window_pupil - mean(search_window_pupil,"omitnan"); 
        
            % Fit data with a polynomial function
            search_window_fit = fit(search_window_sample_vector',double(demean_search_window_pupil'),'poly2');
        
            % Find the last pupil size value of the fitted curve
            fit_value = search_window_fit(length(search_window_pupil));
    
            % Store the fit value
            search_sample_fit_vals = [search_sample_fit_vals; fit_value];
            
            %% Stage 4 - Find pupil events
    
            % Random Event
    
            % If there are previous random events
            if length(all_random_times) >= 1
            
                % Calculate time since the last random event
                time_from_last_random_event = double(search_window_time(end))-double(all_random_times(end));
    
            % The first random event
            else
    
                % Guarantees that the first stimulus triggers accepted random event
                time_from_last_random_event = random_IEI;
    
            end
    
            % Check if random IEI time is exceed
            if time_from_last_random_event >= random_IEI
    
                % Store random event info
                all_random_idx = [all_random_idx; (pupil_sample_num*samples_in_pupil_sample)];
                all_random_times = [all_random_times; search_window_time(end)];
                all_random_pupil = [all_random_pupil; search_window_pupil(end)];
    
                % Add to count
                random_count = random_count + 1;
            
            end
    
            % Peak, Trough, Dilation, and Constriction Events
    
            % Find derivative if there are more than two values (i.e., at least two search windows)
            if length(search_sample_fit_vals) > 1
        
                % Fitting function last value derivative
                diff_fit_vals = diff(search_sample_fit_vals);
        
                % Trough event: (1) Derivative between last and 2nd to last fit >
                % 0 (2) Pupil size of last demeaned sample in search window
                % is less than trough threshold 
                if diff_fit_vals(end) > 0 && demean_search_window_pupil(end) < trough_threshold
        
                    % Store info
                    all_trough_idx = [all_trough_idx; (pupil_sample_num*samples_in_pupil_sample)];
                    all_trough_times = [all_trough_times; search_window_time(end)];
                    all_pupil_event_times = [all_pupil_event_times; search_window_time(end)];
                    all_trough_pupil = [all_trough_pupil; search_window_pupil(end)];
                    all_trough_diff_fit = [all_trough_diff_fit; diff_fit_vals(end)];
    
                    % Count event number
                    trough_count = trough_count + 1;
    
                    % Set found event idx
                    found_event = 4;
                
                % Peak event: (1) Derivative between last and 2nd to last fit <
                % 0 (2) Pupil size of last demeaned sample in search window
                % is greater than peak threshold 
                elseif diff_fit_vals(end) < 0 && demean_search_window_pupil(end) > peak_threshold
        
                    % Store info
                    all_peak_idx = [all_peak_idx; (pupil_sample_num*samples_in_pupil_sample)];
                    all_peak_times = [all_peak_times; search_window_time(end)];
                    all_pupil_event_times = [all_pupil_event_times; search_window_time(end)];
                    all_peak_pupil = [all_peak_pupil; search_window_pupil(end)];
                    all_peak_diff_fit = [all_peak_diff_fit; diff_fit_vals(end)];
        
                    % Count event number
                    peak_count = peak_count + 1;
    
                    % Set found event idx
                    found_event = 2;
    
                % Dilation event: (1) Derivative between last and 2nd to last fit >
                % greater than dilation threshold
                elseif diff_fit_vals(end) > dilation_threshold
        
                    % Store info
                    all_dilation_idx = [all_dilation_idx; (pupil_sample_num*samples_in_pupil_sample)];
                    all_dilation_times = [all_dilation_times; search_window_time(end)];
                    all_dilation_pupil = [all_dilation_pupil; search_window_pupil(end)];
                    all_dilation_diff_fit = [all_dilation_diff_fit; diff_fit_vals(end)];
                    all_pupil_event_times = [all_pupil_event_times; search_window_time(end)];
        
                    % Count event number
                    dilation_count = dilation_count + 1;
    
                    % Set found event idx
                    found_event = 1;
    
                % Constriction event: (1) Derivative between last and 2nd to last fit >
                % greater than constriction threshold
                elseif diff_fit_vals(end) < constriction_threshold
        
                    % Store info
                    all_constriction_idx = [all_constriction_idx; (pupil_sample_num*samples_in_pupil_sample)];
                    all_constriction_times = [all_constriction_times; search_window_time(end)];
                    all_constriction_pupil = [all_constriction_pupil; search_window_pupil(end)];
                    all_constriction_diff_fit = [all_constriction_diff_fit; diff_fit_vals(end)];
                    all_pupil_event_times = [all_pupil_event_times; search_window_time(end)];
    
                    % Count event number
                    constriction_count = constriction_count + 1;
    
                    % Set found event
                    found_event = 3;
    
                % No event found
                else
    
                    % Set found event idx
                    found_event = 0;
    
                end
    
            % No event found
            else
    
                % Set found event idx
                found_event = 0;
        
            end
    
            % Check if an event was found of any type
            if found_event ~= 0
    
                % If a previous pupil event was found and an accepted pupil event
                if length(all_pupil_event_times) > 1 && exist('accepted_pupil_event_times','var')
                    
                    % Calculate time from the last accepted event
                    time_from_last_accepted_event = all_pupil_event_times(end)-accepted_pupil_event_times(end);
    
                % The first pupil phase event
                else
    
                    % Guarantee first stimulus triggers accepted event
                    time_from_last_accepted_event = IEI_jitter_ms;
    
                end
    
                % Check if IEI time is exceeded
                if time_from_last_accepted_event >= IEI_jitter_ms
    
                    % Log accepted event times
                    accepted_pupil_event_times = all_pupil_event_times(end);
    
                    % Store accepted event time index

                    % Dilation
                    if found_event == 1
    
                        accepted_dilation_idx = [accepted_dilation_idx; all_dilation_idx(end)];
    
                    % Peak
                    elseif found_event == 2
                        
                        accepted_peak_idx = [accepted_peak_idx; all_peak_idx(end)];
    
                    % Constriction
                    elseif found_event == 3
                        
                        accepted_constriction_idx = [accepted_constriction_idx; all_constriction_idx(end)];
    
                    % Trough
                    elseif found_event == 4
                        
                        accepted_trough_idx = [accepted_trough_idx; all_trough_idx(end)];
    
                    end
    
                    % Reset variables
                    search_window_pupil = [];
                    search_window_time = [];
                    diff_fit_vals = [];
                    search_sample_fit_vals = [];
    
                % IEI was not exceeded
                else

                    % Note: Will not reset variables, and will continue
                    % looking for an event. 

                    % Warning
                    if warning_flag == 1

                        warning('Event found but IEI time not yet exceeded!')

                    end
    
                end
    
            end
        
        end
    
        % Rename pupil data to affiliate with block
        eval(['Block_',num2str(block),'_pupil_data = all_block_pupil_data;'])
        
        %% Pupil Data Preprocessing
        
        % Pupil conversion value
        % Note: Stublinks expects pupil size values within the normal range 
        % of human pupil size in mm units. The pupil data is collected in
        % pixels with values in the 1000s. The conversion value temporary
        % updates the pupil size data for preprocessing and than restores
        % the original data after preprocessing.
        conversion_val = 1000;
        
        % Update pupil data 
        block_pupil_data = block_pupil_data/conversion_val; 
        
        % Process pupil data
        [processed_pupil_data, blink_data] = Stublinks60(block_pupil_data, sampling_rate);
        
        % Restore pupil data
        block_pupil_data = block_pupil_data*conversion_val;
        processed_pupil_data = processed_pupil_data*conversion_val;

        %% Saccade/Microsaccade Detection
    
        disp('Microsaccade and saccade detection...')
        
        % Calculate velocity/acceleration using 5-sample window 
        % (https://github.com/sj971/neurosci-saccade-detection/blob/master/analyzeEyeData.m)
        % See Engbert and Kliegl, 2003. Denominator accounts for the
        % six sample 'differences' used in numerator (i.e., n-2 to
        % n+2 = 4 samples, n-1 to n+1 = 2 samples).
        
        % Conversion to units of degree 
        
        % Subtract center of diplay screen from gaze data
        gaze_X = block_gaze_X_data-(monitor_width/2);
        gaze_Y = block_gaze_Y_data-(monitor_height/2);
        
        % Convert to mm, using pixel pitch (gotten from the monitor
        % manufacturer's website) = .254mm per pixel
        gaze_X = gaze_X*(pixel_pitch);
        gaze_Y = gaze_Y*(pixel_pitch);
        
        % Convert to degree (spherical coordinates)
        gaze_X = atand(gaze_X/550);
        gaze_Y = atand(gaze_Y/550);
        
        % Check data vector lengths
        if ~isequal(size(blink_data), size(gaze_Y), size(gaze_X))
        
            error('Blink and gaze data size mismatch!')
        
        end
        
        % Remove blink times from data - replace with NaN
        gaze_X_noblink = gaze_X;
        gaze_X_noblink(blink_data) = nan;
    
        gaze_Y_noblink = gaze_Y;
        gaze_Y_noblink(blink_data) = nan;
        
        % Interpolate over blinks
        gaze_X_noblink = naninterp(gaze_X_noblink);
        gaze_Y_noblink = naninterp(gaze_Y_noblink);
        
        % Run EK saccade/microsaccade extraction
        % Input: X and Y gaze data [matrix: 2 x time], sampling rate, microsac threshold, microsac min duration
        [allSac,msdx1,msdy1,stddevL,maddevL] = GetMicrosaccadesEK([gaze_X_noblink;gaze_Y_noblink],sampling_rate,5,5,blink_data);
        
        % Create a saccade data vector
        saccade_data = zeros(1,size(blink_data,2));
        
        % Loop over all saccades
        for n = 1:length(allSac)
        
            % Fill saccade index (1 = saccade event)
            saccade_data(allSac(n,1):allSac(n,2)) = 1;
        
        end
        
        % Get microsaccades only
        microsac_thres = 1;
        
        % Initialize variables
        microsac_idx = [];
        
        % Loop over all saccades
        for sac = 1:size(allSac,1)
        
           % Find if saccade is within microsaccade threshold
           if sqrt(allSac(sac,6)^2 + allSac(sac,7)^2) < microsac_thres
           
               % Add microsaccade event
               microsac_idx = [microsac_idx; sac]; 
        
           end
        
        end
        
        % Store all microsaccade informaiton
        allMicrosac = allSac(microsac_idx,:);
        
        % Create a microsaccade data vector
        microsaccade_data = zeros(1,size(blink_data,2));
        
        % Loop over all microsaccades
        for n = 1:length(allMicrosac)
        
            % Fill microsaccade index (1 = microsaccade event)
            microsaccade_data(allMicrosac(n,1):allMicrosac(n,2)) = 1;
        
        end
   
        %% Cut Pupil, Blink, Saccade, and Microsaccade Event Epochs 

        disp('Preparing to cut epochs ...')

        % Define event types
        event_types = {'random','dilation','peak','constriction','trough'};

        % Loop over event types
        for type = 1:length(event_types)
    
            % Set current type
            current_type = event_types{type};
    
            % Set generically named variable for event type index

            % If random event
            if isequal(current_type,'random')
    
                % Define current idx
                % Note: There are no "accepted" random events.
                all_events_idx = all_random_idx;
                accepted_events_idx = all_random_idx;

            % If non-random event
            else
    
                % Define current idx
                all_events_idx = eval(['all_',current_type,'_idx']);
                accepted_events_idx = eval(['accepted_',current_type,'_idx']);

            end
   
            % Convert index values to milliseconds
            all_events_idx = all_events_idx*ms_per_sample;
            accepted_events_idx = accepted_events_idx*ms_per_sample;

            % Begin cutting epoch

            % If accepted events are found
            if ~isempty(accepted_events_idx)
            
                % Loop over events
                for num_events = 1:length(accepted_events_idx)
                    
                    % Check epoch is inside data interval
                    if (accepted_events_idx(num_events)-half_epoch_duration_ms) < 1 ||...
                            (accepted_events_idx(num_events)+half_epoch_duration_ms) > length(processed_pupil_data)
            
                        % Warning
                        if warning_flag == 1

                            warning('Event epoch oustide of data interval - skipping this event!')

                        end

                        continue
                    
                    else
            
                        % Cut epochs centered on the pupil phase event time
                        pupil_epoch = processed_pupil_data(accepted_events_idx(num_events)-half_epoch_duration_ms:accepted_events_idx(num_events)+half_epoch_duration_ms);
                        blink_epoch = blink_data(accepted_events_idx(num_events)-half_epoch_duration_ms:accepted_events_idx(num_events)+half_epoch_duration_ms);
                        saccade_epoch = saccade_data(accepted_events_idx(num_events)-half_epoch_duration_ms:accepted_events_idx(num_events)+half_epoch_duration_ms);
                        microsaccade_epoch = microsaccade_data(accepted_events_idx(num_events)-half_epoch_duration_ms:accepted_events_idx(num_events)+half_epoch_duration_ms);
                        
                        % Demean pupil data
                        pupil_epoch = pupil_epoch - mean(pupil_epoch,"omitnan"); 

                    end

                    % Store epoch
                    eval(['accepted_',current_type,'_event_pupil_epochs(size(accepted_',current_type,'_event_pupil_epochs,1)+1,:) = pupil_epoch;'])
                    eval(['accepted_',current_type,'_event_blink_epochs(size(accepted_',current_type,'_event_blink_epochs,1)+1,:) = blink_epoch;'])
                    eval(['accepted_',current_type,'_event_saccade_epochs(size(accepted_',current_type,'_event_saccade_epochs,1)+1,:) = saccade_epoch;'])
                    eval(['accepted_',current_type,'_event_microsaccade_epochs(size(accepted_',current_type,'_event_microsaccade_epochs,1)+1,:) = microsaccade_epoch;'])
            
                end
            
            end
    
            % If any event was found (accepted or not)
            if ~isempty(all_events_idx)
            
                % Loop over events
                for num_events = 1:length(all_events_idx)
                    
                    % Check epoch is inside data interval
                    if (all_events_idx(num_events)-half_epoch_duration_ms) < 1 ||...
                            (all_events_idx(num_events)+half_epoch_duration_ms) > length(processed_pupil_data)
            
                        % Warning
                        if warning_flag == 1

                            warning('Event epoch oustide of data interval - skipping this event!')

                        end

                        continue
                    
                    else
                        
                        % Cut epochs centered on the pupil phase event time
                        pupil_epoch = processed_pupil_data(all_events_idx(num_events)-half_epoch_duration_ms:all_events_idx(num_events)+half_epoch_duration_ms);
                        blink_epoch = blink_data(all_events_idx(num_events)-half_epoch_duration_ms:all_events_idx(num_events)+half_epoch_duration_ms);
                        saccade_epoch = saccade_data(all_events_idx(num_events)-half_epoch_duration_ms:all_events_idx(num_events)+half_epoch_duration_ms);
                        microsaccade_epoch = microsaccade_data(all_events_idx(num_events)-half_epoch_duration_ms:all_events_idx(num_events)+half_epoch_duration_ms);
                        
                        % Demean pupil data
                        pupil_epoch = pupil_epoch - mean(pupil_epoch,"omitnan"); 

                    end

                    % Store epoch
                    eval(['all_',current_type,'_event_pupil_epochs(size(all_',current_type,'_event_pupil_epochs,1)+1,:) = pupil_epoch;'])
                    eval(['all_',current_type,'_event_blink_epochs(size(all_',current_type,'_event_blink_epochs,1)+1,:) = blink_epoch;'])
                    eval(['all_',current_type,'_event_saccade_epochs(size(all_',current_type,'_event_saccade_epochs,1)+1,:) = saccade_epoch;'])
                    eval(['all_',current_type,'_event_microsaccade_epochs(size(all_',current_type,'_event_microsaccade_epochs,1)+1,:) = microsaccade_epoch;'])
            
                end
            
            end

        end

        %% Visualize the Block Pupil Data and Events

        % Setup figure
        all_pupil_fig = figure
        hold on
        
        % Figure labels
        title([num2str(subID), ' - Block #',num2str(block)])
        xlabel('Time (ms)')
        ylabel('Pupil size (pixels)')
    
        % Plot block pupil data
        plot(block_pupil_data)

        % Plot all and accepted events
        scatter(all_peak_idx*ms_per_sample,block_pupil_data(all_peak_idx*ms_per_sample),'r')
        scatter(accepted_peak_idx*ms_per_sample,block_pupil_data(accepted_peak_idx*ms_per_sample),'r',"filled")
    
        scatter(all_trough_idx*ms_per_sample,block_pupil_data(all_trough_idx*ms_per_sample),'b')
        scatter(accepted_trough_idx*ms_per_sample,block_pupil_data(accepted_trough_idx*ms_per_sample),'b',"filled")

        scatter(all_dilation_idx*ms_per_sample,block_pupil_data(all_dilation_idx*ms_per_sample),'m')
        scatter(accepted_dilation_idx*ms_per_sample,block_pupil_data(accepted_dilation_idx*ms_per_sample),'m',"filled")
        
        scatter(all_constriction_idx*ms_per_sample,block_pupil_data(all_constriction_idx*ms_per_sample),'c')
        scatter(accepted_constriction_idx*ms_per_sample,block_pupil_data(accepted_constriction_idx*ms_per_sample),'c',"filled")

        % Save figure
        savefig(all_pupil_fig, fullfile(output_dir,['sim_rtPupilPhase_all_block_#',num2str(block),'_pupil_timecourse.fig']))

        close

    end
    
    %% Participant Analyses
    
    % Define event types
    event_types = {'random','dilation','peak','constriction','trough'};

    % Loop over event types
    for type = 1:length(event_types)

        % Set current type
        current_type = event_types{type};

        % Define current variables
        all_pupil_epochs =  eval(['all_',current_type,'_event_pupil_epochs']);
        all_blink_epochs = eval(['all_',current_type,'_event_blink_epochs']);
        all_saccade_epochs = eval(['all_',current_type,'_event_saccade_epochs']);
        all_microsaccade_epochs = eval(['all_',current_type,'_event_microsaccade_epochs']);

        accepted_pupil_epochs = eval(['accepted_',current_type,'_event_pupil_epochs']);
        accepted_blink_epochs = eval(['accepted_',current_type,'_event_blink_epochs']);
        accepted_saccade_epochs = eval(['accepted_',current_type,'_event_saccade_epochs']);
        accepted_microsaccade_epochs = eval(['accepted_',current_type,'_event_microsaccade_epochs']);

        % Find blinks and saccades in critical window
        critical_window = half_epoch_duration_ms-no_blinks_saccades:half_epoch_duration_ms;
        eye_movement_epochs_idx = find(sum(accepted_blink_epochs(:,critical_window),2) +...
            sum(accepted_saccade_epochs(:,critical_window),2)==0);

        % Find the mean across epochs within participant/across blocks
        mean_all_pupil_epochs = nanmean(all_pupil_epochs,1);
        mean_all_blink_epochs = nanmean(all_blink_epochs,1);
        mean_all_saccade_epochs = nanmean(all_saccade_epochs,1);
        mean_all_microsaccade_epochs = nanmean(all_microsaccade_epochs,1);

        mean_accepted_pupil_epochs = nanmean(accepted_pupil_epochs,1);
        mean_accepted_blink_epochs = nanmean(accepted_blink_epochs,1);
        mean_accepted_saccade_epochs = nanmean(accepted_saccade_epochs,1);
        mean_accepted_microsaccade_epochs = nanmean(accepted_microsaccade_epochs,1);
        mean_accepted_pupil_noblinks_nosacs_epochs = nanmean(accepted_pupil_epochs(eye_movement_epochs_idx,:),1);
        mean_accepted_blink_noblinks_nosacs_epochs = nanmean(accepted_blink_epochs(eye_movement_epochs_idx,:),1);
        mean_accepted_saccade_noblinks_nosacs_epochs = nanmean(accepted_saccade_epochs(eye_movement_epochs_idx,:),1);

        % Smooth saccade and microsaccade mean epochs (moving average)
        mean_all_saccade_epochs_smooth = smooth(mean_all_saccade_epochs,100)';
        mean_all_microsaccade_epochs_smooth = smooth(mean_all_microsaccade_epochs,100)';
    
        mean_accepted_saccade_epochs_smooth = smooth(mean_accepted_saccade_epochs,100)';
        mean_accepted_microsaccade_epochs_smooth = smooth(mean_accepted_microsaccade_epochs,100)';
        mean_accepted_saccade_noblinks_nosacs_epochs_smooth = smooth(mean_accepted_saccade_noblinks_nosacs_epochs,100)';

        % Rename generic variable with specific event type name
        eval(['blink_sac_',current_type,'_epochs_idx = eye_movement_epochs_idx;']);

        eval(['mean_all_',current_type,'_event_pupil_epochs = mean_all_pupil_epochs;']);
        eval(['mean_all_',current_type,'_event_blink_epochs = mean_all_blink_epochs;']);
        eval(['mean_all_',current_type,'_event_saccade_epochs = mean_all_saccade_epochs_smooth;']);
        eval(['mean_all_',current_type,'_event_microsaccade_epochs = mean_all_microsaccade_epochs_smooth;']);

        eval(['mean_accepted_',current_type,'_event_pupil_epochs = mean_accepted_pupil_epochs;']);
        eval(['mean_accepted_',current_type,'_event_blink_epochs = mean_accepted_blink_epochs;']);
        eval(['mean_accepted_',current_type,'_event_saccade_epochs = mean_accepted_saccade_epochs_smooth;']);
        eval(['mean_accepted_',current_type,'_event_microsaccade_epochs = mean_accepted_microsaccade_epochs_smooth;']);
        eval(['mean_accepted_',current_type,'_event_pupil_noblinks_nosacs_epochs = mean_accepted_pupil_noblinks_nosacs_epochs;']);
        eval(['mean_accepted_',current_type,'_event_blink_noblinks_nosacs_epochs = mean_accepted_blink_noblinks_nosacs_epochs;']);
        eval(['mean_accepted_',current_type,'_event_saccade_noblinks_nosacs_epochs = mean_accepted_saccade_noblinks_nosacs_epochs_smooth;']);

    end

    % Save variables
    cd(output_dir)
    save sim_rtPupilPhase_results.mat all*idx accepted*idx all_*epochs accepted_*epochs mean_* blink_sac* trough_threshold_array peak_threshold_array ...
        dilation_threshold_array constriction_threshold_array
    
    %% Save Stats
    
    % Navigate to output folder
    cd(output_dir)
    
    % Open text file
    stat_file = fopen('sim_rtPupilPhase_Fixation_Task_Stats.txt','w');
    fprintf(stat_file,'%s\n\r\n',['*** Subject ID - ',subID,' ***']);
    
    % *** Stats Text File ***
    fprintf(stat_file,'%s\n\r\n','All Dilations #');
    fprintf(stat_file,'%f\n\r',size(all_dilation_event_pupil_epochs,1));
    
    fprintf(stat_file,'%s\n\r\n','All Peaks #');
    fprintf(stat_file,'%f\n\r',size(all_peak_event_pupil_epochs,1));
    
    fprintf(stat_file,'%s\n\r\n','All Constrictions #');
    fprintf(stat_file,'%f\n\r',size(all_constriction_event_pupil_epochs,1));
    
    fprintf(stat_file,'%s\n\r\n','All Troughs #');
    fprintf(stat_file,'%f\n\r',size(all_trough_event_pupil_epochs,1));
    
    fprintf(stat_file,'%s\n\r\n','All Randoms #');
    fprintf(stat_file,'%f\n\r',size(all_random_event_pupil_epochs,1));

    fprintf(stat_file,'%s\n\r\n','Accepted Dilations #');
    fprintf(stat_file,'%f\n\r',size(accepted_dilation_event_pupil_epochs,1));
    
    fprintf(stat_file,'%s\n\r\n','Accepted Peaks #');
    fprintf(stat_file,'%f\n\r',size(accepted_peak_event_pupil_epochs,1));
    
    fprintf(stat_file,'%s\n\r\n','Accepted Constrictions #');
    fprintf(stat_file,'%f\n\r',size(accepted_constriction_event_pupil_epochs,1));
    
    fprintf(stat_file,'%s\n\r\n','Accepted Troughs #');
    fprintf(stat_file,'%f\n\r',size(accepted_trough_event_pupil_epochs,1));
    
    fprintf(stat_file,'%s\n\r\n','Accepted Randoms #');
    fprintf(stat_file,'%f\n\r',size(accepted_random_event_pupil_epochs,1));
    
    % Close text file
    fclose(stat_file);
    
    %% Visualize Event Data

    % Plot all accepted event epoch timecourses
    
    % Timevector
    timevector = -half_epoch_duration_ms:half_epoch_duration_ms;
        
    % Define event types to plot
    event_types = {'peak','trough','dilation','constriction','random'};
    
    % Color values
    % Note: The order of the color list corresponds with the event list
    color_values = {'r','b','m','c','g'};
    
    % Loop over event types
    for type = 1:length(event_types)
    
        % Define current data
        pupil_epoch = eval(['accepted_',event_types{type},'_event_pupil_epochs']);

        % Setup figure
        epoch_fig = figure
        hold on
    
        % Setup figure labels
        title([subID,' Accepted ',event_types{type},' Pupil Phase Event Epochs'])
        ylabel('Pupil Size (pixels)')
        xlabel('Time (ms)')
        
        % Axis limits
        xlim([-half_epoch_duration_ms, half_epoch_duration_ms])
        ylim([-1000, 1000])
        
        % Loop over epochs
        for epoch = 1:size(pupil_epoch,1)
            
            % Plot individual epoch
            plot(timevector, pupil_epoch(epoch,:),color_values{type})
        
        end
    
        % Plot reference lines
        stim_time = plot([0 0],[-half_epoch_duration_ms, half_epoch_duration_ms],'k')
        zero_line = plot([-half_epoch_duration_ms, half_epoch_duration_ms], [0, 0], 'k')
    
        % Save figure
        savefig(epoch_fig, fullfile(output_dir,['sim_rtPupilPhase_pupil_',event_types{type},'_epoch_timecourses.fig']))
    
        close

    end
    
    % Plot mean accepted event epoch timecourses

    % Pupil size figure
    
    % Setup figure
    pupil_fig = figure
    hold on
    
    % Setup figure labels
    title([subID,' Mean Accepted Pupil Phase Events'])
    ylabel('Pupil Size (pixels)')
    xlabel('Time (ms)')
    
    % Axis limits
    xlim([-half_epoch_duration_ms, half_epoch_duration_ms])
    ylim([-200, 200])
    
    % Plot reference line
    stim_time = plot([0 0],[-half_epoch_duration_ms, half_epoch_duration_ms],'k')
    zero_line = plot([-half_epoch_duration_ms, half_epoch_duration_ms], [0, 0], 'k')
    
    % Plot pupil timecourses
    plot(timevector, mean_accepted_peak_event_pupil_epochs,'r','LineWidth',2)
    plot(timevector, mean_accepted_trough_event_pupil_epochs,'b','LineWidth',2) 
    plot(timevector, mean_accepted_rising_event_pupil_epochs,'m','LineWidth',2)
    plot(timevector, mean_accepted_falling_event_pupil_epochs,'c','LineWidth',2)
    plot(timevector, mean_accepted_random_event_pupil_epochs,'g','LineWidth',2)
    
    % Save figure
    savefig(pupil_fig, fullfile(output_dir,'sim_rtPupilPhase_accepted_pupil_mean_timecourse.fig'))

    close

    % Blink occurrence figure
    
    % Setup figure
    blink_fig = figure
    hold on
    
    % Setup figure labels
    title([subID,' Mean Accepted Pupil Phase Events'])
    ylabel('Blink Fraction')
    xlabel('Time (ms)')
    
    % Axis limits
    xlim([-half_epoch_duration_ms, half_epoch_duration_ms])
    ylim([0, 1])
    
    % Plot reference line
    stim_time = plot([0 0],[0,1],'k')
    
    % Plot blink timecourses
    plot(timevector, mean_accepted_peak_event_blink_epochs,'r','LineWidth',2)
    plot(timevector, mean_accepted_trough_event_blink_epochs,'b','LineWidth',2)
    plot(timevector, mean_accepted_rising_event_blink_epochs,'m','LineWidth',2)
    plot(timevector, mean_accepted_falling_event_blink_epochs,'c','LineWidth',2)
    plot(timevector, mean_accepted_random_event_blink_epochs,'g','LineWidth',2)
    
    % Save figure
    savefig(blink_fig, fullfile(output_dir,'sim_rtPupilPhase_accepted_blink_mean_timecourse.fig'))

    close

end

toc