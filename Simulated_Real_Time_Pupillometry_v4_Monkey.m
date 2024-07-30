%% Simulated Real Time Pupillometry - Monkey Data Set

% This script simulates running a rtPupilPhase on a previously acquired
% monkey pupillometry data set. The rtPupilPhase method parameters have
% been adjusted to accomodate the monkey pupil size data eccentricities,
% however, the majority are identical to those parameters used in humans.

% Written by: Sharif I. Kronemer
% Last Modified: 1/5/2024

clear all

tic

%% Directories and Paths

% Root path
root_path = pwd; 

% Add paths (Note: Paths are added that houses functions used for pupil
% preprocessing)
addpath(fullfile(root_path,'utils'))

% Output directory 
output_dir = fullfile(root_path,'analysis', 'subject_analysis','monkey');

% Make output directory
if ~exist(output_dir)
   mkdir(output_dir)
end

%% Parameters

% Dictionary:
% ms = milliseconds
% IEI = inter-event interval
% num = number
% idx = index

% *** Subject and Recording Parameters ***

% Define the monkeys to test
subject_list = {'Monkey_1','Monkey_2'};

% Define pupillometry sampling rate
sampling_rate = 1000; % in Hz
ms_per_sample = round(1000/sampling_rate); 

% Minimal trial duration ms
min_trial_duration_ms = 5000;

% *** rtPupilPhase Parameters ***

% Pupil sample parameters
pupil_sample_duration_ms = 100; % in milliseconds
samples_in_pupil_sample = round(pupil_sample_duration_ms/ms_per_sample);

% Random event parameters
% Note: The number of random events specified and the trial duration will
% determine the random IEI. Also, note that the number of random events
% selected will have implications on the number of pupil phase events.
num_random_events = 1;
random_IEI = min_trial_duration_ms/num_random_events; % in milliseconds

% Baseline window duration for setting new pupil size and derivative thresholds
baseline_window_ms = 5000;%2500;%500; % in milliseconds
baseline_window_samples = round(baseline_window_ms/ms_per_sample); % in samples

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
IEI_jitter_ms = 3000;%1500; % in milliseconds
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

% Run analysis from scratch (R) or load previous data (L)
approach = 'R';

% 50% of the epoch length to be extracted
half_epoch_duration_ms = 500; % in milliseconds
half_epoch_duration_ms = round(half_epoch_duration_ms/ms_per_sample); % in samples

% Visualize trial data and pupil phase events (y = yes; n = no)
graph_trial_data = 'n';

%% Begin rtPupilPhase Simulation

% Loop over subjects
for monkey = 1:length(subject_list)

    % Run analyses from scratch
    if isequal(approach, 'R')

        disp(['Running monkey #',num2str(monkey),'...'])
    
        % Initialize variables
        % Note: "accepted" epochs store the timecourses for detected events
        % that exceed that IEI. The results reported in Kronemer et al., 2024
        % are from the accepted epochs. 
        accepted_peak_events_pupil_epochs = [];
        accepted_trough_events_pupil_epochs = [];
        accepted_constriction_events_pupil_epochs = [];
        accepted_dilation_events_pupil_epochs = [];
        accepted_random_events_pupil_epochs = [];
        
        % Setup variables to count the number of events
        peak_count = 0;
        trough_count = 0;
        dilation_count = 0;
        constriction_count = 0;
        random_count = 0;
        
        % Total tested time and test trial count
        tested_time = 0;
        tested_trial_count = 0;
    
        % Data directory
        data_dir = fullfile(root_path,'data', 'monkey',['Monkey_',num2str(monkey)]);
    
        %% Load and Process Pupil Data
        
        % Note: Monkey pupillometry data previously converted to MATLAB mat format
        
        % Find mat files - data filename
        data_files = dir(fullfile(data_dir,'*.mat'));
    
        % Loop over data files
        for files = 1%:size(data_files, 1)
    
            disp(['Running data file #',num2str(files),'...'])
    
            % Define filename
            filename = data_files(files).name;
    
            % Load data
            load(fullfile(data_dir,filename))
        
            % Adjust data vector dimensions
            if size(pupil,1) < size(pupil,2)
        
                % Reorient the pupil data
                pupil = pupil';
        
            end
    
            % Count the number trials
            num_trials = size(pupil,1);
        
            % Loop over trials
            for trial = 1:num_trials
            
                disp(['Running trial #', num2str(trial),' ...'])
            
                % Define pupil data and time data
                trial_pupil_data = pupil{trial,1}';
                trial_time_data = 1:length(trial_pupil_data);
    
                % Skip short trials
                if length(trial_pupil_data) < min_trial_duration_ms
            
                    continue
            
                end
    
                %% Stublink EyeLink Preprocessing
                
                % Pupil conversion value
                conversion_val = 3;
                
                % Update pupil data 
                trial_pupil_data = trial_pupil_data*conversion_val; 
                
                % Process data
                [processed_pupil_data, blink_data] = Stublinks60(trial_pupil_data, sampling_rate);
                %figure; hold on; plot(trial_pupil_data); plot(processed_pupil_data); plot(blink_data)
    
                % Restore the pupil data to pixel units 
                trial_pupil_data = trial_pupil_data/conversion_val;
                processed_pupil_data = processed_pupil_data/conversion_val;
    
                % Add the total time tested
                tested_time = tested_time + length(trial_pupil_data);
                tested_trial_count = tested_trial_count + 1;
    
                %% Run Simulation
                
                % Initialize variables
            
                % Search window
                search_window_pupil = [];
                search_window_time = [];
                
                % Baseline window
                baseline_window_pupil_data = [];
            
                % Pupil event sample #
                all_peak_idx = [];
                all_trough_idx = [];
                all_dilation_idx = [];
                all_constriction_idx = [];
                all_random_idx = [];
                accepted_peak_idx = [];
                accepted_trough_idx = [];
                accepted_dilation_idx = [];
                accepted_constriction_idx = [];
            
                % Pupil event pupil size
                all_peak_pupil = [];
                all_trough_pupil = [];
                all_dilation_pupil = [];
                all_constriction_pupil = [];
                all_random_pupil = [];
            
                % Pupil event times
                all_peak_times = [];
                all_trough_times = [];
                all_dilation_times = [];
                all_constriction_times = [];
                all_random_times = [];
                all_pupil_event_times = [];
                    
                % Model fit paramters
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
                for pupil_sample_num = 1:length(trial_pupil_data)/samples_in_pupil_sample - 1
                
                    %% Stage 1: Fill pupil sample
                
                    % If the first pupil sample
                    if pupil_sample_num == 1
                
                        % If first pupil sample, start from the first index of data matrix
                        current_pupil_sample = trial_pupil_data(1:pupil_sample_num*samples_in_pupil_sample);
                        current_pupil_sample_time = trial_time_data(1:pupil_sample_num*samples_in_pupil_sample);
    
                        % Blink sample idx 
                        current_blink_sample = blink_data(1:pupil_sample_num*samples_in_pupil_sample);
    
                    % If 2nd or later pupil sample
                    else
                
                        % Count from +1 sample data from the previous pupil sample to last sample of the current pupil sample
                        current_pupil_sample = trial_pupil_data(((pupil_sample_num-1)*samples_in_pupil_sample)+1:pupil_sample_num*samples_in_pupil_sample);
                        current_pupil_sample_time = trial_time_data(((pupil_sample_num-1)*samples_in_pupil_sample)+1:pupil_sample_num*samples_in_pupil_sample);
                        
                        % Blink sample idx
                        current_blink_sample = blink_data(((pupil_sample_num-1)*samples_in_pupil_sample)+1:pupil_sample_num*samples_in_pupil_sample);
    
                    end
    
                    % Replace blink periods with NaN ("not a number")
                    current_pupil_sample(current_blink_sample == 1) = nan;
            
                    %% Stage 2: Create Search Window
                
                    % Append pupil sample to search window
                    search_window_pupil = [search_window_pupil, current_pupil_sample];
                    search_window_time = [search_window_time, current_pupil_sample_time];
                
                    % Store all pupil data
                    %all_block_pupil_data = [all_block_pupil_data, current_pupil_sample];
            
                    % Store pupil data for baseline interval to reset threshold values
                    baseline_window_pupil_data = [baseline_window_pupil_data, current_pupil_sample];
                
                    % Find the pupil phase event thresholds
                    
                    % If the number of samples in the baseline window exceeds the minimum
                    % samples in the basline window, attempt to update the
                    % thresholds
                    if length(baseline_window_pupil_data) > baseline_window_samples
                
                        % Check that less than 50% of values are NaN
                        if sum(isnan(baseline_window_pupil_data)) < length(baseline_window_pupil_data)*0.5
                
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
                        
                        % If more than 50% of values in baseline window are NaN
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
                
                    % Setup fitting sample vector for all buffers
                    search_window_sample_vector = 1:length(search_window_pupil);
                
                    % Demean search window
                    demean_search_window_pupil = search_window_pupil - mean(search_window_pupil,"omitnan"); 
    
                    % Fit data with a polynomial function
                    search_window_fit = fit(search_window_sample_vector',double(demean_search_window_pupil'),'poly2');
                
                    % Find the last value of the fitting curve
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
                
                            % Store trough info
                            all_trough_idx = [all_trough_idx; (pupil_sample_num*samples_in_pupil_sample)];
                            all_trough_times = [all_trough_times; search_window_time(end)];
                            all_pupil_event_times = [all_pupil_event_times; search_window_time(end)];
                            all_trough_pupil = [all_trough_pupil; search_window_pupil(end)];
                            all_trough_diff_fit = [all_trough_diff_fit; diff_fit_vals(end)];
            
                            % Count event number
                            trough_count = trough_count + 1;
            
                            % Set found event
                            found_event = 4;
                                
                        % Peak event: (1) Derivative between last and 2nd to last fit <
                        % 0 (2) Pupil size of last demeaned sample in search window
                        % is greater than peak threshold 
                        elseif diff_fit_vals(end) < 0 && demean_search_window_pupil(end) > peak_threshold
                
                            % Store trough info
                            all_peak_idx = [all_peak_idx; (pupil_sample_num*samples_in_pupil_sample)];
                            all_peak_times = [all_peak_times; search_window_time(end)];
                            all_pupil_event_times = [all_pupil_event_times; search_window_time(end)];
                            all_peak_pupil = [all_peak_pupil; search_window_pupil(end)];
                            all_peak_diff_fit = [all_peak_diff_fit; diff_fit_vals(end)];
                
                            % Count event number
                            peak_count = peak_count + 1;
            
                            % Set found event
                            found_event = 2;
            
                        % Dilation event: (1) Derivative between last and 2nd to last fit >
                        % greater than dilation threshold
                        elseif diff_fit_vals(end) > dilation_threshold
                
                            % Store trough info
                            all_dilation_idx = [all_dilation_idx; (pupil_sample_num*samples_in_pupil_sample)];
                            all_dilation_times = [all_dilation_times; search_window_time(end)];
                            all_dilation_pupil = [all_dilation_pupil; search_window_pupil(end)];
                            all_dilation_diff_fit = [all_dilation_diff_fit; diff_fit_vals(end)];
                            all_pupil_event_times = [all_pupil_event_times; search_window_time(end)];
                
                            % Count event number
                            dilation_count = dilation_count + 1;
            
                            % Set found event
                            found_event = 1;
            
                        % Constriction event: (1) Derivative between last and 2nd to last fit >
                        % greater than constriction threshold
                        elseif diff_fit_vals(end) < constriction_threshold
                
                            % Store trough info
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

                        % Note: A unique adaptation made for the monkey
                        % data set is to consider the IEI within pupil
                        % phase events. This adjustment was made because
                        % the monkey pupil data segments are short and
                        % applying the IEI across pupil phase events would
                        % effectively allow ~1-2 pupil phase events per
                        % trial.
            
                        % If a previous pupil event was found
                        if length(all_pupil_event_times) > 1
                                
                            % Dilation
                            if found_event == 1
            
                                % Previous accepted events
                                if ~isempty(accepted_dilation_idx)
    
                                    % Calculate time from the last accepted event
                                    time_from_last_accepted_event = all_dilation_idx(end)-accepted_dilation_idx(end);
    
                                else
    
                                    % Guarantee accepted event
                                    time_from_last_accepted_event = IEI_jitter_ms;
    
                                end
                            
                            % Peak
                            elseif found_event == 2

                                % Previous accepted events
                                if ~isempty(accepted_peak_idx)
                                       
                                    % Calculate time from the last accepted event
                                    time_from_last_accepted_event = all_peak_idx(end)-accepted_peak_idx(end);
    
                                else
    
                                    % Guarantee accepted event
                                    time_from_last_accepted_event = IEI_jitter_ms;
    
                                end
    
                            % Constriction
                            elseif found_event == 3
                                
                                % Previous accepted events
                                if ~isempty(accepted_constriction_idx)
    
                                    % Calculate time from the last accepted event
                                    time_from_last_accepted_event = all_constriction_idx(end)-accepted_constriction_idx(end);
        
                                else 
    
                                    % Guarantee accepted event
                                    time_from_last_accepted_event = IEI_jitter_ms;
        
                                end
    
                            % Trough
                            elseif found_event == 4

                                % Previous accepted events
                                if ~isempty(accepted_trough_idx)
                                 
                                    % Calculate time from the last accepted event
                                    time_from_last_accepted_event = all_trough_idx(end)-accepted_trough_idx(end);
        
                                else 
                                    
                                    % Guarantee accepted event
                                    time_from_last_accepted_event = IEI_jitter_ms;
    
                                end
    
                            end
    
                        else
            
                            % Guarantee accepted event
                            time_from_last_accepted_event = IEI_jitter_ms;
            
                        end
            
                        % Check if IEI time is exceeded
                        if time_from_last_accepted_event >= IEI_jitter_ms
            
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
            
                %% Visualize the Trial Pupil Data

                % Graph trial data
                if isequal(graph_trial_data, 'y')

                    % Setup figure
                    figure
                    hold on
                    
                    % Figure labels
                    title(['Monkey #',num2str(monkey),' - Trial #',num2str(trial)])
                    xlabel('Time (ms)')
                    ylabel('Pupil size (arbitrary)')
                
                    % Plot pupil data
                    plot(trial_pupil_data)
    
                    % Plot all and accepted events
                    scatter(all_peak_idx,trial_pupil_data(all_peak_idx),'r')
                    scatter(accepted_peak_idx,trial_pupil_data(accepted_peak_idx),'r',"filled")
                
                    scatter(all_trough_idx,trial_pupil_data(all_trough_idx),'b')
                    scatter(accepted_trough_idx,trial_pupil_data(accepted_trough_idx),'b',"filled")
                    
                    scatter(all_dilation_idx,trial_pupil_data(all_dilation_idx),'m')
                    scatter(accepted_dilation_idx,trial_pupil_data(accepted_dilation_idx),'m',"filled")
                    
                    scatter(all_constriction_idx,trial_pupil_data(all_constriction_idx),'c')
                    scatter(accepted_constriction_idx,trial_pupil_data(accepted_constriction_idx),'cyan',"filled")
        
                    close 

                end
    
                %% Cut Pupil Size Event Epochs
            
                disp('Preparing to cut epochs ...')
           
                % Define event types
                event_types = {'random','dilation','peak','constriction','trough'};

                % Loop over event times
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
            
                    % Convert idx to ms (Note: At 1000Hz sampling rate
                    % there is no converstion)

                    % Begin cutting epoch
        
                    % If accepted events are found
                    if ~isempty(accepted_events_idx)
                        
                        % Loop over events
                        for num_events = 1:length(accepted_events_idx)
                            
                            % Check epoch is inside data interval
                            if (accepted_events_idx(num_events)-half_epoch_duration_ms) < 1 ||...
                                    (accepted_events_idx(num_events)+half_epoch_duration_ms) > length(trial_pupil_data)
    
                                % Warning
                                if warning_flag == 1
    
                                    warning('Event epoch oustide of data interval - skipping this event!')
    
                                end
    
                                continue
                            
                            else
                    
                                % Cut epochs centered on the pupil phase event time
                                pupil_epoch = processed_pupil_data(accepted_events_idx(num_events)-half_epoch_duration_ms:accepted_events_idx(num_events)+half_epoch_duration_ms);
    
                            end
                    
                            % Demean pupil data
                            pupil_epoch = pupil_epoch - mean(pupil_epoch,"omitnan"); 
    
                            % Detrend pupil data (Note: This is necessary in
                            % monkey data set because of a slow frequency
                            % drift from the beginning to end of an epoch)
                            pupil_epoch = detrend(pupil_epoch);
                            %figure; hold on; plot(pupil_epoch)

                            % Store epoch
                            eval(['accepted_',current_type,'_events_pupil_epochs(size(accepted_',current_type,'_events_pupil_epochs,1)+1,:) = pupil_epoch;'])

                        end
                    
                    end
            
                end
            
            end
    
        end
    
        % Average across all epochs
        accepted_monkey_mean_peak_event_pupil_epochs(monkey,:) = nanmean(accepted_peak_events_pupil_epochs,1);
        accepted_monkey_mean_trough_event_pupil_epochs(monkey,:) = nanmean(accepted_trough_events_pupil_epochs,1);
        accepted_monkey_mean_constriction_event_pupil_epochs(monkey,:) = nanmean(accepted_constriction_events_pupil_epochs,1);
        accepted_monkey_mean_dilation_event_pupil_epochs(monkey,:) = nanmean(accepted_dilation_events_pupil_epochs,1);
        accepted_monkey_mean_random_event_pupil_epochs(monkey,:) = nanmean(accepted_random_events_pupil_epochs,1);

        % Save variables
        cd(output_dir)
        save(['sim_rtPupilPhase_monkey_',num2str(monkey),'_results.mat'], 'all*', 'accepted*', 'tested*',...
            'trough_threshold_array', 'peak_threshold_array', 'dilation_threshold_array', 'constriction_threshold_array')
    
        %% Visualize data
        
        % Timevector
        timevector = -half_epoch_duration_ms:half_epoch_duration_ms;
        
        % All Epoch Timecourses
        
        % Event types to plot
        event_types = {'peak','trough','dilation','constriction','random'};
        
        % Color values
        color_values = {'r','b','m','c','g'};
        
        % Loop over events
        for type = 1:length(event_types)
        
            % Define current data
            pupil_epoch = eval(['accepted_',event_types{type},'_events_pupil_epochs']);
        
            % Figure
            epoch_fig = figure
            hold on
        
            % Labels
            title(['Monkey #',num2str(monkey),' ',event_types{type}])
            ylabel('Demeaned Pupil Size (arbitrary)')
            xlabel('Time (ms)')
            
            % Axis limits
            xlim([-half_epoch_duration_ms, half_epoch_duration_ms])
            ylim([-0.5, 0.5])
            
            % Loop over epochs
            for epoch = 1:size(pupil_epoch,1)
               
               % Plot individual epoch
               plot(timevector, pupil_epoch(epoch,:),color_values{type})
            
            end
        
            % Reference Line
            stim_time = plot([0 0],[-half_epoch_duration_ms, half_epoch_duration_ms],'k')
            zero_line = plot([-half_epoch_duration_ms, half_epoch_duration_ms], [0, 0], 'k')
        
            % Save figure
            savefig(epoch_fig, fullfile(output_dir,['sim_rtPupilPhase_monkey_',num2str(monkey),'_pupil_',event_types{type},'_epoch_timecourse.fig']))
        
            close
    
        end
    
        % Save Stats
        
        % Navigate to output folder
        cd(output_dir)
        
        % Open text file
        stat_file = fopen(['sim_rtPupilPhase_Monkey_#',num2str(monkey),'_Stats.txt'],'w');
        fprintf(stat_file,'%s\n\r\n',['*** Monkey #',num2str(monkey),' ***']);
        
        % *** Stats Text File ***
        fprintf(stat_file,'%s\n\r\n','All Dilations #');
        fprintf(stat_file,'%f\n\r',dilation_count);
        
        fprintf(stat_file,'%s\n\r\n','All Peaks #');
        fprintf(stat_file,'%f\n\r',peak_count);
        
        fprintf(stat_file,'%s\n\r\n','All Constrictions #');
        fprintf(stat_file,'%f\n\r',constriction_count);
        
        fprintf(stat_file,'%s\n\r\n','All Troughs #');
        fprintf(stat_file,'%f\n\r',trough_count);
        
        fprintf(stat_file,'%s\n\r\n','All Randoms #');
        fprintf(stat_file,'%f\n\r',random_count);
    
        fprintf(stat_file,'%s\n\r\n','Accepted Dilations #');
        fprintf(stat_file,'%f\n\r',size(accepted_dilation_events_pupil_epochs,1));
        
        fprintf(stat_file,'%s\n\r\n','Accepted Peaks #');
        fprintf(stat_file,'%f\n\r',size(accepted_peak_events_pupil_epochs,1));
        
        fprintf(stat_file,'%s\n\r\n','Accepted Constrictions #');
        fprintf(stat_file,'%f\n\r',size(accepted_constriction_events_pupil_epochs,1));
        
        fprintf(stat_file,'%s\n\r\n','Accepted Troughs #');
        fprintf(stat_file,'%f\n\r',size(accepted_trough_events_pupil_epochs,1));
        
        fprintf(stat_file,'%s\n\r\n','Accepted Randoms #');
        fprintf(stat_file,'%f\n\r',size(accepted_random_events_pupil_epochs,1));
        
        % Close text file
        fclose(stat_file);

    % Load previously saved data
    elseif isequal(approach, 'L')

        disp(['Loading monkey #',num2str(monkey),' data...'])

        % Save variables
        cd(output_dir)
        load(['sim_rtPupilPhase_monkey_',num2str(monkey),'_results.mat'])

        % Average across all epochs
        accepted_monkey_mean_peak_event_pupil_epochs(monkey,:) = nanmean(accepted_peak_events_pupil_epochs,1);
        accepted_monkey_mean_trough_event_pupil_epochs(monkey,:) = nanmean(accepted_trough_events_pupil_epochs,1);
        accepted_monkey_mean_constriction_event_pupil_epochs(monkey,:) = nanmean(accepted_constriction_events_pupil_epochs,1);
        accepted_monkey_mean_dilation_event_pupil_epochs(monkey,:) = nanmean(accepted_dilation_events_pupil_epochs,1);
        accepted_monkey_mean_random_event_pupil_epochs(monkey,:) = nanmean(accepted_random_events_pupil_epochs,1);

    else

        error('Approach method not found!')

    end

end

%% Group-Level Statistics

% Z-score pupil data
zscore_group_peak_pupil = zscore(accepted_monkey_mean_peak_event_pupil_epochs,[],2);
zscore_group_trough_pupil = zscore(accepted_monkey_mean_trough_event_pupil_epochs,[],2);
zscore_group_dilation_pupil = zscore(accepted_monkey_mean_dilation_event_pupil_epochs,[],2);
zscore_group_constriction_pupil = zscore(accepted_monkey_mean_constriction_event_pupil_epochs,[],2);
zscore_group_random_pupil = zscore(accepted_monkey_mean_random_event_pupil_epochs,[],2);

% Average across monkeys
zscore_group_mean_peak_pupil = nanmean(zscore_group_peak_pupil,1);
zscore_group_mean_trough_pupil = nanmean(zscore_group_trough_pupil,1);
zscore_group_mean_dilation_pupil = nanmean(zscore_group_dilation_pupil,1);
zscore_group_mean_constriction_pupil = nanmean(zscore_group_constriction_pupil,1);
zscore_group_mean_random_pupil = nanmean(zscore_group_random_pupil,1);

group_mean_peak_pupil = nanmean(accepted_monkey_mean_peak_event_pupil_epochs,1);
group_mean_trough_pupil = nanmean(accepted_monkey_mean_trough_event_pupil_epochs,1);
group_mean_dilation_pupil = nanmean(accepted_monkey_mean_dilation_event_pupil_epochs,1);
group_mean_constriction_pupil = nanmean(accepted_monkey_mean_constriction_event_pupil_epochs,1);
group_mean_random_pupil = nanmean(accepted_monkey_mean_random_event_pupil_epochs,1);

% Standard deviation across monkeys
zscore_group_SEM_peak_pupil = std(zscore_group_peak_pupil,0,1)/sqrt(size(zscore_group_peak_pupil,1));
zscore_group_SEM_trough_pupil = std(zscore_group_trough_pupil,0,1)/sqrt(size(zscore_group_trough_pupil,1));
zscore_group_SEM_dilation_pupil = std(zscore_group_dilation_pupil,0,1)/sqrt(size(zscore_group_dilation_pupil,1));
zscore_group_SEM_constriction_pupil = std(zscore_group_constriction_pupil,0,1)/sqrt(size(zscore_group_constriction_pupil,1));
zscore_group_SEM_random_pupil = std(zscore_group_random_pupil,0,1)/sqrt(size(zscore_group_random_pupil,1));

group_SEM_peak_pupil = std(accepted_monkey_mean_peak_event_pupil_epochs,0,1)/sqrt(size(accepted_monkey_mean_peak_event_pupil_epochs,1));
group_SEM_trough_pupil = std(accepted_monkey_mean_trough_event_pupil_epochs,0,1)/sqrt(size(accepted_monkey_mean_trough_event_pupil_epochs,1));
group_SEM_dilation_pupil = std(accepted_monkey_mean_dilation_event_pupil_epochs,0,1)/sqrt(size(accepted_monkey_mean_dilation_event_pupil_epochs,1));
group_SEM_constriction_pupil = std(accepted_monkey_mean_constriction_event_pupil_epochs,0,1)/sqrt(size(accepted_monkey_mean_constriction_event_pupil_epochs,1));
group_SEM_random_pupil = std(accepted_monkey_mean_random_event_pupil_epochs,0,1)/sqrt(size(accepted_monkey_mean_random_event_pupil_epochs,1));

% Save Data

% Save variables
cd(output_dir)
save sim_rtPupilPhase_results.mat zscore* group* 

%% Plot Z-score Pupil Group Timecourses

% Figure parameters
ymin = -2.5;
ymax = 2.5;
xmin = -1500;
xmax = 1500;

% Timevector
timevector = -half_epoch_duration_ms:half_epoch_duration_ms;

% Mean Pupil Phase Timecourses

% Pupil Figure
pupil_fig = figure
hold on

% Labels
title('Simulated rtPupilPhase Timecourses - Z-score Pupil')
ylabel('Z-scored Pupil Size')
xlabel('Time (ms)')

% Axis limits
xlim([-half_epoch_duration_ms, half_epoch_duration_ms])
ylim([ymin, ymax])
%yticks([-2.5,-2,-1,0,1,2,2.5])

% Reference Line
stim_time = plot([0 0],[-half_epoch_duration_ms, half_epoch_duration_ms],'k')
zero_line = plot([-half_epoch_duration_ms, half_epoch_duration_ms], [0, 0], 'k')

% Mean timecourse
plot(timevector, zscore_group_mean_peak_pupil,'r','LineWidth',2)
plot(timevector, zscore_group_mean_trough_pupil,'b','LineWidth',2)
plot(timevector, zscore_group_mean_dilation_pupil,'m','LineWidth',2)
plot(timevector, zscore_group_mean_constriction_pupil,'c','LineWidth',2)
plot(timevector, zscore_group_mean_random_pupil,'g','LineWidth',2)

% Error timecourse
plot(timevector,zscore_group_mean_peak_pupil + zscore_group_SEM_peak_pupil,'r')
plot(timevector,zscore_group_mean_trough_pupil + zscore_group_SEM_trough_pupil,'b')
plot(timevector,zscore_group_mean_dilation_pupil + zscore_group_SEM_dilation_pupil,'m')
plot(timevector,zscore_group_mean_constriction_pupil + zscore_group_SEM_constriction_pupil,'c')
plot(timevector,zscore_group_mean_random_pupil + zscore_group_SEM_random_pupil,'g')

plot(timevector,zscore_group_mean_peak_pupil - zscore_group_SEM_peak_pupil,'r')
plot(timevector,zscore_group_mean_trough_pupil - zscore_group_SEM_trough_pupil,'b')
plot(timevector,zscore_group_mean_dilation_pupil - zscore_group_SEM_dilation_pupil,'m')
plot(timevector,zscore_group_mean_constriction_pupil - zscore_group_SEM_constriction_pupil,'c')
plot(timevector,zscore_group_mean_random_pupil - zscore_group_SEM_random_pupil,'g')

% Save figure
savefig(pupil_fig, fullfile(output_dir,'sim_rtPupilPhase_monkey_group_zscore_pupil_timecourse.fig'))

%% Plot Z-score Pupil Monkey Timecourses

% Plot dimensions
xlim_ms = 500;
xlim_samples = round(xlim_ms/ms_per_sample);
ymin = -3;
ymax = 3;

% Pupil events
pupil_events = {'dilation','peak','constriction','trough','random'};

% Color values
color_values = {'m','r','c','b','g'};

% Timevector
timevector = -half_epoch_duration_ms:half_epoch_duration_ms;

% Loop over subjects
for event = 1:length(pupil_events)

    % Define current data
    pupil_epoch = eval(['zscore_group_',pupil_events{event},'_pupil']);
    group_mean_pupil = eval(['zscore_group_mean_',pupil_events{event},'_pupil']);

    % Setup figure
    group_fig = figure
    hold on
    
    % Figure labels
    title(['Simulated rtPupilPhase Timecourses - ',pupil_events{event},' - Z-score Pupil'])
    ylabel('Z-score Pupil Size')
    xlabel('Time (ms)')
    
    % Axis limits
    xlim([-xlim_samples, xlim_samples])
    ylim([ymin ymax])
    
    % Plot reference lines
    plot([-xlim_samples, xlim_samples],[0 0],'k')
    plot([0 0],[ymin ymax],'k')

    % Loop over mice
    for monkey = 1:size(pupil_epoch,1)
        
        % Plot mouse timecourse
        plot(timevector, pupil_epoch(monkey,:),color_values{event})
    
    end

    % Mean timecourses
    plot(timevector, group_mean_pupil,color_values{event},'LineWidth',4)

    % Save figure
    savefig(group_fig,fullfile(output_dir,['sim_rtPupilPhase_',pupil_events{event},'_monkey_subject_zscore_pupil_timecourse.fig']))

end


%% Plot Demeaned Pupil Group Timecourses

% Figure parameters
ymin = -0.035;
ymax = 0.035;
xmin = -1500;
xmax = 1500;

% Timevector
timevector = -half_epoch_duration_ms:half_epoch_duration_ms;

% Mean Pupil Phase Timecourses

% Setup figure
pupil_fig = figure
hold on

% Setup figure labels
title('Simulated rtPupilPhase Timecourses - Demeaned Pupil')
ylabel('Demeaned Pupil Size (arbitrary)')
xlabel('Time (ms)')

% Axis limits
xlim([-half_epoch_duration_ms, half_epoch_duration_ms])
ylim([ymin, ymax])

% Plot reference lines
stim_time = plot([0 0],[-half_epoch_duration_ms, half_epoch_duration_ms],'k')
zero_line = plot([-half_epoch_duration_ms, half_epoch_duration_ms], [0, 0], 'k')

% Mean timecourse
plot(timevector, group_mean_peak_pupil,'r','LineWidth',2)
plot(timevector, group_mean_trough_pupil,'b','LineWidth',2)
plot(timevector, group_mean_dilation_pupil,'m','LineWidth',2)
plot(timevector, group_mean_constriction_pupil,'c','LineWidth',2)
plot(timevector, group_mean_random_pupil,'g','LineWidth',2)

% Error timecourse
plot(timevector,group_mean_peak_pupil + group_SEM_peak_pupil,'r')
plot(timevector,group_mean_trough_pupil + group_SEM_trough_pupil,'b')
plot(timevector,group_mean_dilation_pupil + group_SEM_dilation_pupil,'m')
plot(timevector,group_mean_constriction_pupil + group_SEM_constriction_pupil,'c')
plot(timevector,group_mean_random_pupil + group_SEM_random_pupil,'g')

plot(timevector,group_mean_peak_pupil - group_SEM_peak_pupil,'r')
plot(timevector,group_mean_trough_pupil - group_SEM_trough_pupil,'b')
plot(timevector,group_mean_dilation_pupil - group_SEM_dilation_pupil,'m')
plot(timevector,group_mean_constriction_pupil - group_SEM_constriction_pupil,'c')
plot(timevector,group_mean_random_pupil - group_SEM_random_pupil,'g')

% Save figure
savefig(pupil_fig, fullfile(output_dir,'sim_rtPupilPhase_monkey_group_demeaned_pupil_timecourse.fig'))

toc