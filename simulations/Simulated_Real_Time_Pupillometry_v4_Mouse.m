%% Simulated Real Time Pupil Phase - Mouse Data Set

% This script simulates running a rtPupilPhase on a previously acquired
% mouse pupillometry data set. The rtPupilPhase method parameters are 
% identical to those parameters used in humans.

% Written by: Sharif I. Kronemer
% Last Modified: 1/4/2024

clear all

tic

%% Directories and Paths

% Root path
root_path = pwd;

% Add paths
addpath(fullfile(root_path,'utils'))

% Data directory
data_dir = fullfile(root_path,'data','mouse');

% Output directory 
output_dir = fullfile(root_path,'analysis','subject_analysis','mouse');

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

% Select animals and sessions to run (columns = mice; rows = runs)  
test_sessions = [{'cm124_1_runB';'cm124_1_runC';'cm124_1_runH';'cm124_2_runB';'cm124_2_runD';'cm124_2_runE';'cm124_2_runF';'cm124_2_runG'},...
{'cm125_1_runB';'cm125_1_runC';'cm125_1_runD';'cm125_1_runE';'cm125_3_runB';'cm125_3_runC';'cm125_3_runD';'cm125_3_runE'},...
{'cm126_3_runB';'cm126_3_runC';'cm126_3_runH';'cm126_3_runI';'cm126_3_runJ';'cm126_6_runB';'cm126_6_runC';'cm126_6_runD'},...
{'cm127_1_runB';'cm127_1_runC';'cm127_1_runD';'cm127_1_runE';'cm127_1_runF';'cm127_1_runG';'cm127_1_runH';'cm127_1_runJ'},...
{'cm128_2_runB';'cm128_2_runC';'cm128_2_runD';'cm128_2_runE';'cm128_2_runF';'cm128_2_runG';'cm128_2_runH';'cm128_2_runI'}];

% Define pupillometry sampling rate
sampling_rate = 20; % in Hz
ms_per_sample = round(1000/sampling_rate); 

% Number of samples in block data
block_duration_samples = 11980;

% *** rtPupilPhase Parameters ***

% Pupil sample parameters
pupil_sample_duration_ms = 100; % in milliseconds
samples_in_pupil_sample = round(pupil_sample_duration_ms/ms_per_sample);

% Random event parameters
% Note: The number of random events specified and the block duration will
% determine the random IEI. Also, note that the number of random events
% selected will have implications on the number of pupil phase events.
num_random_events = 10;
random_IEI = block_duration_samples/num_random_events; % in milliseconds

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
IEI_jitter_samples = round(IEI_jitter_ms/ms_per_sample); % in samples

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
half_epoch_duration_ms = 2500;
half_epoch_duration_samples = round(half_epoch_duration_ms/ms_per_sample);

% Visualize trial data and pupil phase events (y = yes; n = no)
graph_trial_data = 'n';

%% Load Pupillometry Data

% Note: Mouse pupillometry data downloaded from https://www.sciencedirect.com/science/article/pii/S2211124723005387

% Find mat files - data filename
data_files = dir(fullfile(data_dir,'*.mat'));

% If one mat file found
if size(data_files, 1) == 1

    % Define filename
    filename = data_files.name;

% If multiple mat files found
elseif size(data_files, 1) > 1

    error('More than one pupil data file found!')

end

% Load data
load(fullfile(data_dir,filename))

% Count the number mouse-sessions
mouse_sessions = fieldnames(B2);

%% Begin rtPupilPhase Simulation

% Loop over mice
for mouse = 1:size(test_sessions,2)

    % Run analysis from scratch
    if isequal(approach,'R')

        disp(['Running mouse session ', num2str(mouse),'...'])
    
        % Find the index for the selected mouse-session above in the data list
        run_idx = find(ismember(mouse_sessions,test_sessions(:,mouse)));
    
        % Initialize variables
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
    
        % Loop over sessions
        for num = 1:size(test_sessions(:,mouse),1)
        
            % Define current mouse run index number
            current_run = run_idx(num);
        
            disp(['Running ', mouse_sessions{current_run},' ...'])
        
            % Extract pupil, whisking, and locomotion data
            run_pupil_data = eval(['B2.',mouse_sessions{current_run},'.pupil'])';
            run_whisk_data = eval(['B2.',mouse_sessions{current_run},'.whisk'])';
            run_loco_data = eval(['B2.',mouse_sessions{current_run},'.rotf'])';
    
            % Time data
            run_time_data = 1:length(run_pupil_data);

            % Find locomotion periods
            run_loco_idx = (getrunningpulses2(run_loco_data',0.1,20,100,"off"))';

            % Data size check
            if ~isequal(size(run_loco_idx,2),size(run_loco_data,2),size(run_pupil_data,2),size(run_time_data,2))

                error("Data vector sizes are not equal!")

            end

            %% Run Simulation
            
            % Initialize variables
        
            % Search window
            search_window_pupil = [];
            search_window_time = [];
            
            % Pupil data variables
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
            all_constriction_diff_fit = [];
            all_dilation_diff_fit = [];
            
            % Pupil event threshold
            peak_threshold_array = [];
            trough_threshold_array = [];
            dilation_threshold_array = [];
            constriction_threshold_array = [];
            
            % Loop over pupil sample windows
            for pupil_sample_num = 1:length(run_pupil_data)/samples_in_pupil_sample - 1
            
                %% Stage 1: Fill pupil sample
            
                % If the first pupil sample
                if pupil_sample_num == 1
            
                    % If first pupil samplewindow, start from the first index of data matrix
                    current_pupil_sample = run_pupil_data(1:pupil_sample_num*samples_in_pupil_sample);
                    current_pupil_sample_time = run_time_data(1:pupil_sample_num*samples_in_pupil_sample);

                    % Locomotion idx 
                    current_loco_sample = run_loco_idx(1:pupil_sample_num*samples_in_pupil_sample);
                    
                % If 2nd or later pupil sample
                else
            
                    % Count from +1 sample data from the previous pupil sample to last sample of the current pupil sample
                    current_pupil_sample = run_pupil_data(((pupil_sample_num-1)*samples_in_pupil_sample)+1:...
                        pupil_sample_num*samples_in_pupil_sample);
                    current_pupil_sample_time = run_time_data(((pupil_sample_num-1)*samples_in_pupil_sample)+1:...
                        pupil_sample_num*samples_in_pupil_sample);
                    
                    % Locomotion idx
                    current_loco_sample = run_loco_idx(((pupil_sample_num-1)*samples_in_pupil_sample)+1:pupil_sample_num*samples_in_pupil_sample);

                end
               
                % Replace 0s in pupil sample with NaN ("not a number")
                current_pupil_sample(current_pupil_sample <= 0) = nan;

                % Replace locomotion periods with NaN ("not a number")
                current_pupil_sample(current_loco_sample == 1) = nan;
                
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
            
                    % Check that not all values are NaN
                    if sum(~isnan(baseline_window_pupil_data)) > 2

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
                    
                    % If too many NaN values in the baseline window
                    else
            
                        % Reset baseline window
                        baseline_window_pupil_data = [];
            
                    end
            
                end
                
                % Reset the entire search window if NaN is found (i.e. blink or locomotion event)
                if any(isnan(search_window_pupil))
            
                    % Warning
                    if warning_flag == 1

                        warning('Blink/locomotion event detected - skipping search window!')
            
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

                        warning('Skipping search window becuase too long!')
            
                    end

                    % Reset all_buffers
                    search_window_pupil = [];
                    search_window_time = [];
                    diff_fit_vals = [];
                    search_sample_fit_vals = [];
            
                    % Skip to the next buffer
                    continue
                    
                end
            
                %% Stage 3 - Model search window pupil data with polynomial fit
            
                % Confirm there are at least 3 samples in search_window
                % Note: This is added because the sampling rate and pupil
                % sample size means that there can be a search_window
                % smaller than 3 samples that will not be accepted by fit
                if length(search_window_time) > 2
        
                    % Setup fitting sample vector
                    search_window_sample_vector = 1:length(search_window_pupil);
                
                    % Demean search window
                    demean_search_window_pupil = search_window_pupil - nanmean(search_window_pupil); 
                
                    % Fit data with a polynomial function
                    search_window_fit = fit(search_window_sample_vector',double(demean_search_window_pupil'),'poly2');
                
                    % Find the last pupil size value of the fitted curve
                    fit_value = search_window_fit(length(search_window_pupil));
            
                    % Store the fit value
                    search_sample_fit_vals = [search_sample_fit_vals; fit_value];
        
                else
        
                    continue
        
                end
                
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
        
                    % Store random info
                    all_random_idx = [all_random_idx; (pupil_sample_num*samples_in_pupil_sample)];
                    all_random_times = [all_random_times; search_window_time(end)];
                    all_random_pupil = [all_random_pupil; search_window_pupil(end)];
        
                    % Add to count
                    random_count = random_count + 1;
                
                end
        
                % Peak, Trough, Dilation, and Constriction Events
        
                % Find diff if there are more than two values (i.e., at least two search windows)
                if length(search_sample_fit_vals) > 1
            
                    % Fitting function last value diff
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
        
                        % Set found event idx
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
                        time_from_last_accepted_event = IEI_jitter_samples;
        
                    end
        
                    % Check if IEI time is exceeded
                    if time_from_last_accepted_event >= IEI_jitter_samples
        
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
        
            %% Visualize the Trial Pupil Data

            % Graph trial data
            if isequal(graph_trial_data, 'y')

                % Setup figure
                all_pupil_fig = figure
                hold on
                
                % Setup figure labels
                title(mouse_sessions{current_run})
                xlabel('Time (sampling rate 20Hz)')
                ylabel('Pupil size (uV)')
            
                % Plot data
                plot(run_pupil_data)

                % Plot all and accepted events
                scatter(all_peak_idx,run_pupil_data(all_peak_idx),'r')
                scatter(accepted_peak_idx,run_pupil_data(accepted_peak_idx),'r',"filled")
            
                scatter(all_trough_idx,run_pupil_data(all_trough_idx),'b')
                scatter(accepted_trough_idx,run_pupil_data(accepted_trough_idx),'b',"filled")
                
                scatter(all_dilation_idx,run_pupil_data(all_dilation_idx),'m')
                scatter(accepted_dilation_idx,run_pupil_data(accepted_dilation_idx),'m',"filled")
                
                scatter(all_constriction_idx,run_pupil_data(all_constriction_idx),'c')
                scatter(accepted_constriction_idx,run_pupil_data(accepted_constriction_idx),'cyan',"filled")
    
                % Save figure
                savefig(all_pupil_fig, fullfile(output_dir,['sim_rtPupilPhase_',mouse_sessions{current_run},'_all_block_pupil_timecourse.fig']))
    
                close
        
            end

            %% Cut Pupil Event Epochs
        
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
        
                % Begin cutting epoch
    
                % If accepted events are found
                if ~isempty(accepted_events_idx)
                
                    % Loop over events
                    for num_events = 1:length(accepted_events_idx)
                        
                        % Check epoch is inside data interval
                        if (accepted_events_idx(num_events)-half_epoch_duration_samples) < 1 || (accepted_events_idx(num_events)+half_epoch_duration_samples) > length(run_pupil_data)
                    
                            % Warning
                            if warning_flag == 1

                                warning('Event epoch oustide of data interval - skipping!')
                            
                            end

                            continue
                        
                        else
                
                            % Cut epochs centered on the pupil phase event time
                            pupil_epoch = run_pupil_data(accepted_events_idx(num_events)-half_epoch_duration_samples:accepted_events_idx(num_events)+half_epoch_duration_samples);
                            
                            % Demean pupil data
                            pupil_epoch = pupil_epoch - mean(pupil_epoch,"omitnan"); 
                
                        end
                
                        % Store epoch
                        eval(['accepted_',current_type,'_events_pupil_epochs(size(accepted_',current_type,'_events_pupil_epochs,1)+1,:) = pupil_epoch;'])
                
                    end
                
                end
        
            end
        
        end
    
        % Average across all epochs and add to mouse matrix
        accepted_mice_mean_peak_event_pupil_epochs(mouse,:) = nanmean(accepted_peak_events_pupil_epochs,1);
        accepted_mice_mean_trough_event_pupil_epochs(mouse,:) = nanmean(accepted_trough_events_pupil_epochs,1);
        accepted_mice_mean_constriction_event_pupil_epochs(mouse,:) = nanmean(accepted_constriction_events_pupil_epochs,1);
        accepted_mice_mean_dilation_event_pupil_epochs(mouse,:) = nanmean(accepted_dilation_events_pupil_epochs,1);
        accepted_mice_mean_random_event_pupil_epochs(mouse,:) = nanmean(accepted_random_events_pupil_epochs,1);
    
        % Save variables
        cd(output_dir)
        save(['sim_rtPupilPhase_mouse_',num2str(mouse),'_results.mat'], 'accepted*', 'all*',...
            'trough_threshold_array', 'peak_threshold_array', 'dilation_threshold_array', 'constriction_threshold_array')
    
        %% Save Stats
        
        % Navigate to output folder
        cd(output_dir)
        
        % Open text file
        stat_file = fopen(['sim_rtPupilPhase_Mouse_#',num2str(mouse),'_Stats.txt'],'w');
        fprintf(stat_file,'%s\n\r\n',['*** Mouse #',num2str(mouse),' ***']);
        
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

        %% Visualize Event Data
        
        % Plot all accepted event epoch timecourses

        % Timevector
        timevector = -half_epoch_duration_samples:half_epoch_duration_samples;
                
        % Define event types to plot
        event_types = {'peak','trough','dilation','constriction','random'};
        
        % Color values
        % Note: The order of the color list corresponds with the event list
        color_values = {'r','b','m','c','g'};
        
        % Loop over event types
        for type = 1:length(event_types)
        
            % Define current data
            pupil_epoch = eval(['accepted_',event_types{type},'_events_pupil_epochs']);
        
            % Setup figure
            epoch_fig = figure
            hold on
        
            % Setup figure labels
            title(['Mouse #',num2str(mouse),' ',event_types{type}],'Interpreter','none')
            ylabel('Demeaned Pupil Size (pixels)')
            xlabel('Time (ms)')
            
            % Axis limits
            xlim([-half_epoch_duration_samples, half_epoch_duration_samples])
            ylim([-3, 3])
            
            % Loop over epochs
            for epoch = 1:size(pupil_epoch,1)
                
               % Plot individual epoch
               plot(timevector, pupil_epoch(epoch,:),color_values{type})
            
            end
        
            % Plot reference lines
            stim_time = plot([0 0],[-half_epoch_duration_samples, half_epoch_duration_samples],'k')
            zero_line = plot([-half_epoch_duration_samples, half_epoch_duration_samples], [0, 0], 'k')
        
            % Save figure
            savefig(epoch_fig, fullfile(output_dir,['sim_rtPupilPhase_mouse_',num2str(mouse),'_pupil_',event_types{type},'_epoch_timecourse.fig']))
        
            close

        end

    % Load previously saved data
    elseif isequal(approach, 'L')

        disp(['Loading mouse #',num2str(mouse),' data...'])

        % Save variables
        cd(output_dir)
        load(['sim_rtPupilPhase_mouse_',num2str(mouse),'_results.mat'])

        % Average across all epochs
        accepted_mice_mean_peak_event_pupil_epochs(mouse,:) = nanmean(accepted_peak_events_pupil_epochs,1);
        accepted_mice_mean_trough_event_pupil_epochs(mouse,:) = nanmean(accepted_trough_events_pupil_epochs,1);
        accepted_mice_mean_constriction_event_pupil_epochs(mouse,:) = nanmean(accepted_constriction_events_pupil_epochs,1);
        accepted_mice_mean_dilation_event_pupil_epochs(mouse,:) = nanmean(accepted_dilation_events_pupil_epochs,1);
        accepted_mice_mean_random_event_pupil_epochs(mouse,:) = nanmean(accepted_random_events_pupil_epochs,1);

    else

        error('Approach method not found!')

    end

end

%% Group-Level Statistics

% Z-score pupil data
zscore_group_peak_pupil = zscore(accepted_mice_mean_peak_event_pupil_epochs,[],2);
zscore_group_trough_pupil = zscore(accepted_mice_mean_trough_event_pupil_epochs,[],2);
zscore_group_rising_pupil = zscore(accepted_mice_mean_dilation_event_pupil_epochs,[],2);
zscore_group_falling_pupil = zscore(accepted_mice_mean_constriction_event_pupil_epochs,[],2);
zscore_group_random_pupil = zscore(accepted_mice_mean_random_event_pupil_epochs,[],2);

% Average across mice
zscore_group_mean_peak_pupil = nanmean(zscore_group_peak_pupil,1);
zscore_group_mean_trough_pupil = nanmean(zscore_group_trough_pupil,1);
zscore_group_mean_rising_pupil = nanmean(zscore_group_rising_pupil,1);
zscore_group_mean_falling_pupil = nanmean(zscore_group_falling_pupil,1);
zscore_group_mean_random_pupil = nanmean(zscore_group_random_pupil,1);

group_mean_peak_pupil = nanmean(accepted_mice_mean_peak_event_pupil_epochs,1);
group_mean_trough_pupil = nanmean(accepted_mice_mean_trough_event_pupil_epochs,1);
group_mean_rising_pupil = nanmean(accepted_mice_mean_dilation_event_pupil_epochs,1);
group_mean_falling_pupil = nanmean(accepted_mice_mean_constriction_event_pupil_epochs,1);
group_mean_random_pupil = nanmean(accepted_mice_mean_random_event_pupil_epochs,1);

% Standard deviation across mice
zscore_group_SEM_peak_pupil = std(zscore_group_peak_pupil,0,1)/sqrt(size(zscore_group_peak_pupil,1));
zscore_group_SEM_trough_pupil = std(zscore_group_trough_pupil,0,1)/sqrt(size(zscore_group_trough_pupil,1));
zscore_group_SEM_rising_pupil = std(zscore_group_rising_pupil,0,1)/sqrt(size(zscore_group_rising_pupil,1));
zscore_group_SEM_falling_pupil = std(zscore_group_falling_pupil,0,1)/sqrt(size(zscore_group_falling_pupil,1));
zscore_group_SEM_random_pupil = std(zscore_group_random_pupil,0,1)/sqrt(size(zscore_group_random_pupil,1));

group_SEM_peak_pupil = std(accepted_mice_mean_peak_event_pupil_epochs,0,1)/sqrt(size(accepted_mice_mean_peak_event_pupil_epochs,1));
group_SEM_trough_pupil = std(accepted_mice_mean_trough_event_pupil_epochs,0,1)/sqrt(size(accepted_mice_mean_trough_event_pupil_epochs,1));
group_SEM_rising_pupil = std(accepted_mice_mean_dilation_event_pupil_epochs,0,1)/sqrt(size(accepted_mice_mean_dilation_event_pupil_epochs,1));
group_SEM_falling_pupil = std(accepted_mice_mean_constriction_event_pupil_epochs,0,1)/sqrt(size(accepted_mice_mean_constriction_event_pupil_epochs,1));
group_SEM_random_pupil = std(accepted_mice_mean_random_event_pupil_epochs,0,1)/sqrt(size(accepted_mice_mean_random_event_pupil_epochs,1));

% Save Data

% Save variables
cd(output_dir)
save sim_rtPupilPhase_results.mat zscore* group* accepted_mice*

%% Plot Z-score Pupil Mouse Timecourses

% Plot dimensions
xlim_ms = 1500;
xlim_samples = round(xlim_ms/ms_per_sample);

% Pupil events
pupil_events = {'rising','peak','falling','trough','random'};

% Color values
color_values = {'m','r','c','b','g'};

% Timevector
timevector = -half_epoch_duration_samples:half_epoch_duration_samples;

% Loop over subjects
for event = 1:length(pupil_events)

    % Define current data
    pupil_epoch = eval(['zscore_group_',pupil_events{event},'_pupil']);
    group_mean_pupil = eval(['zscore_group_mean_',pupil_events{event},'_pupil']);

    % Setup figure
    group_fig = figure
    hold on

    % Figure labels
    title('Simulated rtPupilPhase Timecourses - Z-score Pupil')
    ylabel('Z-score Pupil Size (uV)')
    xlabel('Time (ms)')
    
    % Axis limits
    xlim([-xlim_samples, xlim_samples])
    ylim([-3 3])
    
    % Plot reference lines
    plot([-xlim_samples, xlim_samples],[0 0],'k')
    plot([0 0],[-3 3],'k')
    
    % Loop over mice
    for mouse = 1:size(pupil_epoch,1)
        
        % Plot mouse timecourse
        plot(timevector, pupil_epoch(mouse,:),color_values{event})
    
    end
    
    % Mean timecourses
    plot(timevector, group_mean_pupil,color_values{event},'LineWidth',4)

    % Save figure
    savefig(group_fig,fullfile(output_dir,['sim_rtPupilPhase_',pupil_events{event},'_mouse_subject_zscore_pupil_timecourse.fig']))

end

%% Plot Z-score Pupil Group Timecourses

% Plot dimensions
xlim_ms = 1500;
xlim_samples = round(xlim_ms/ms_per_sample);

% Timevector
timevector = -half_epoch_duration_samples:half_epoch_duration_samples;

% Setup figure
pupil_fig = figure
hold on

% Setup figure labels
title('Simulated rtPupilPhase Timecourses - Z-score Pupil')
ylabel('Z-score Pupil Size (uV)')
xlabel('Time (ms)')
    
% Axis limits
xlim([-xlim_samples, xlim_samples])
ylim([-2.5, 2.5])

% Plot reference lines
stim_time = plot([0 0],[-xlim_samples, xlim_samples],'k')
zero_line = plot([-xlim_samples, xlim_samples], [0, 0], 'k')

% Mean timecourse
plot(timevector, zscore_group_mean_peak_pupil,'r','LineWidth',2)
plot(timevector, zscore_group_mean_trough_pupil,'b','LineWidth',2)
plot(timevector, zscore_group_mean_rising_pupil,'m','LineWidth',2)
plot(timevector, zscore_group_mean_falling_pupil,'c','LineWidth',2)
plot(timevector, zscore_group_mean_random_pupil,'g','LineWidth',2)

% Error timecourse
plot(timevector,zscore_group_mean_peak_pupil + zscore_group_SEM_peak_pupil,'r')
plot(timevector,zscore_group_mean_trough_pupil + zscore_group_SEM_trough_pupil,'b')
plot(timevector,zscore_group_mean_rising_pupil + zscore_group_SEM_rising_pupil,'m')
plot(timevector,zscore_group_mean_falling_pupil + zscore_group_SEM_falling_pupil,'c')
plot(timevector,zscore_group_mean_random_pupil + zscore_group_SEM_random_pupil,'g')

plot(timevector,zscore_group_mean_peak_pupil - zscore_group_SEM_peak_pupil,'r')
plot(timevector,zscore_group_mean_trough_pupil - zscore_group_SEM_trough_pupil,'b')
plot(timevector,zscore_group_mean_rising_pupil - zscore_group_SEM_rising_pupil,'m')
plot(timevector,zscore_group_mean_falling_pupil - zscore_group_SEM_falling_pupil,'c')
plot(timevector,zscore_group_mean_random_pupil - zscore_group_SEM_random_pupil,'g')

% Save figure
savefig(pupil_fig, fullfile(output_dir,'sim_rtPupilPhase_mouse_group_zscore_pupil_timecourse.fig'))

%% Plot Demeaned Pupil Group Timecourses

% Plot dimensions
xlim_ms = 1500;
xlim_samples = round(xlim_ms/ms_per_sample);

% Timevector
timevector = -half_epoch_duration_samples:half_epoch_duration_samples;

% Setup figure
pupil_fig = figure
hold on

% Setup figure labels
title('Simulated rtPupilPhase Timecourses - Demeaned Pupil')
ylabel('Demeaned Pupil Size (uV)')
xlabel('Time (ms)')
    
% Axis limits
xlim([-xlim_samples, xlim_samples])
ylim([-1, 1])

% Plot reference lines
stim_time = plot([0 0],[-xlim_samples, xlim_samples],'k')
zero_line = plot([-xlim_samples, xlim_samples], [0, 0], 'k')

% Mean timecourse
plot(timevector, group_mean_peak_pupil,'r','LineWidth',2)
plot(timevector, group_mean_trough_pupil,'b','LineWidth',2)
plot(timevector, group_mean_rising_pupil,'m','LineWidth',2)
plot(timevector, group_mean_falling_pupil,'c','LineWidth',2)
plot(timevector, group_mean_random_pupil,'g','LineWidth',2)

% Error timecourse
plot(timevector,group_mean_peak_pupil + group_SEM_peak_pupil,'r')
plot(timevector,group_mean_trough_pupil + group_SEM_trough_pupil,'b')
plot(timevector,group_mean_rising_pupil + group_SEM_rising_pupil,'m')
plot(timevector,group_mean_falling_pupil + group_SEM_falling_pupil,'c')
plot(timevector,group_mean_random_pupil + group_SEM_random_pupil,'g')

plot(timevector,group_mean_peak_pupil - group_SEM_peak_pupil,'r')
plot(timevector,group_mean_trough_pupil - group_SEM_trough_pupil,'b')
plot(timevector,group_mean_rising_pupil - group_SEM_rising_pupil,'m')
plot(timevector,group_mean_falling_pupil - group_SEM_falling_pupil,'c')
plot(timevector,group_mean_random_pupil - group_SEM_random_pupil,'g')

% Save figure
savefig(pupil_fig, fullfile(output_dir,'sim_rtPupilPhase_mouse_group_demeaned_pupil_timecourse.fig'))
