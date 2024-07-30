function [microsaccades, msdx,msdy,stddev, maddev] = GetMicrosaccadesEK(EyeDeg, SAMPLING, VFAC, MINDUR, blink_data)

%% INPUTS

%EyeDeg is 2xN long matrix of gaze location data (in degrees, Left and
%Right eye are input separately

%SAMPLING is the hertz of the data. Ex: 1000 for 1000Hz

%VFAC is the constant that will be multiplied against the E&K metric to
%form the microsaccade threshold (ellipse). Typically 3-5 (in the 2003
%paper E&K used 5), in the 2006 (E&M used 4). Lower = more sensitive

%MINDUR is the minimum duration of the microsaccade in INDICES, it's not
%milliseconds

%% OUTPUTS
% microsaccades: this the Nx7 structure described in the
% MicrosaccadeAnalysis2.m
% msdx = the E&K metrix for x dimension for this EYE / TRIAL
% msdy = the E&K metrix for x dimension for this EYE / TRIAL
% stddev = standard deviation of the EYE / TRIAL's velocity
% maddev = median absolute deviation of the EYE / TRIAL's velocity

% E&K algo detects ALL saccades, including large ones (> 1 deg)
% Microsaccade detection = Which data points exceed msdx * VFAC for more
% than MINDUR. (Of course msdy is also considered). 

% This code was adapted from code online. The current provenance has been lost; if this is your code, 
% please contact the authors and we will attribute it. 

% Reorient data and add a ones column
EyeDeg = EyeDeg';
EyeDeg = [ones(size(EyeDeg,1),1) EyeDeg]; %quirk of microsacc_plugin format

% Find the number of times
total_time = length(EyeDeg);

% Setup velocity matrix with 3 columns
velocity_matrix = zeros(total_time,3);

% Loop over time points and add 1s to velocity matrix
for time=1:total_time
    
    velocity_matrix(time,1)= EyeDeg(time,1);
    
end

%% Convert to velocity with a 6-sample smoothing

% Loop over the time points of eye gaze
for time=2:total_time-1
    
    % Don't consider times at the very beginning or end of data matrix
    % because smoothing want work
    if time>=3 & time<=total_time-2 %same as Engbert/Kliegel
    
        % Simultaneously applied to the X and Y gaze data
        velocity_matrix(time,2:3) = SAMPLING/6*[EyeDeg(time+2,2)+EyeDeg(time+1,2)-EyeDeg(time-1,2)-EyeDeg(time-2,2) EyeDeg(time+2,3)+EyeDeg(time+1,3)-EyeDeg(time-1,3)-EyeDeg(time-2,3)];
    
    end

end

% vel = sqrt(v(:,2).^2 + v(:,3).^2);

%% Call Engbert/Kliegel algorithm 

% Find saccades/microsaccades (Input: x/y gaze info, x/y velocity info,
% threshold, min duration threshold)
[microsaccades,~,msdx,msdy,stddev,maddev] = microsacc_plugin(EyeDeg(:,2:3),velocity_matrix(:,2:3),VFAC,MINDUR, blink_data);

end