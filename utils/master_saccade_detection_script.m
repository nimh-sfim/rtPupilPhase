%% Find Saccades and Microssacades in EyeLink Data

% Written by: Sharif I. Kronemer with open access scripts
% Last Modified: 9/22/2023

%% Saccade detection

% Calculate velocity/acceleration using 5-sample window 
% (https://github.com/sj971/neurosci-saccade-detection/blob/master/analyzeEyeData.m)
% See Engbert and Kliegl, 2003. Denominator accounts for the
% six sample 'differences' used in numerator (i.e., n-2 to
% n+2 = 4 samples, n-1 to n+1 = 2 samples).

% Conversion to Degrees
% subtract center of screen (1024/2, 768/2)
% Y data is inverted
gaze_X = gaze_X_org-512;
gaze_Y = gaze_Y_org-384;

% Convert to mm, using pixel pitch (gotten from the monitor
% manufacturer's website) = .254mm per pixel
gaze_X = gaze_X*(.254);
gaze_Y = gaze_Y*(.254);

% Convert to deg (spherical coordinates)
gaze_X = atand(gaze_X/550);
gaze_Y = atand(gaze_Y/550);

% APPROACH 1

% Check vector length
if ~isequal(size(blink_data), size(gaze_Y), size(gaze_X))

    error('Blink and gaze data size mismatch!')

end

% Remove blink times from data
gaze_X_noblink = gaze_X;
gaze_X_noblink(blink_data) = nan;
gaze_Y_noblink = gaze_Y;
gaze_Y_noblink(blink_data) = nan;

% Interpolate over blinks
gaze_X_noblink = naninterp(gaze_X_noblink);
gaze_Y_noblink = naninterp(gaze_Y_noblink);

% EK Method (Input: X and Y gaze data [matrix: 2 x time], sampling rate, microsac threshold, microsac min duration)
% Note: allSac ouput contents
[allSac,msdx1,msdy1,stddevL,maddevL] = GetMicrosaccadesEK([gaze_X_noblink;gaze_Y_noblink],1000,5,5,blink_data);

% Create a saccade data vector
saccade_data = zeros(1,size(blink_data,2));

% Loop over all saccades
for n = 1:length(allSac)

    saccade_data(allSac(n,1):allSac(n,2)) = 1;

end

% Get microsaccades only
microsac_thres = 1;

% Initialize variables
microsac_idx = [];

% Loop over all saccades
for j = 1:size(allSac,1)

   % Find if saccade is within microsac threshold
   if sqrt(allSac(j,6)^2 + allSac(j,7)^2) < microsac_thres
   
       microsac_idx = [microsac_idx; j]; 

   end

end

% Store all microsaccade informaiton
allMicrosac = allSac(microsac_idx,:);

% Create a microsaccade data vector
microsaccade_data = zeros(1,size(blink_data,2));

% Loop over all microsaccades
for n = 1:length(allMicrosac)

    microsaccade_data(allMicrosac(n,1):allMicrosac(n,2)) = 1;

end

% Draw Figure

start = 75000;
stop = 90000;

figure;
hold on
ylim([0 1])
plot(saccade_data(start:stop),'b')
plot(microsaccade_data(start:stop),'c')
plot(blink_data(start:stop),'r')
plot(gaze_X(start:stop),'g')
plot(gaze_Y(start:stop),'y')