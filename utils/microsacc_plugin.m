function [sac, radius, msdx, msdy, stddev, maddev] = microsacc_plugin(x,vel,VFAC,MINDUR,blink_data)
%(x,vel,VFAC,MINDUR,MSDX,MSDY,blink_data)
%-------------------------------------------------------------------
%
%  FUNCTION microsacc.m
%  Detection of monocular candidates for microsaccades;
%  Please cite: Engbert, R., & Mergenthaler, K. (2006) Microsaccades 
%  are triggered by low retinal image slip. Proceedings of the National 
%  Academy of Sciences of the United States of America, 103: 7192-7197.
%
%  (Version 2.1, 03 OCT 05)
%  
%-------------------------------------------------------------------
%
%  INPUT:
%
%  x(:,1:2)         position vector
%  vel(:,1:2)       velocity vector
%  VFAC             relative velocity threshold
%  MINDUR           minimal saccade duration

%  Optional inputs (modification by OD for plugin:)
%  MSDX             threshold for x-component (horizontal)
%  MSDY             threshold for y-component (vertical)

%  OUTPUT:
%
%  sac(1:num,1)   onset of saccade
%  sac(1:num,2)   end of saccade
%  sac(1:num,3)   peak velocity of saccade (vpeak)
%  sac(1:num,4)   horizontal component     (dx)
%  sac(1:num,5)   vertical component       (dy)
%  sac(1:num,6)   horizontal amplitude     (dX)
%  sac(1:num,7)   vertical amplitude       (dY)
%
%
%-------------------------------------------------------------------
% Note by olaf.dimigen@hu-berlin.de (2012):
% This version of microsacc() has been modified for use with the plugin. 
% msdx and msdy can now be provided as function inputs. Their computation 
% is outsourced into new function velthresh(), which contains code that was 
% formerly part of microsacc(). Helpful to compute threholds globally 
% (over all data epochs of a subject) or to apply fixed thresholds.
%-------------------------------------------------------------------

% // changes to original microsacc()
% // olaf.dimigen@hu-berlin.de
% https://github.com/olafdimigen/eye-eeg

% Compute the thresholds
%if nargin < 5

[msdx msdy stddev maddev] = velthresh(vel); %original:   [msdx msdy] = velthreshold(vel);

%else
%    msdx = MSDX;
%    msdy = MSDY;

%end

% // end of changes

radiusx = VFAC*msdx;
radiusy = VFAC*msdy;
radius = [radiusx radiusy];

% Compute test criterion: ellipse equation
test = (vel(:,1)/radiusx).^2 + (vel(:,2)/radiusy).^2;

% Remove all blink periods from consideration as a saccade
if isequal(size(blink_data,2),size(test,1))

    test(blink_data) = 0;

else

    error('Cannot apply blink exclusion!')

end

% Find index of test greater than 1
indx = find(test>1);

% Determine saccades
% Note: the code will look for continuous time points from the indx
% variable and identify the start and end of them and confirm that breach
% the min number of time points.

% Initialize varialbles
N = length(indx); 
sac = [];
nsac = 0;
dur = 1;
a = 1;
k = 1;

% Check all instances of N
while k<N
    
    % If the distance between indices is 1
    if indx(k+1)-indx(k)==1
        
        dur = dur + 1; 
    
    else

        % Check if the dur breaches threshold
        if dur>=MINDUR
            
            % Count # saccades
            nsac = nsac + 1;
            
            % Add saccade instance start/stop
            b = k;
            sac(nsac,:) = [indx(a) indx(b)];
        
        end
        
        % Update/reset variables
        a = k+1;
        dur = 1;

    end
    
    % Update k value
    k = k + 1;

end

% Check for minimum duration for last indx event
if dur>=MINDUR
    % Count # saccades
    nsac = nsac + 1;

    % Add saccade instance start/stop
    b = k;
    sac(nsac,:) = [indx(a) indx(b)];

end

% Compute peak velocity, horizonal and vertical components
% Note: sac variable columns are as follows: Onset time, offset time, peak
% velocity, x vector, 
% Loop over saccade events
for s=1:nsac

    % Onset and offset
    a = sac(s,1); 
    b = sac(s,2); 

    % Saccade peak velocity (vpeak)
    vpeak = max( sqrt( vel(a:b,1).^2 + vel(a:b,2).^2 ) );

    % Store result
    sac(s,3) = vpeak;

    % Saccade vector direction (dx,dy)

    % Find start and end X gaze degree
    dx = x(b,1)-x(a,1); 

    % Find start and end Y gaze degree
    dy = x(b,2)-x(a,2); 

    % Store result
    sac(s,4) = dx;
    sac(s,5) = dy;

    % Saccade amplitude (dX,dY)
    
    % Find the start to stop saccade time range
    i = sac(s,1):sac(s,2);

    % Find the min/max position in saccade period X and Y gaze degree
    [minx, ix1] = min(x(i,1));
    [maxx, ix2] = max(x(i,1));
    [miny, iy1] = min(x(i,2));
    [maxy, iy2] = max(x(i,2));

    % Calculate amplitude
    dX = sign(ix2-ix1)*(maxx-minx);
    dY = sign(iy2-iy1)*(maxy-miny);

    % Store result
    sac(s,6:7) = [dX dY];

end
