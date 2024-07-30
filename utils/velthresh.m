function [msdx msdy stddev maddev ] = velthresh(vel)
% VELTHRESH - compute median-based standard deviation (SD) estimators
% of velocity time series
%
%-------------------------------------------------------------------
% INPUT:
%
% vel: n-by-2 matrix with hor. [vel(:,1)] and vert. [vel(:,2)] eye velocity
%
% OUTPUT:
%
% msdx: median-based SD estimator for horizontal eye velocity
% msdy: median-based SD estimator for vertical eye velocity
%-------------------------------------------------------------------
%
% Note by olaf.dimigen@hu-berlin.de (OD): The code below was originally 
% part of the microsacc() function by Engbert et al. but was outsourced 
% into a stand-alone function for usage with the plugin. This allows to
% compute velocity thresholds globally (based on all epochs of a given 
% participant) or to work with fixed rather than relative thresholds.

% Code taken from : https://github.com/olafdimigen/eye-eeg

% compute threshold
msdx = sqrt( nanmedian(vel(:,1).^2) - (nanmedian(vel(:,1)))^2 );
msdy = sqrt( nanmedian(vel(:,2).^2) - (nanmedian(vel(:,2)))^2 );
if msdx<realmin
    msdx = sqrt( nanmean(vel(:,1).^2) - (nanmean(vel(:,1)))^2 );
    if msdx<realmin
        error('msdx<realmin in microsacc.m');
        error('msdx<realmin in vecthresh.m. Did you exclude blinks/missing data before saccade detection?'); % // added by O.D.
    end
end
if msdy<realmin
    msdy = sqrt( nanmean(vel(:,2).^2) - (nanmean(vel(:,2)))^2 );
    if msdy<realmin
        %error('msdy<realmin in microsacc.m');
        error('msdy<realmin in vecthresh.m. Did you exclude blinks/missing data before saccade detection?'); % // added by O.D.
    end
end

%% Added by Jackson
stddev = [nanstd(vel(:,1)) nanstd(vel(:,2))];
maddev = [mad(vel(:,1),1) mad(vel(:,2),1)];

end