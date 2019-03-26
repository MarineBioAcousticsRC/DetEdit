function [nb,eb,sb,bd] = calculate_bouts(clickTimes,p)

% calculate_bouts.m

% Takes detection times and parameters and outputs start, end, duration and
% number of bout

% Inputs:
%
%   clickTimes - An Nx1 vector of detection times, where N is the number of
%           detection times.
%
%   p - A struct with parameter settings.
%
%
%
% Outputs:
%
%   nb - Number of bouts.
%
%   eb - An Nx1 vector of end times of bouts, where N is the number of
%       bouts.
%
%   sb - An Nx1 vector of start times of bouts.
%
%   bd - An Nx1 vector of duration of bouts in days.



% find edges (start and end times) of bouts or sessions
dt = diff(clickTimes)*24*60*60; % calculate time between detections in seconds
gt = p.gth*60*60;    % gap time in sec
I = find(dt>gt);  % find start of gaps
sb = [clickTimes(1);clickTimes(I+1)];   % start time of bout
eb = [clickTimes(I);clickTimes(end)];   % end time of bout
dd = clickTimes(end)-clickTimes(1);     % deployment duration [d]
nb = length(sb);        % number of bouts
bd = (eb - sb);      % duration of bout in days

% find bouts longer than the minimum
if ~isempty(p.minBout)
    bdI = find(bd > (p.minBout / (60*60*24)));
    bd = bd(bdI);
    sb = sb(bdI);
    eb = eb(bdI);
    nb = length(sb);
end

% limit the length of a bout
blim = p.ltsaMax/24;       % 6 hr bout length limit in days
ib = 1;
while ib <= nb
    bd = (eb - sb);   %duration bout in days
    if (bd(ib) > blim)      % find long bouts
        nadd = ceil(bd(ib)/blim) - 1; % number of bouts to add
        for imove = nb : -1: (ib +1)
            sb(imove+nadd)= sb(imove);
        end
        for iadd = 1 : 1: nadd
            sb(ib+iadd) = sb(ib) + blim*iadd;
        end
        for imove = nb : -1 : ib
            eb(imove+nadd) = eb(imove);
        end
        for iadd = 0 : 1 : (nadd - 1)
            eb(ib+iadd) = sb(ib) + blim*(iadd+1);
        end
        nb = nb + nadd;
        ib = ib + nadd;
    end
    ib = ib + 1;
end

disp(['Number Bouts : ',num2str(nb)])
