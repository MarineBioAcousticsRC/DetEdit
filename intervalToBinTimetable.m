function T = intervalToBinTimetable(Start,End,p) 

% intervalToBinTimetable.m

% Simple script to get bin times of specified time from given intervals

Startnum = datenum(Start);
timevec = datevec(Startnum);
binStartEffort = datetime([timevec(:,1:4), floor(timevec(:,5)/p.binDur)*p.binDur, ...
    zeros(length(timevec),1)]);

Endnum = datenum(End);
timevec = datevec(Endnum);
binEndEffort = datetime([timevec(:,1:4), floor(timevec(:,5)/p.binDur)*p.binDur, ...
    zeros(length(timevec),1)]);

tbin = [];
effortBin = [];
effortSec = [];
for i = 1: length(Start)
    currInterval = (binStartEffort(i):minutes(5):binEndEffort(i))';
    if ~isempty(currInterval)
    bin = ones(length(currInterval),1);
    sec = ones(length(currInterval),1)*p.binDur*60;
    % check if first and last bin are complete
    if (binStartEffort(i)-Start(i)) ~= 0
       bin(1) =  (minutes(5)-(Start(i)-binStartEffort(i)))/minutes(5);
       sec(1) = seconds(minutes(5)-(Start(i)-binStartEffort(i)));
    end
    
    if End(i) > binEndEffort(i)
       bin(end) = (End(i) - binEndEffort(i))/minutes(5);
       sec(end) = seconds(End(i) - binEndEffort(i));
    end
    tbin = [tbin;currInterval];
    effortBin = [effortBin;bin];
    effortSec = [effortSec;sec];
    
    end
end

T = timetable(tbin,effortBin,effortSec);