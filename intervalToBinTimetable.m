function T = intervalToBinTimetable(Start,End,p) 


Startnum = datenum(Start);
timevec = datevec(Startnum);
stv = size(timevec);
binStartEffort = datetime([timevec(:,1:4), floor(timevec(:,5)/p.binDur)*p.binDur, ...
    zeros(stv(1),1)]);

Endnum = datenum(End);
timevec = datevec(Endnum);
binEndEffort = datetime([timevec(:,1:4), floor(timevec(:,5)/p.binDur)*p.binDur, ...
    zeros(stv(1),1)]);

tbin = [];
for i = 1: length(Start)
    tbin = [tbin;(binStartEffort(i):minutes(5):binEndEffort(i))']; 
end

T = timetable(tbin);
T.bin = ones(height(T),1);