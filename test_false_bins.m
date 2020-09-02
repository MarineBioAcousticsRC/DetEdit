function zTD = test_false_bins(zTD)
%%% Checks a window in a session and asks if it's true or false
% moved into subroutine kf 10/4/2016
global dHANDLES dPARAMS p

dPARAMS.lab = input('Enter label to test: ');
thisLabField = sprintf('Label_%d', dPARAMS.lab);

if zTD(dPARAMS.k).(thisLabField)(1) == 0
    disp(['No test clicks for ',thisLabField, ' in this session']);
    dPARAMS.k = dPARAMS.k+1;
    return
end
if zTD(dPARAMS.k).(thisLabField)(2) == -1
    disp(['WARNING: False click test not yet completed for ',thisLabField,', cannot carry out false bin test']);
    dPARAMS.k = dPARAMS.k+1;
    return
end

% (Might be nice to highlight the clicks in question on fig 51 also?)
% Plot all test clicks in session
hold(dHANDLES.LTSAsubs(1),'on')
for inxfd = 1 : zTD(dPARAMS.k).(thisLabField)(1)
    testTime = dPARAMS.xt{1,dPARAMS.lab}(inxfd);
    xH201a = plot(dHANDLES.LTSAsubs(1),testTime,dPARAMS.xPP{1,dPARAMS.lab}(inxfd),...
        'ro','MarkerSize',10,'UserData',dPARAMS.xt{1,dPARAMS.lab}(inxfd));
end
hold(dHANDLES.LTSAsubs(1),'off')
hold(dHANDLES.LTSAsubs(3),'on')
for inxfd = 1 : zTD(dPARAMS.k).(thisLabField)(1)
    testTime = dPARAMS.xt{1,dPARAMS.lab}(inxfd);
    clickInBoutIdx = find(dPARAMS.t==testTime);
    
    xH201b = plot(dHANDLES.LTSAsubs(3),testTime,dPARAMS.dt(clickInBoutIdx),...
        'ro','MarkerSize',10,'UserData',dPARAMS.xt{1,dPARAMS.lab}(inxfd));
end
hold(dHANDLES.LTSAsubs(3),'off')
% figure(201)
% subplot(3,1,1)  % top panel RL vs Time
% hold(dHANDLES.LTSAsubs(1),'on')
% for inxfd = 1 : zTD(dPARAMS.k).(thisLabField)(1)
%     plot(dPARAMS.xt{1,dPARAMS.lab}(inxfd),dPARAMS.xPP{1,dPARAMS.lab}(inxfd),...
%         'ro','MarkerSize',10,'UserData',dPARAMS.xt{1,dPARAMS.lab}(inxfd));
% end
% hold off
testBins = length(vertcat(dPARAMS.binCX{:,dPARAMS.lab}));
prompt = (['SESSION:',num2str(dPARAMS.k),', # Test Detect Bins: ',...
    num2str(testBins),', # False Bins: ']);
pzTD = input(prompt);

% if the input is out of bounds, say so, and ask again
while ~(pzTD >= 0 && pzTD <= length(dPARAMS.binCX))
    disp(' Number False Bins Out of Bounds')
    pzTD = input(prompt);
end

zTD(dPARAMS.k).(thisLabField)(4) = pzTD;
zTD(dPARAMS.k).(thisLabField)(3) = testBins;
dPARAMS.k = dPARAMS.k+1;


