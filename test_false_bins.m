function zTD = test_false_bins(k,zTD,xt,xPP,binCX)
%%% Checks a window in a session and asks if it's true or false
% moved into subroutine kf 10/4/2016

% (Might be nice to highlight the clicks in question on fig 51 also?)
% Plot all test clicks in session
figure(201)
subplot(3,1,1)  % top panel RL vs Time
hold on
for inxfd = 1 : zTD(k,1)
    plot(xt(inxfd),xPP(inxfd),...
        'ro','MarkerSize',10,'UserData',xt(inxfd));
end
hold off
prompt = (['SESSION:',num2str(k),' #Test Detect Bins: ',...
    num2str(length(binCX)),' #False Bins: ']);
pzTD = input(prompt);

% if the input is out of bounds, say so, and ask again
while ~(pzTD >= 0 && pzTD <= length(binCX))
    disp(' Number False Bins Out of Bounds')
    pzTD = input(prompt);
end

zTD(k,4) = pzTD;
zTD(k,3) = length(binCX);

