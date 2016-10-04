function [zFD,zTD] = test_false_dets(XFD,k,zTD,zFD,xt,xPP,ixfd,ft,trueSpec,...
    flowt,fimint,fimaxt,csn,csp,loadMSP,specploton)
%%% Checks each test click in a session and asks if it's true or false
% moved into subroutine kf 10/4/2016

if (~isempty(XFD))
    zTD(k,2) = 0;
    for inxfd = 1 : zTD(k,1)
        figure(201)
        subplot(3,1,1)  % top panel RL vs Time
        hold on
        testTimes = xt(inxfd);
        plot(testTimes,xPP(inxfd),...
            'ro','MarkerSize',10,'UserData',testTimes);
        hold off;
        disp(['Showing #: ',num2str(inxfd),' click. Press ''z'' to reject']);
        if (specploton == 1)
            figure(50)  % add click to spec plot in BLACK
            % meanspec;
            plot(ft,trueSpec(fimint:fimaxt),'Linewidth',4);
            hold on;
            if loadMSP
                csnTemp = csn(ixfd(XFD(inxfd)),:);
                cspTemp = csp(ixfd(XFD(inxfd)),:);
            else
                cspTemp = inFileMat.MSP(keepers(ixfd(XFD(inxfd))),:);
                csnTemp = inFileMat.MSN(keepers(ixfd(XFD(inxfd))),:);
            end
            % make low freq part = 0
            tempSPEC = norm_spec(cspTemp, flowt,fimint,fimaxt);
            plot(ft,tempSPEC(fimint:fimaxt),'k','Linewidth',4);
            hold off;
            %
            figure(52) % add click to waveform plot in BLACK
            plot(csnTemp,'k');
            hold off;
            pause
            cc = get(52,'CurrentCharacter');
            if (strcmp(cc,'z'))
                zTD(k,2) = zTD(k,2) + 1;
                zFD = [zFD; xt(inxfd)]; % add to FD
            end
        end
    end
    disp([' Tested: ',num2str(zTD(k,1)),' False: ',...
        num2str(zTD(k,2))]);
end