function [zFD,zTD] = test_false_dets(XFD,k,zTD,zFD,xt,xPP,testClickIdx,ft,trueSpec,...
    flowt,fimint,fimaxt,csn,csp,loadMSP,specploton,keepers,inFileMat, hA201)
%%% Checks each test click in a session and asks if it's true or false
% moved into subroutine kf 10/4/2016


zTD(k,2) = 0;
for inxfd = 1 : zTD(k,1)
    axes(hA201(1))
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
            csnTemp = csn(testClickIdx(XFD(inxfd)),:);
            cspTemp = csp(testClickIdx(XFD(inxfd)),:);
        else
            cspTemp = inFileMat.MSP(keepers(testClickIdx(XFD(inxfd))),:);
            csnTemp = inFileMat.MSN(keepers(testClickIdx(XFD(inxfd))),:);
        end
        % make low freq part = 0
        tempSPEC = norm_spec(cspTemp, flowt,fimint,fimaxt);
        plot(ft,tempSPEC(fimint:fimaxt),'k','Linewidth',4);
        hold off;
        %
        figure(52) % add click to waveform plot in BLACK
        plot(csnTemp,'k');
        hold off;
                
        figure(52)
        hold on2
        plot(pxmsp(yell),xmpp(yell),'ro','MarkerSize',10,...
            'LineWidth',2)
        hold off
        
        figure(53)
        hold on
        plot(pxmsp,freq) % true ones in blue
        hold off
        
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
