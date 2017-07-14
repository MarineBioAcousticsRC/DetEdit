function [yell,zFD,zID,zMD,bFlag] = brush_color(hFig,cc,zFD,zID,zMD,colorTab,t)
% get brushed data and figure out what to do based on color

yell = [];

hBrush = findall(hFig,'tag','Brushing');

% get x and x info from brushed data
brushDataX = get(hBrush, {'Xdata'});
brushDataY = get(hBrush, {'Ydata'});
% don't understand why nan values appear, create index of valid points here

if ~isempty(brushDataX)
    brushID = ~isnan(brushDataX{1,1});
    
    bFlag = 1;
    brushDataX = brushDataX{1,1}(brushID);
    brushDataY = brushDataY{1,1}(brushID);
    
    % get color info
    brushColor = get(hBrush, {'Color'});
    brushColor = round(brushColor{1,1}.*100)./100;
    % get vector of dates associated with these points
    markerDates = get(findall(hFig),'UserData');
    filledMarkers = markerDates(~cellfun('isempty',markerDates));
    brushDate = filledMarkers{end}(brushID);
    
    if ~isempty(brushDate)
        if isequal(brushColor,[1,0,0]) || strcmp(cc,'r');
            % Red paintbrush = False Detections
            disp(['Number of False Detections = ',num2str(length(brushDataX))])
            % make sure you're not flagging something outside this session
            [newFD,~] = intersect(t, brushDate);
            zFD = [zFD; newFD];   % cummulative False Detection matrix
        elseif isequal(brushColor,[1,1,0]) || strcmp(cc,'y');
            % yellow paintbrush = give time of detection
            disp(['       Year              Month          Day           Hour', ...
                '          Min          Sec']);
            disp(['Datevec ',num2str(datevec(brushDate(1)))]);
            [~,yell] = intersect(t, brushDate);
            
        elseif isequal(brushColor,[0,0,0]) || strcmp(cc,'i');
            % black paintbrush = % set back to true
            disp(['Number of Detections Selected = ',num2str(length(brushDate(1)))])
            if exist('zFD','var')
                % make sure you're not flagging something outside this session
                [newTD,~] = intersect(t, brushDate);
                [zFD,iC] = setdiff(zFD,newTD);
                if ~isempty(zID) % clear brushed set from ID set
                    [~,zIDkeep] = setdiff(zID(:,1),newTD);
                    zID = zID(zIDkeep,:);
                end
                if ~isempty(zMD) % clear brushed set from ID set
                    [~,zMDkeep] = setdiff(zMD(:,1),newTD);
                    zMD = zMD(zMDkeep,:);
                end
                if ~isempty(zFD) % clear brushed set from FD set
                    [~,zFDkeep] = setdiff(zFD(:,1),newTD);
                    zFD = zFD(zFDkeep,:);
                end
                disp(['Remaining Number of False Detections = ',num2str(length(iC))])
                % save(fn2,'zFD')
            end
        elseif isequal(brushColor,[0,1,0]) || strcmp(cc,'m');
            % green paintbrush = misidentified Detections
            disp(['Number of Mid-ID Detections = ',num2str(length(brushDataX))])
            % make sure you're not flagging something outside this session
            [newMD,~] = intersect(t, brushDate);
            % remove new MD from ID
            if ~isempty(zID) % clear brushed set from ID set
            	[~,zIDkeep] = setdiff(zID(:,1),newMD);
            	zID = zID(zIDkeep,:);
            end
            % lowest priority, could alternatively remove things from 
            zMD = [zMD; newMD];   % cummulative False Detection matrix
        else
            % find the brush color
            tfCompare = find(sum(bsxfun(@eq,colorTab,brushColor),2)==3);
            if ~isempty(tfCompare)
                specID = tfCompare;
                % write to ID file
                disp(['Number of ID Detections = ',num2str(length(brushDate))])
                [newIDtimes,~] = intersect(t, brushDate);
                spLabels = ones(size(newIDtimes)).*specID;
                newID = [newIDtimes,spLabels];
                zID = [zID;newID];
                % remove new ID from MD
                if ~isempty(zMD) % clear brushed set from ID set
                    [~,zMDkeep] = setdiff(zMD,zID(:,1));
                    zMD = zMD(zMDkeep,:);
                end
            else
                bFlag = 0;
            end
        end
    end
    if ~isempty(zID) % remove any ID from FD
        [~,zIDkeep2] = setdiff(zID(:,1),zFD);
        zID = zID(zIDkeep2,:);
    end
    if ~isempty(zMD) % remove any MD from FD
        [~,zMDkeep2] = setdiff(zMD(:,1),zFD);
        zMD = zMD(zMDkeep2,:);
    end
else
    bFlag = 0;
end
