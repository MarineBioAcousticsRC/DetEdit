function [yell,zFD,zID,bFlag] = brush_color(hFig,cc,zFD,zID,colorTab,t)
% get brushed data and figure out what to do based on color

yell = [];

hBrushAll = get(hFig,'Children');

% find axes that correspond to the selected points
iloop = 0;
idx = 0;
while iloop == 0
    idx = idx + 1;
    iloop = hBrushAll(idx).Marker == '.';
end
hBrush = hBrushAll(idx);

% get indices of brushed data
brushID = find(get(hBrush, 'BrushData'));

if ~isempty(brushID)
    
    bFlag = 1;
    
    % get color info from brush
    h = get(brush);
    brushColor = round(h.Color*100)./100;
    % get vector of dates associated with these points
    markerDates = hBrush.UserData;
    brushDate = markerDates(brushID);
    
    if ~isempty(brushDate)
        if isequal(brushColor,[1,0,0]) || strcmp(cc,'r')
            % Red paintbrush = False Detections
            disp(['Number of False Detections = ',num2str(length(brushDate))])
            % make sure you're not flagging something outside this session
            [newFD,~] = intersect(t, brushDate);
            zFD = [zFD; newFD];   % cummulative False Detection matrix
        elseif isequal(brushColor,[1,1,0]) || strcmp(cc,'y')
            % yellow paintbrush = give time of detection
            disp(['       Year              Month          Day           Hour', ...
                '          Min          Sec']);
            disp(['Datevec ',num2str(datevec(brushDate(1)))]);
            [~,yell] = intersect(t, brushDate);
            
        elseif isequal(brushColor,[0,5,0]) || strcmp(cc,'g')||...
            isequal(brushColor,[0,0,0]) || strcmp(cc,'i')
        
            disp(['Number of Detections Selected = ',num2str(length(brushDate))])
            if exist('zFD','var')
                % make sure you're not flagging something outside this session
                [newTD,~] = intersect(t, brushDate);
                [zFD,iC] = setdiff(zFD,newTD);
                if ~isempty(zID) % clear brushed set from ID set
                    [~,zIDkeep] = setdiff(zID(:,1),newTD);
                    zID = zID(zIDkeep,:);
                end
                if ~isempty(zFD) % clear brushed set from FD set
                    [~,zFDkeep] = setdiff(zFD(:,1),newTD);
                    zFD = zFD(zFDkeep,:);
                end
                disp(['Remaining Number of False Detections = ',num2str(length(iC))])
                % save(fn2,'zFD')
            end
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
            else
                bFlag = 0;
            end
        end
    end
    if ~isempty(zID) % remove any ID from FD
        [~,zIDkeep2] = setdiff(zID(:,1),zFD);
        zID = zID(zIDkeep2,:);
    end
else
    bFlag = 0;
end
