function [yell,zFD,zID,bFlag] = brush_color(hFig,cc,zFD,zID,colorTab,t)

% brush_color.m

% Get brushed data from figure and based on specified color code modify data

% Inputs
%   hFig - Figure handle
%
%   cc - Keyboard shortcut label
%       A string with a keyboard shortcut to label data as:
%       'r' - False detection
%       'y' - Highlight selected data to show features (shown in black)
%       'g' - True detection
%       'u' - Update data according to current color brush selection. The
%       other keyboard shortcuts do not work if brush color is:
%           red - False detection
%           yellow - Highlight selected data to show features (shown in black)
%           green - True detection
%
%   zFD - An [N x 1] vector of detection times labeled as false detections,
%         where N is the number of detections.
%
%   zID - An [N x 2] matrix of detection labeled as ID detections, where
%         1st column contains the detection times and 2nd column an ID  
%         number associated with the colorTab
%
%   colorTab - Color code for classification - ID signal types
%              [191, 191,   0] type 1 green
%              [191,   0, 191] type 2 purple
%              [  0, 127,   0] type 3 dark-green
%              [  0, 191, 191] type 4 light-blue
%              [ 20,  43, 140] type 5 dark-blue
%              [218, 179, 255] type 6 pale-purple
%              [255, 214,   0] type 7 yellow
%              [222, 125,   0] type 8 orange
%              [255, 153, 199] type 9 pink
%              [153,  51,   0] type 10  brown
%
%   t - An [N x 1] vector of detection times from current window session.
%
%
%
% Output:
%
%   yell - An [N x 1] vector of indices of highlighted detection times,
%          where N is the number of detections.
%
%   zFD - An [N x 1] vector of detection times labeled as false detections,
%         where N is the number of detections.
%
%   zID - An [N x 2] matrix of detection labeled as ID detections, where
%         1st column contains the detection times and 2nd column an ID  
%         number associated with the colorTab
%
%   bFlag - Logical,track if brush is on record 

yell = [];

% Find brushed axes
hBrushAll = get(hFig,'Children');
selLine = arrayfun(@(x) hBrushAll(x).Marker == '.', 1:length(hBrushAll), 'UniformOutput', false);
if iscell(selLine)
    lenAx = cell2mat(cellfun(@length,selLine,'UniformOutput',false));
    selLine(lenAx > 1) = {false};
    selLine = cell2mat(selLine);
end
brushedLine = find(selLine,1,'last');
hBrush = hBrushAll(brushedLine);

% Get indices of brushed data
brushID = find(get(hBrush, 'BrushData'));

if ~isempty(brushID)
    
    bFlag = 1;
    
    % Get color info from brush
    h = get(brush);
    brushColor = round(h.Color*100)./100;
    
    % Get vector of dates associated with these points
    markerDates = hBrush.UserData;
    brushDate = markerDates(brushID);
    
    if ~isempty(brushDate)
        
        if isequal(brushColor,[1,0,0]) || strcmp(cc,'r')
            % Red paintbrush or 'r' = False Detections
            disp(['Number of False Detections = ',num2str(length(brushDate))])
            % Add false detections to FD matrix
            [newFD,~] = intersect(t, brushDate);
            zFD = [zFD; newFD]; 
            
        elseif isequal(brushColor,[0,1,0]) || strcmp(cc,'g')
            % Green paintbrush or 'g' = True Detections
            disp(['Number of Detections Selected = ',num2str(length(brushDate))])
            if exist('zFD','var')
                
                % Make sure you're not flagging something outside this session
                [newTD,~] = intersect(t, brushDate);
                [zFD,iC] = setdiff(zFD,newTD);
                
                % Clear brushed set from ID set
                if ~isempty(zID) 
                    [~,zIDkeep] = setdiff(zID(:,1),newTD);
                    zID = zID(zIDkeep,:);
                end
                
                % Clear brushed set from FD set
                if ~isempty(zFD) 
                    [~,zFDkeep] = setdiff(zFD(:,1),newTD);
                    zFD = zFD(zFDkeep,:);
                end
                disp(['Remaining Number of False Detections = ',num2str(length(iC))])
            end    
            
        elseif isequal(brushColor,[1,1,0]) %|| strcmp(cc,'y')
            % Yellow paintbrush or 'y' = Highlight Detections
            disp(['Start time selected data: ',datestr(brushDate(1),'dd-mm-yyyy HH:MM:SS.FFF')]);
            [~,yell] = intersect(t, brushDate);
            
        else
            % Find color code to assign ID signal type
            tfCompare = find(sum(bsxfun(@eq,colorTab,brushColor),2)==3);
            if ~isempty(tfCompare)
                specID = tfCompare;
                
                % Write to ID file
                disp(['Number of ID Detections = ',num2str(length(brushDate))])
                [newIDtimes,~] = intersect(t, brushDate);
                
                % check if already brushed
                if ~isempty(zID)
                    [~,oldIDtimes] = intersect(zID(:,1), newIDtimes);
                    if ~isempty(oldIDtimes)
                        % remove any old labels for these times
                        zID(oldIDtimes, :) = [];
                    end
                end

                spLabels = ones(size(newIDtimes)).*specID;
                newID = [newIDtimes,spLabels];
                zID = [zID;newID];
            else
                bFlag = 0;
            end
        end
    end
    
    % Remove any ID from FD
    if ~isempty(zID) 
        [~,zIDkeep2] = setdiff(zID(:,1),zFD);
        zID = zID(zIDkeep2,:);
    end
else
    bFlag = 0;
end
