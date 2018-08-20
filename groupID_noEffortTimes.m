% select directory where ship files are located
shipDir = 'E:\metadata_reduced';
shipTimesDir = 'H:\ShipTimes'; % directory where to save ship times .mat files

IDDir = 'H:\TPWS';
IDTimesDir = 'H:\IDTimes'; % directory where to save ID times .mat files
maxDetEdit = 4; % number of TPWS folders (i.e. TPWS4 is 4)
saveTable = 'H:\Pm_Effort.xls'; % directory where to save excel file with effort times

%% write ship file times
%run join_IDs to group noe effort times
%join_IDs
writeSFilesTimes(shipDir,shipTimesDir);
writeIDTimes(IDDir,IDTimesDir,maxDetEdit);

%% Get a list of all the files in the start directory
shipList = cellstr(ls(shipTimesDir));
shipfiles = shipList(3:end); % exclude dots

IDList = cellstr(ls(IDTimesDir));
IDfiles = IDList(3:end); % exclude dots

%% get start end dates of disks
[edgeffort,latLongs, depl, site] = gofmx_dates;

strdepl = num2str(depl,'%02d');
locations = strcat(site,strdepl);
nameLoc = unique(locations);
%% Extract effort times
effTable = table();
for n = 1:length(nameLoc)
    ieff = find(contains(locations,nameLoc(n)));
    loc = nameLoc{n};
    s = loc(1:2);
    % load ship times
    iS = contains(shipfiles,loc);
    iID = contains(IDfiles,loc);
    if sum(iS) == 0
        loc = strcat(loc(1:2),'_',loc(3:4)); % some sites have this case
        iS = contains(shipfiles,loc);
        iID = contains(IDfiles,loc);
    end
    
    if any([iS;iID]) % only if there is ship and ID files
        if any(iS)
            ship = load(fullfile(shipTimesDir,shipfiles{iS}));
        else
            ship.times = [];
        end
        
        if any(iID)
            ID = load(fullfile(IDTimesDir,IDfiles{iID}));
        else 
            ID.times = []; 
        end
        
        % check if multiple edgeffort times (some disk are corrupted)
        if length(ieff) > 1
            addInterval = [edgeffort(ieff(1:end-1),2) + datenum(0,0,0,0,0,1/1000), edgeffort(ieff(2:end),1)];
            % group all times
            times = [ship.times; ID.times; addInterval];
        else
            % group all times
            times = [ship.times; ID.times];
        end
        
        % group overlapping intervals of no effort
        times = groupoverlaps(times);

        %{
        figure
        nline  = repmat((1:size(times,1))',1,2);
        plot(times',nline','bo-')
        title(loc)
        hold on
        nline  = repmat((1:size(ship.times,1))',1,2);
        plot(ship.times',nline','r.-')
        nline  = repmat((1:size(ID.times,1))',1,2);
        plot(ID.times',nline','g.-')
        xlim([edgeffort(ieff(1),1) edgeffort(ieff(end),2)])
        %}
        
        % get times between no effort intervals
        effort = [];
        effort(:,1) = [edgeffort(ieff(1),1); times(:,2) + datenum(0,0,0,0,0,1/1000)]; % add a ms
        effort(:,2) = [times(:,1) - datenum(0,0,0,0,0,1/1000); edgeffort(ieff(end),2)]; % extract a ms
        
        %{
        figure
        nline  = repmat((1:size(effort,1))',1,2);
        plot(effort',nline','bo-')
        title(loc)
        hold on
        nline  = repmat((1:size(times,1))',1,2);
        plot(times',nline','r.-') 
        %}
        
        S = cell(length(effort),1);
        LOC = cell(length(effort),1);
        lat = cell(length(effort),1);
        log = cell(length(effort),1);
        
        S(:) = {s};
        LOC(:) = {loc};
        lat(:) = {latLongs(n,1)};
        log(:) = {latLongs(n,2)};
        
        effTable = [effTable; table(S,LOC,lat,log,m2xdate(effort(:,1)),m2xdate(effort(:,2)))];
    end
    
end
effTable.Properties.VariableNames = {'Sites','Deployments','Latitude','Longitude','StartEffort','EndEffort'};

% save effort times in excel file
writetable(effTable,saveTable)

disp('Effort times saved')
