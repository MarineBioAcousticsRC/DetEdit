% select directory where ship files are located
shipDir = 'E:\metadata';
shipTimesDir = 'E:\ShipTimes'; % directory where to save ship times .mat files

IDDir = 'E:\TPWS\TPWS2\';
IDTimesDir = 'E:\IDTimes'; % directory where to save ID times .mat files

saveTable = 'E:\Pm_Effort.xls'; % directory where to save excel file with effort times
%% write ship file times
% writeSFilesTimes(shipDir,shipTimesDir);
% writeIDTimes(IDDir,shipTimesDir);

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
    
    if sum([iS;iID]) ~= 0
        ship = load(fullfile(shipTimesDir,shipfiles{iS}));
        ID = load(fullfile(IDTimesDir,IDfiles{iID}));
        
        % check if multiple multiple edgeffort times (some disk are corrupted)
        if length(ieff) > 1
            addInterval = [edgeffort(ieff(1:end-1),2) + datenum(0,0,0,0,0,1), edgeffort(ieff(2:end),1)];
            % group all times
            times = [ship.times; ID.times; addInterval];
        else
            % group all times
            times = [ship.times; ID.times];
        end
        
        % group overlapping intervals of no effort
        times = groupoverlaps(times);
        
        % get times between no effort intervals
        effort = [];
        effort(:,1) = [edgeffort(ieff(1),1); times(:,2) + datenum(0,0,0,0,0,1)]; % add a second
        effort(:,2) = [times(:,1) - datenum(0,0,0,0,0,1); edgeffort(ieff(end),2)]; % extract a second
        
        figure
        nline  = repmat((1:size(effort,1))',1,2);
        plot(effort',nline','o-')
        title(loc)
        
        S = cell(length(effort),1);
        LOC = cell(length(effort),1);
        lat = cell(length(effort),1);
        log = cell(length(effort),1);
        
        S(:) = {s};
        LOC(:) = {loc};
        lat(:) = {latLongs(n,1)};
        log(:) = {latLongs(n,2)};
        
        effTable = [effTable; table(S,LOC,lat,log,effort(:,1),effort(:,2))];
    end
    
end
effTable.Properties.VariableNames = {'Sites','Deployments','Latitude','Longitude','StartEffort','EndEffort'};

% save effort times in excel file
writetable(effTable,saveTable)

disp('Effort times saved')
