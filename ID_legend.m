function ID_legend(hID,p)

figure(hID);clf
% set up column sizes
rowCount = size(p.colorTab,1);
if isempty(p.mySpID)
    colCount = 2;
    w = [.75,.25];
else
    colCount = 3;
    w = [.3,.2,.5];
end

h = 1/(rowCount+1);% +1 to save room for title

% Set up window placement & size on screen
defaultPos=[0,0.2,0.03*colCount,.02*rowCount];
% open and setup figure window
h10 = set(hID, ...
    'NumberTitle','off', ...
    'Name','ID Legend',...
    'Units','normalized',...
    'MenuBar','none',...
    'Position',defaultPos);

% Write heading
labelStr = 'ID Color Legend';
txtPos = [0 1-h 1 h];
h10handles.headTxt = uicontrol(hID, ...
    'Style','text', ...
    'Units','normalized', ...
    'Position',txtPos, ...
    'String',labelStr, ...
    'FontUnits','normalized',...
    'FontSize',.75,...
    'FontWeight','bold');

% Put each color on a row
for iColor = 1:rowCount
    
    patchPos=[0, 1-((iColor+1)*h), w(1), h];
    h10handles.colorbox{iColor} = uicontrol(hID,...
        'Style','text',...
        'Units','normalized',...
        'Position',patchPos,...
        'BackgroundColor',p.colorTab(iColor,:),...
        'String','',...
        'FontUnits','normalized');
    
    labelPos=[w(1), 1-((iColor+1)*h), w(2), h];
    h10handles.colorLabel{iColor} = uicontrol(hID,...
        'Style','text',...
        'Units','normalized',...
        'Position',labelPos,...
        'String',num2str(iColor),...
        'FontUnits','normalized',...
        'FontSize',.5,...
        'FontWeight','bold');
    
    if ~isempty(p.mySpID)
        % if species labels are provided, add them as 3rd column
        if iColor<= length(p.mySpID)
            spPos=[sum(w(1:2)), 1-((iColor+1)*h), w(3), h];
            h10handles.spLabel{iColor} = uicontrol(hID,...
                'Style','text',...
                'Units','normalized',...
                'Position',spPos,...
                'String',p.mySpID{iColor},...
                'FontUnits','normalized',...
                'FontSize',.5,...
                'FontWeight','bold',...
                'HorizontalAlignment','left');
            
        end
    end
end
