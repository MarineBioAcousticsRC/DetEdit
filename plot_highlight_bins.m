function plot_highlight_bins(kind)

% plot_highlight_bins
%
% Shared by figure54 (kind = 'spectra'), figure55 (kind = 'waveform') and
% figure56 (kind = 'ici').
%
% Average spectra (or, for waveforms, the single most typical click: the
% one closest to the median waveform) of the yellow-brushed clicks,
% split into received level bins (dB pp from the top LTSA panel), one panel
% per bin in a near-square grid, loudest bin first (top left).
%
% In each panel, clicks are averaged per predicted label (zID col 2) and
% drawn in that label's color (unlabeled = default blue, false = red).
% When p.plotCorrect_incorrect is on, each label is further split by
% ground truth (zID col 4):
%   solid  - ground truth matches the prediction
%   dotted - ground truth differs from the prediction (drawn on top); one
%            line per ground truth label, each with its own dotted pattern /
%            width (still in the predicted label's color), and the legend
%            says which ground truth it was
%   dashed - no ground truth for those clicks
%
% The window only opens once something is highlighted; it is cleared when
% nothing is highlighted.

global dPARAMS dHANDLES p

binWidth = 5; % dB per bin
waveSpacing = 1.2; % vertical offset between stacked waveforms / ICI histograms
iciNBins = 25; % ICI histogram bins between 0 and p.dtHi

switch kind
    case 'spectra'
        figField = 'specBinFig'; figNum = 54; dataSet = 'cspJ';
        figName = 'Highlighted spectra by received level';
        whatText = 'spectra';
        xText = 'Frequency (kHz)'; yText = 'Normalized RL';
    case 'waveform'
        figField = 'waveBinFig'; figNum = 55; dataSet = 'csnJ';
        figName = 'Highlighted waveforms by received level';
        whatText = 'median-click waveforms';
        xText = sprintf('Samples (%0.1fms @ %0.0fkHz)',...
            size(dPARAMS.csnJ,2)/p.sampleRate,p.sampleRate);
        yText = 'Normalized amplitude';
    case 'ici'
        figField = 'iciBinFig'; figNum = 56; dataSet = 'dt';
        figName = 'Highlighted ICI by received level';
        whatText = 'ICI histograms';
        xText = 'Time between detections (s)';
        yText = 'Normalized count';
    otherwise
        error('plot_highlight_bins: kind must be ''spectra'', ''waveform'' or ''ici''')
end

figExists = isfield(dHANDLES,figField) && ~isempty(dHANDLES.(figField)) && ...
    ishghandle(dHANDLES.(figField));
haveSelection = ~isempty(dPARAMS.yell) && p.specploton && ~isempty(dPARAMS.(dataSet));

if ~haveSelection
    if figExists
        clf(dHANDLES.(figField))
        uicontrol(dHANDLES.(figField),'Style','text','Units','normalized',...
            'Position',[.1 .45 .8 .1],'FontSize',11,...
            'String',['Highlight clicks (yellow brush) to see their ',whatText,' by received level']);
    end
    return
end

if ~figExists
    oldFig = get(0,'CurrentFigure');
    dHANDLES.(figField) = figure(figNum);
    set(dHANDLES.(figField),'name',figName,...
        'NumberTitle','off','KeyPressFcn',@keyAction,...
        'KeyReleaseFcn',@keyRelease)
    if ~isempty(oldFig) && ishghandle(oldFig)
        set(0,'CurrentFigure',oldFig) % don't steal focus from the main windows
    end
end
fig = dHANDLES.(figField);
clf(fig)

yell = dPARAMS.yell(:);
nClicks = length(yell);
RLy = dPARAMS.RL(yell); RLy = RLy(:);

% --- what each highlighted click is: predicted code + ground truth match
% predCode: label code, 0 = unlabeled, 99 = false
% matchType: 1 = correct, 2 = no ground truth, 3 = wrong, 0 = n/a
predCode = zeros(nClicks,1);
matchType = zeros(nClicks,1);
gtSel = zeros(nClicks,1); % ground truth code, only kept for wrong ones
isFD = false(nClicks,1);
if dPARAMS.ff2
    isFD = ismember(yell,dPARAMS.K2);
end
predCode(isFD) = 99;
if dPARAMS.ff3
    [isID,idLoc] = ismember(yell,dPARAMS.K3);
    isID(isID) = dPARAMS.labelConfIdx(idLoc(isID)); % low confidence = unlabeled
    isID = isID & ~isFD;
    predCode(isID) = dPARAMS.spCodeSet(idLoc(isID));
    if dPARAMS.showGT
        gt = dPARAMS.gtCode(idLoc(isID));
        mt = 3*ones(size(gt));
        mt(gt == predCode(isID)) = 1;
        mt(isnan(gt)) = 2;
        matchType(isID) = mt;
        gtWrong = zeros(size(gt));
        gtWrong(mt == 3) = gt(mt == 3);
        gtSel(isID) = gtWrong;
    end
end

% --- level bins, highest on top
binLow = floor(RLy/binWidth)*binWidth;
bins = sort(unique(binLow),'descend');
nBins = length(bins);

defaultColor = get(groot,'defaultAxesColorOrder');
lineStyles = {'-','-','--',':'}; % index = matchType+1
% wrong lines: within each predicted label in a panel, the 1st ground
% truth gets pattern 1, the 2nd pattern 2, ... (the legend names the ground
% truth). Solid (correct) and thick dashed (no ground truth) aren't used.
wrongStyles = {':',1.2; '-.',1.2; ':',2.4; '-.',2.4; '--',0.8};
styleNames = {'','',' (no gt)',' (wrong)'};

% layout: near-square grid of panels (each keeps the wide shape of the
% main spectra/waveform plots), filled left to right, top to bottom, loudest first
left = .07; right = .02; bottom = .08; top = .90; hGap = .05; vGap = .09;
targetAspect = 1.3; % width/height of the main spectra plot's axes
figUnits = get(fig,'Units'); set(fig,'Units','pixels');
figPos = get(fig,'Position'); set(fig,'Units',figUnits);
% pick the grid that gives the biggest panels once each panel is held to
% the target shape (centered in its grid cell)
bestArea = -inf;
for tryCols = 1:nBins
    tryRows = ceil(nBins/tryCols);
    [w,h] = panel_size(tryCols,tryRows);
    if w > 0 && h > 0 && w*h > bestArea
        bestArea = w*h; nCols = tryCols; nRows = tryRows;
    end
end
cellW = (1 - left - right - hGap*(nCols-1))/nCols;
cellH = (top - bottom - vGap*(nRows-1))/nRows;
[axW,axH] = panel_size(nCols,nRows);
annotation(fig,'textbox',[0 .93 1 .06],'String',...
    sprintf('%d highlighted clicks: normalized %s by received level (%g dB bins)',...
    nClicks,whatText,binWidth),'EdgeColor','none','HorizontalAlignment','center',...
    'FontWeight','bold','FontSize',10);

for iB = 1:nBins
    iRow = ceil(iB/nCols);
    iCol = iB - (iRow-1)*nCols;
    cellX = left + (iCol-1)*(cellW+hGap);
    cellY = top - iRow*cellH - (iRow-1)*vGap;
    ax = axes('Parent',fig,'Units','normalized','Position',...
        [cellX + (cellW-axW)/2, cellY + (cellH-axH)/2, axW, axH]);
    hold(ax,'on')
    inBin = binLow == bins(iB);

    % groups of (prediction, match type, ground truth if wrong) in this
    % bin; wrong ones last (on top)
    groups = unique([predCode(inBin), matchType(inBin), gtSel(inBin)],'rows');
    [~,order] = sort(groups(:,2) == 3);
    groups = groups(order,:);
    hLines = []; lineNames = {};
    nWrongSoFar = containers.Map('KeyType','double','ValueType','double');

    for iG = 1:size(groups,1)
        thisCode = groups(iG,1);
        thisMatch = groups(iG,2);
        thisGT = groups(iG,3);
        idx = inBin & predCode == thisCode & matchType == thisMatch & gtSel == thisGT;
        if thisCode == 0
            thisColor = defaultColor(1,:);
        elseif thisCode == 99
            thisColor = [1,0,0];
        else
            thisColor = p.colorTab(thisCode,:);
        end
        thisName = label_name(thisCode);
        if strcmp(kind,'spectra')
            xData = dPARAMS.ft;
            yData = norm_spec_simple(dPARAMS.cspJ(yell(idx),:),...
                dPARAMS.fimint,dPARAMS.fimaxt);
        else
            % each waveform on its own row (first group on top), like
            % the main waveform window, so they don't overlap
            % show ONE real click: the one whose normalized waveform is
            % closest to the group's sample-by-sample median waveform
            grpWav = norm_wav(dPARAMS.csnJ(yell(idx),:));
            medWav = median(grpWav,1);
            [~,iMed] = min(sum((grpWav - medWav).^2,2));
            yData = grpWav(iMed,:) + waveSpacing*(size(groups,1) - iG);
            xData = 1:length(yData);
        end
        if strcmp(kind,'ici')
            % histogram of each click's time to the next detection (same
            % ICI as the LTSA bottom panel), scaled to a peak of 1 and
            % stacked one group per row like the waveforms
            clk = yell(idx);
            clk = clk(clk <= numel(dPARAMS.dt));
            edges = linspace(0,p.dtHi,iciNBins+1);
            cnt = histc(min(dPARAMS.dt(clk(:)),p.dtHi),edges); % values above dtHi go in the last bin
            cnt = cnt(:).';
            cnt(end-1) = cnt(end-1) + cnt(end); % histc's last bin is "== last edge"
            cnt = cnt(1:end-1);
            if max(cnt) > 0
                cnt = cnt/max(cnt);
            end
            [xData,yData] = stairs(edges,[cnt,cnt(end)]);
            xData = xData(:).'; 
            yData = yData(:).' + waveSpacing*(size(groups,1) - iG);
        end
        if thisMatch == 3
            % wrong: dotted-type line; k-th wrong ground truth for this
            % predicted label in this panel gets the k-th pattern/width
            if isKey(nWrongSoFar,thisCode)
                nWrongSoFar(thisCode) = nWrongSoFar(thisCode) + 1;
            else
                nWrongSoFar(thisCode) = 1;
            end
            iStyle = mod(nWrongSoFar(thisCode)-1,size(wrongStyles,1))+1;
            hLines(end+1) = plot(ax,xData,yData,'Color',thisColor,...
                'LineStyle',wrongStyles{iStyle,1},...
                'LineWidth',wrongStyles{iStyle,2}); %#ok<AGROW>
            lineNames{end+1} = sprintf('%s (n=%d, gt=%s)',thisName,sum(idx),...
                label_name(thisGT)); %#ok<AGROW>
        else
            hLines(end+1) = plot(ax,xData,yData,'Color',thisColor,...
                'LineStyle',lineStyles{thisMatch+1},'LineWidth',2); %#ok<AGROW>
            lineNames{end+1} = sprintf('%s%s (n=%d)',thisName,styleNames{thisMatch+1},sum(idx)); %#ok<AGROW>
        end
    end
    hold(ax,'off')

    if strcmp(kind,'spectra')
        xlim(ax,[p.fLow,p.fHi])
    elseif strcmp(kind,'ici')
        xlim(ax,[0,p.dtHi])
    else
        xlim(ax,[1,size(dPARAMS.csnJ,2)])
    end
    grid(ax,'on')
    if strcmp(kind,'spectra')
        ylim(ax,[0,1.05]); set(ax,'YTick',[0,.5,1])
    else
        ylim(ax,[-.05, waveSpacing*(size(groups,1)-1) + 1.05])
        set(ax,'YTick',[]) % stacked rows: amplitude ticks don't apply
    end
    title(ax,sprintf('%g-%g dB',bins(iB),bins(iB)+binWidth))
    legend(ax,hLines,lineNames,'Location','northeast','FontSize',7)
    if iCol == 1
        ylabel(ax,yText)
    end
    if iB + nCols > nBins % nothing below this panel
        xlabel(ax,xText)
    end
end

    function name = label_name(code)
        % display name for a label code (legend text)
        if code == 0
            name = 'Unlabeled';
        elseif code == 99
            name = 'False';
        elseif iscell(p.mySpID) && code <= numel(p.mySpID) && ...
                ischar(p.mySpID{code}) && ~isempty(p.mySpID{code})
            name = p.mySpID{code};
        else
            name = num2str(code);
        end
    end

    function [w,h] = panel_size(nC,nR)
        % panel size (normalized) for an nC x nR grid, held to targetAspect
        w = (1 - left - right - hGap*(nC-1))/nC;
        h = (top - bottom - vGap*(nR-1))/nR;
        if (w*figPos(3))/(h*figPos(4)) > targetAspect
            w = targetAspect*h*figPos(4)/figPos(3); % too wide: narrow it
        else
            h = w*figPos(3)/(targetAspect*figPos(4)); % too tall: shorten it
        end
    end
end
