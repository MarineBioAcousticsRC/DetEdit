function LabCertainty_Eval
global dHANDLES dPARAMS p fNameList

dPARAMS.lab = [];

A6 = exist(strrep(fNameList.TPWS,'TPWS1','labCert'));
if (A6 ~= 2)
    labelCertainty = {};
    RL = {};
    for i = 1:length(p.mySpID)
        labelCertainty{1,i} = [];
        RL{1,i} = [];
        RL{2,i} = [];
    end
    save(strrep(fNameList.TPWS,'TPWS1','labCert'),'labelCertainty','RL','p');    % create new labCert file
    disp(' Making new labCert file');
else
    disp(' Loading labCert file');
    load(strrep(fNameList.TPWS,'TPWS1','labCert'),'labelCertainty','RL')
end

while isempty(dPARAMS.lab)
    dPARAMS.lab = input('Enter label to test: ');
    if ~isempty(dPARAMS.lab) & ~ismember(dPARAMS.lab,1:length(p.mySpID))
        disp('WARNING: Entry not allowed');
        dPARAMS.lab = [];
    end
end

% check for test bins in this session
tbidx = dPARAMS.binTimes(dPARAMS.ixtb{1,dPARAMS.lab}) > dPARAMS.sb(dPARAMS.k)-datenum([0,0,0,0,p.binDur,0])...
    & dPARAMS.binTimes(dPARAMS.ixtb{1,dPARAMS.lab}) < dPARAMS.eb(dPARAMS.k);
if any(tbidx)
    
    BT = dPARAMS.binTimes(dPARAMS.ixtb{1,dPARAMS.lab}(tbidx)); %start times of test bins in this session
    spCodeSet = dPARAMS.spCodeSet(dPARAMS.labelConfIdx); % labels in this session
    labConfSet = dPARAMS.labelConf(dPARAMS.labelConfIdx); % label confidences in this session
    
    % if data for any of these bins already exists in labelCertainty, remove it
    if ~isempty(labelCertainty{1,dPARAMS.lab})
        prevEntries = ismember(labelCertainty{1,dPARAMS.lab}(:,1),BT);
        labelCertainty{1,dPARAMS.lab}(prevEntries,:) = [];
    end
    
    for i = 1:length(BT) % for each test bin
        z = [];
        % identify indices, times, & labels of labeled clicks in this bin
        cidx = find(dPARAMS.tID >= BT(i) & dPARAMS.tID < (BT(i)+p.binDur/(24*60)));
        if ~isempty(cidx)
            % demarcate bin on TL vs. Time and ICI vs. Time plots
            hold(dHANDLES.LTSAsubs(1),'on')
            hold(dHANDLES.LTSAsubs(3),'on')
            if BT(i)<dPARAMS.sb(dPARAMS.k)
                xH201a = plot(dHANDLES.LTSAsubs(1),repmat(dPARAMS.sb(dPARAMS.k),81,1),100:1:180,...
                    'r','LineWidth',2);
                xH201c = plot(dHANDLES.LTSAsubs(3),repmat(dPARAMS.sb(dPARAMS.k),201,1),0:0.01:2,...
                    'r','LineWidth',2);
            else
                xH201a = plot(dHANDLES.LTSAsubs(1),repmat(BT(i),81,1),100:1:180,...
                    'r','LineWidth',2);
                xH201c = plot(dHANDLES.LTSAsubs(3),repmat(BT(i),201,1),0:0.01:2,...
                    'r','LineWidth',2);
            end
            if BT(i)+p.binDur/(24*60)>dPARAMS.eb(dPARAMS.k)
                xH201b = plot(dHANDLES.LTSAsubs(1),repmat(dPARAMS.eb(dPARAMS.k),...
                    81,1),100:1:180,'r','LineWidth',2);
                xH201d = plot(dHANDLES.LTSAsubs(3),repmat(dPARAMS.eb(dPARAMS.k),...
                    201,1),0:0.01:2,'r','LineWidth',2);
            else
                xH201b = plot(dHANDLES.LTSAsubs(1),repmat(BT(i)+p.binDur/(24*60),...
                    81,1),100:1:180,'r','LineWidth',2);
                xH201d = plot(dHANDLES.LTSAsubs(3),repmat(BT(i)+p.binDur/(24*60),...
                    201,1),0:0.01:2,'r','LineWidth',2);
            end
            hold(dHANDLES.LTSAsubs(1),'off')
            hold(dHANDLES.LTSAsubs(3),'off')
            
            tIDbin = dPARAMS.tID(cidx);
            spCodeSetbin = spCodeSet(cidx);
            labConfSetbin = labConfSet(cidx);
            
            % in case more than one spectrum in this bin received this label,
            % average confidences corresponding to this label
            conf = mean(labConfSetbin(spCodeSetbin==dPARAMS.lab));
            
            % find clicks with this label in session vars
            [~, kInd, ~] = intersect(dPARAMS.t,tIDbin(spCodeSetbin==dPARAMS.lab));
            maxRL = max(dPARAMS.RL(kInd));
            meanRL = mean(dPARAMS.RL(kInd));
            
            % average clicks with this label and plot them
            % calculate mean spectrum
            specs = dPARAMS.cspJ(kInd,:);
            mspec = mean(specs,1);
            minSpec = min(mspec,[],2);
            Spec = mspec - repmat(minSpec,1,size(mspec,2));
            mxSpec = max(Spec,[],2);
            Spec = Spec./repmat(mxSpec,1,size(mspec,2));
            meanSpec = 20*log10(nanmean(10.^Spec./20,1));
            meanMin = meanSpec-min(meanSpec);
            specNorm = meanMin./max(meanMin);
            % calculate ICI dist
            ICI = histcounts(diff(dPARAMS.t(kInd)*60*60*24),0:.01:1);
            % calculate mean waveform
            meanWav = norm_wav(mean(dPARAMS.csnJ(kInd,:),1));
            
            % plot
            figure(99);clf
            subplot(1,3,1)
            plot(dPARAMS.fmsp,specNorm);
            grid on
            xlim([dPARAMS.fmsp(1) dPARAMS.fmsp(end)]);
            xticks(10:10:90)
            xticklabels({'','20','','40','','60','','80'});
            xlabel('Frequency (kHz)');
            ylabel('Normalized Amplitude');
            subplot(1,3,2)
            histogram('BinEdges',0:.01:1,'BinCounts',ICI);
            grid on
            xticks([0 0.25 0.5 0.75 1]);
            xlabel('ICI (s)');
            ylabel('Counts');
            title(['Bin ',num2str(i),', Label ',num2str(dPARAMS.lab),': ',p.mySpID(dPARAMS.lab).Name]);
            subplot(1,3,3)
            plot(meanWav);
            xlabel('Sample Number');
            ylabel('Normalized Amplitude');
            
            while isempty(z)
                % ask for user input
                z = input('Are these clicks correctly labeled? Enter 1 for yes, 0 for no: ');
                if ~isempty(z) && z~=1 && z~=0
                    fprintf('WARNING: Entry not allowed\n');
                    z = [];
                elseif isempty(z)
                    fprintf('WARNING: Entry not allowed\n');
                end
            end
            labelCertainty{1,dPARAMS.lab} = [labelCertainty{1,dPARAMS.lab}(:,:);...
                [BT(i),conf,length(kInd),z]];
            RL{1,dPARAMS.lab} = [RL{1,dPARAMS.lab}(:,:); maxRL];
            RL{2,dPARAMS.lab} = [RL{2,dPARAMS.lab}(:,:); meanRL];
            
            % remove bin indicators from time series plots
            xH201a.XData = [];xH201a.YData = [];
            xH201b.XData = [];xH201b.YData = [];
            xH201c.XData = [];xH201c.YData = [];
            xH201d.XData = [];xH201d.YData = [];
        end
    end
    
    % remove redundancies, sort and save labelCertainty
    [~,uniqueID] = unique(labelCertainty{1,dPARAMS.lab}(:,1));
    labelCertainty{1,dPARAMS.lab} = labelCertainty{1,dPARAMS.lab}(uniqueID,:);
    [labelCertainty{1,dPARAMS.lab}, sortInd] = sortrows(labelCertainty{1,dPARAMS.lab});
    RL{1,dPARAMS.lab} = RL{1,dPARAMS.lab}(uniqueID);
    RL{1,dPARAMS.lab} = RL{1,dPARAMS.lab}(sortInd);
    RL{2,dPARAMS.lab} = RL{2,dPARAMS.lab}(uniqueID);
    RL{2,dPARAMS.lab} = RL{2,dPARAMS.lab}(sortInd);
    save(strrep(fNameList.TPWS,'TPWS1','labCert'),'labelCertainty','RL','p');
    
    % advance to next bout with test bins for this label (if there is one)
    if dPARAMS.k < dPARAMS.nb
        next_tbidx = find(dPARAMS.binTimes(dPARAMS.ixtb{1,dPARAMS.lab}) >= dPARAMS.eb(dPARAMS.k),1,'first');
        next_tbTime = dPARAMS.binTimes(dPARAMS.ixtb{1,dPARAMS.lab}(next_tbidx));
        if ~isempty(next_tbTime)
            next_sessionIdx = find(dPARAMS.sb-datenum([0,0,0,0,p.binDur,0])<=next_tbTime,1,'last');
            dPARAMS.k = next_sessionIdx;
        else
            fprintf(['No more test bins for label ',num2str(p.mySpID(dPARAMS.lab).zID_Label),' before end of file\n']);
            dPARAMS.k = 1;
        end
    elseif dPARAMS.k == dPARAMS.nb
        disp('End of file');
        dPARAMS.k = 1;
    end
    
else
    
    if ~isempty(labelCertainty{1,dPARAMS.lab})
        % remove redundancies, sort and save labelCertainty
        [~,uniqueID] = unique(labelCertainty{1,dPARAMS.lab}(:,1));
        labelCertainty{1,dPARAMS.lab} = labelCertainty{1,dPARAMS.lab}(uniqueID,:);
        [labelCertainty{1,dPARAMS.lab},sortInd] = sortrows(labelCertainty{1,dPARAMS.lab});
        RL{1,dPARAMS.lab} = RL{1,dPARAMS.lab}(uniqueID);
        RL{1,dPARAMS.lab} = RL{1,dPARAMS.lab}(sortInd);
        RL{2,dPARAMS.lab} = RL{2,dPARAMS.lab}(uniqueID);
        RL{2,dPARAMS.lab} = RL{2,dPARAMS.lab}(sortInd);
        save(strrep(fNameList.TPWS,'TPWS1','labCert'),'labelCertainty','RL','p');
    end
    
    % advance to next bout with test bins for this label (if there is one)
    if dPARAMS.k < dPARAMS.nb
        fprintf(['No test bins for label ',num2str(p.mySpID(dPARAMS.lab).zID_Label),' in this session, moving to next session with test bins for this label\n']);
        next_tbidx = find(dPARAMS.binTimes(dPARAMS.ixtb{1,dPARAMS.lab}) >= dPARAMS.eb(dPARAMS.k),1,'first');
        next_tbTime = dPARAMS.binTimes(dPARAMS.ixtb{1,dPARAMS.lab}(next_tbidx));
        if ~isempty(next_tbTime)
            next_sessionIdx = find(dPARAMS.sb-datenum([0,0,0,0,p.binDur,0])<=next_tbTime,1,'last');
            dPARAMS.k = next_sessionIdx;
        else
            fprintf(['No more test bins for label ',num2str(p.mySpID(dPARAMS.lab).zID_Label),' before end of file\n']);
            dPARAMS.k = 1;
        end
    elseif dPARAMS.k == dPARAMS.nb
        disp('End of file');
        dPARAMS.k = 1;
    end
    
end
dPARAMS.lab = [];
end