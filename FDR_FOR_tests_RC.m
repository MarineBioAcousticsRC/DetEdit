function FDR_FOR_tests_RC
global dHANDLES dPARAMS p zTD cMat

% check for test bins in this session
btidx = dPARAMS.binTimes(dPARAMS.ixtb) >= dPARAMS.sb(dPARAMS.k) & ...
    dPARAMS.binTimes(dPARAMS.ixtb) < dPARAMS.eb(dPARAMS.k);
if any(btidx)
    
    for i = 2:size(zTD,2)
        zTD{dPARAMS.k,i}(3:6) = 0;
    end
    BT = dPARAMS.binTimes(dPARAMS.ixtb(btidx)); %start times of test bins
    
    for i = 1:length(BT) % for each test bin
        FN = [];
        labs = [];
        z = [];
        a = [];
        b = [];
        % demarcate bin on time series plots
        hold(dHANDLES.LTSAsubs(1),'on')
        xH201a = plot(dHANDLES.LTSAsubs(1),repmat(BT(i),61,1),110:1:170,...
            'r','LineWidth',2);
        xH201b = plot(dHANDLES.LTSAsubs(1),repmat(BT(i)+p.binDur/(24*60),...
            61,1),110:1:170,'r','LineWidth',2);
        hold(dHANDLES.LTSAsubs(1),'off')
        hold(dHANDLES.LTSAsubs(3),'on')
        xH201c = plot(dHANDLES.LTSAsubs(3),repmat(BT(i),101,1),0:0.01:1,...
            'r','LineWidth',2);
        xH201d = plot(dHANDLES.LTSAsubs(3),repmat(BT(i)+p.binDur/(24*60),...
            101,1),0:0.01:1,'r','LineWidth',2);
        hold(dHANDLES.LTSAsubs(3),'off')
        
        % identify indices of labeled clicks in this bin
        tID = dPARAMS.tID(dPARAMS.labelConfIdx);
        spCodeSet = dPARAMS.spCodeSet(dPARAMS.labelConfIdx);
        cidx = find(tID >= BT(i) & tID < (BT(i)+p.binDur/(24*60)));
        tIDbin = tID(cidx);
        spCodeSetbin = spCodeSet(cidx);
        
        if ~isempty(cidx)
            % identify unique labels in this bin
            labs = unique(spCodeSetbin);
            
            for j = 1:length(labs) % for each label in this bin
                z = [];
                a = [];
                % which clicks have this label?
                lidx = find(spCodeSetbin==labs(j));
                
                if length(lidx) > 25 %don't evaluate labels with very few clicks
                    % add to running tally of FDR test bins for this label
                    zTD{dPARAMS.k,labs(j)+1}(3) = zTD{dPARAMS.k,labs(j)+1}(3)+1;

                    % find clicks with this label in session vars
                    [~, kInd, ~] = intersect(dPARAMS.t,tIDbin(lidx));
                    
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
                    xlim([dPARAMS.fmsp(1) dPARAMS.fmsp(end)]);
                    xlabel('Frequency (kHz)');
                    ylabel('Normalized Amplitude');
                    subplot(1,3,2)
                    histogram('BinEdges',0:.01:1,'BinCounts',ICI);
                    grid on
                    xticks([0 0.25 0.5 0.75 1]);
                    xlabel('ICI (s)');
                    ylabel('Counts');
                    title(['Bin ',num2str(i),', Label ',num2str(labs(j)),': ',p.mySpID(labs(j)).Name]);
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
                        elseif ~isempty(z) && z==1 % record true positives for this label
                            for k = 1:length(p.mySpID)
                                cMat{labs(j),k}(1) = cMat{labs(j),k}(1)+1;
                            end
                            cMat{labs(j),labs(j)}(2) = cMat{labs(j),labs(j)}(2)+1;
                        elseif ~isempty(z) && z==0 % record false positives for this label
                            zTD{dPARAMS.k,labs(j)+1}(4) = zTD{dPARAMS.k,labs(j)+1}(4)+1;
                            
                            while isempty(a)
                                % check if also a false negative for another label
                                a = input('Enter label # these click should have been assigned to, or enter 0 if uncertain: ');
                                
                                if ~isempty(a) && ~ismember(a,1:size(zTD,2)-1) && a~=0
                                    fprintf('WARNING: Entry not allowed\n');
                                    a = [];
                                    % count as FN only if this label doesn't already exist in this bin
                                elseif ~isempty(a) && a > 0 && ~ismember(a,labs)
                                    zTD{dPARAMS.k,a+1}(5) = zTD{dPARAMS.k,a+1}(5)+1;
                                    zTD{dPARAMS.k,a+1}(6) = zTD{dPARAMS.k,a+1}(6)+1;
                                    FN = [FN,a];
                                    
                                    for k = 1:length(p.mySpID)
                                        cMat{a,k}(1) = cMat{a,k}(1)+1;
                                    end
                                    cMat{a,labs(j)}(2) = cMat{a,labs(j)}(2)+1;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % identify indices of unlabeled clicks in this bin
        m = find(dPARAMS.t(:,1)>= BT(i) & dPARAMS.t(:,1) < BT(i)+p.binDur/(24*60));
        [noLab, ia] = setdiff(dPARAMS.t(m),tIDbin);
        
        if ~isempty(noLab) && length(noLab) > 5
            % average unlabeled clicks
            % calculate mean spectrum
            specs = dPARAMS.cspJ(m(ia),:);
            mspec = mean(specs,1);
            minSpec = min(mspec,[],2);
            Spec = mspec - repmat(minSpec,1,size(mspec,2));
            mxSpec = max(Spec,[],2);
            Spec = Spec./repmat(mxSpec,1,size(mspec,2));
            meanSpec = 20*log10(nanmean(10.^Spec./20,1));
            meanMin = meanSpec-min(meanSpec);
            specNorm = meanMin./max(meanMin);            
            % calculate ICI dist
            ICI = histcounts(diff(dPARAMS.t(m(ia))*60*60*24),0:.01:1);
            % calculate mean waveform
            meanWav = norm_wav(mean(dPARAMS.csnJ(m(ia),:),1));
            
            % plot average of unlabeled clicks in this bin
            figure(99);clf
            subplot(1,3,1)
            plot(dPARAMS.fmsp,specNorm);
            xlim([dPARAMS.fmsp(1) dPARAMS.fmsp(end)]);
            xlabel('Frequency (kHz)');
            ylabel('Normalized Amplitude');
            subplot(1,3,2)
            histogram('BinEdges',0:.01:1,'BinCounts',ICI);
            grid on
            xticks([0 0.25 0.5 0.75 1]);
            xlabel('ICI (s)');
            ylabel('Counts');
            title(['Bin ',num2str(i),': Unlabeled Clicks']);
            subplot(1,3,3)
            plot(meanWav);
            xlabel('Sample Number');
            ylabel('Normalized Amplitude');
            
            while isempty(b)
                % check if false negative for any label
                b = input('Enter label # these click should have been assigned to, or enter 0 if uncertain: ');
                if ~isempty(b) && ~ismember(b,1:size(zTD,2)-1) && b~=0
                    fprintf('WARNING: Entry not allowed\n');
                    b = [];
                    % count as FN only if this label doesn't already exist in this bin
                elseif ~isempty(b) && b > 0 && ~ismember(b,labs) && ~ismember(b,FN)
                    zTD{dPARAMS.k,b+1}(5) = zTD{dPARAMS.k,b+1}(5)+1;
                    zTD{dPARAMS.k,b+1}(6) = zTD{dPARAMS.k,b+1}(6)+1;
                    FN = [FN,b];
                    
                    for k = 1:length(p.mySpID)
                        cMat{b,k}(1) = cMat{b,k}(1)+1;
                    end
                end
            end
        end
        
        if isempty(cidx) && isempty(noLab)
            disp('No clicks in this bin, moving on to next');
        end
        
        % record lack of false negatives for all other labels for this bin
        if ~isempty(labs)
            q = setdiff(1:size(zTD,2)-1,labs);
        else
            q = 1:size(zTD,2)-1;
        end
        r = setdiff(q,FN); % don't double count bin if FN has already been identified for some label(s)
        
        for l = 1:length(r)
            zTD{dPARAMS.k,r(l)+1}(5) = zTD{dPARAMS.k,r(l)+1}(5)+1;
        end
        
        % remove bin indicators from time series plots
        xH201a.XData = [];xH201a.YData = [];
        xH201b.XData = [];xH201b.YData = [];
        xH201c.XData = [];xH201c.YData = [];
        xH201d.XData = [];xH201d.YData = [];
        
    end
    % advance to next bout (if there is one)
    if dPARAMS.k < dPARAMS.nb
        dPARAMS.k = dPARAMS.k + 1;
    else
        disp('End of file');
    end
    
else
    for i = 2:size(zTD,2)
        zTD{dPARAMS.k,i}(3:6) = 0;
    end
    if dPARAMS.k < dPARAMS.nb
        disp('No test bins in this session, moving to next session');
        dPARAMS.k = dPARAMS.k + 1;
    else
        disp('End of file');
    end
end

end