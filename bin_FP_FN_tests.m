function bin_FP_FN_tests
global dHANDLES dPARAMS p fpfnTD cMat

dPARAMS.lab = []; % ask which label to test
while isempty(dPARAMS.lab)
    dPARAMS.lab = input('Enter label to test: ');
    if ~isempty(dPARAMS.lab) & ~ismember(dPARAMS.lab,1:length(p.mySpID))
        disp('WARNING: Entry not allowed');
        dPARAMS.lab = [];
    end
end

% check for test bins in this session
btidx = dPARAMS.binTimes(dPARAMS.ftb{1,dPARAMS.lab}) >= dPARAMS.sb(dPARAMS.k) & ...
    dPARAMS.binTimes(dPARAMS.ftb{1,dPARAMS.lab}) < dPARAMS.eb(dPARAMS.k);
if any(btidx)
    
    %     for i = 2:size(zTD,2)
    %         zTD{dPARAMS.k,i}(3:6) = 0;
    %     end
    BT = dPARAMS.binTimes(dPARAMS.ftb{1,dPARAMS.lab}(btidx)); %start times of test bins
    
    for i = 1:length(BT) % for each test bin
        %FN = [];
        labs = [];
        z = [];
        a = [];
        b = [];
        specNorm = [];
        ICI = [];
        meanWav = [];
        conf = [];
        
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
        spCodeSet = dPARAMS.spCodeSet(dPARAMS.labelConfIdx); % labels in this session
        labConfSet = dPARAMS.labelConf(dPARAMS.labelConfIdx); % label confidences in this session
        cidx = find(tID >= BT(i) & tID < (BT(i)+p.binDur/(24*60)));
        tIDbin = tID(cidx);
        spCodeSetbin = spCodeSet(cidx);
        labConfSetbin = labConfSet(cidx);
        
        if ~isempty(cidx)
            % identify unique labels in this bin
            labs = unique(spCodeSetbin);
            
            %for j = 1:length(labs) % for each label in this bin  %%%%%%%%%%%%%%%
            if ismember(dPARAMS.lab,labs) % if selected label is in this bin
                
                % which clicks have this label?
                lidx = find(spCodeSetbin==dPARAMS.lab);
                % in case more than one spectrum in this bin received this label,
                % average confidences corresponding to this label
                conf = mean(labConfSetbin(lidx));
                
                %if length(lidx) > 25 %don't evaluate labels with very few clicks
                % add to running tally of FP test bins for this label
                fpfnTD{1,dPARAMS.lab} = [fpfnTD{1,dPARAMS.lab}(:,:);[BT(i),1,-1,conf,length(lidx)]];
                
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
                grid on
                xticks(10:10:90)
                xticklabels({'','20','','40','','60','','80'});
                xlim([dPARAMS.fmsp(1) dPARAMS.fmsp(end)]);
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
                    elseif ~isempty(z) && z==1 % record true positives for this label
                        fpfnTD{1,dPARAMS.lab}(end,3) = 1;
                        for k = 1:length(p.mySpID)
                            cMat{dPARAMS.lab,k}(1) = cMat{dPARAMS.lab,k}(1)+1;
                        end
                        cMat{dPARAMS.lab,dPARAMS.lab}(2) = cMat{dPARAMS.lab,dPARAMS.lab}(2)+1;
                        
                    elseif ~isempty(z) && z==0 % record false positives for this label
                        fpfnTD{1,dPARAMS.lab}(end,3) = 0;
                        
                        while isempty(a)
                            % check if also a false negative for another label
                            a = input('Enter label # these click should have been assigned to, or enter 0 if uncertain: ');
                            
                            if ~isempty(a) && ~ismember(a,1:size(zTD,2)-1) && a~=0
                                fprintf('WARNING: Entry not allowed\n');
                                a = [];
                                % count as FN only if this label doesn't already exist in this bin
                            elseif ~isempty(a) && a > 0 && ~ismember(a,labs)
                                fpfnTD{1,a} = [fpfnTD{1,a};BT(i),0,0,-1,-1];
                                %FN = [FN,a];
                                for k = 1:length(p.mySpID)
                                    cMat{a,k}(1) = cMat{a,k}(1)+1;
                                end
                                cMat{a,dPARAMS.lab}(2) = cMat{a,dPARAMS.lab}(2)+1;
                            end
                        end
                    end
                end
                %end
            else % if selected label is not in this bin, check for false negative
                for j = 1:size(labs) % calculate mean spec, ICI dist, & mean waveform for each label
                    specNorm = [];
                    ICI = [];
                    meanWav = [];
                    conf = [];
                    numClicks = [];
                    % which clicks have this label?
                    lidx = find(spCodeSetbin==labs(j));
                    % in case more than one spectrum in this bin received this label,
                    % average confidences corresponding to this label
                    conf = [conf;mean(labConfSetbin(lidx))];
                    numClicks = [numClicks;length(lidx)];
                    %                 % add to running tally of FP test bins for this label
                    %                 fpfnTD{1,dPARAMS.lab} = [fpfnTD{1,dPARAMS.lab};BT(i),1,-1,conf,length(lidx)];
                    
                    % find clicks with this label in session vars
                    [~, kInd, ~] = intersect(dPARAMS.t,tIDbin(lidx));
                    
                    % calculate click params
                    specs = dPARAMS.cspJ(kInd,:);
                    mspec = mean(specs,1);
                    minSpec = min(mspec,[],2);
                    Spec = mspec - repmat(minSpec,1,size(mspec,2));
                    mxSpec = max(Spec,[],2);
                    Spec = Spec./repmat(mxSpec,1,size(mspec,2));
                    meanSpec = 20*log10(nanmean(10.^Spec./20,1));
                    meanMin = meanSpec-min(meanSpec);
                    spNorm = meanMin./max(meanMin);
                    specNorm = [specNorm;spNorm];
                    % calculate ICI dist
                    ici = histcounts(diff(dPARAMS.t(kInd)*60*60*24),0:.01:1);
                    ICI = [ICI;ici];
                    % calculate mean waveform
                    wav = norm_wav(mean(dPARAMS.csnJ(kInd,:),1));
                    meanWav = [meanWav;wav];
                end
                % calculate mean spec, ICI dist, & mean waveform for any
                % unlabeled clicks in this bin
                m = find(dPARAMS.t(:,1)>= BT(i) & dPARAMS.t(:,1) < BT(i)+p.binDur/(24*60));
                [noLab, ia] = setdiff(dPARAMS.t(m),tIDbin);
                if ~isempty(noLab)
                    labs = [labs;0];
                    conf = [conf;NA];
                    numClicks = [numClicks;length(noLab)];
                    
                    % calculate unlabeled click params
                    specs = dPARAMS.cspJ(m(ia),:);
                    mspec = mean(specs,1);
                    minSpec = min(mspec,[],2);
                    Spec = mspec - repmat(minSpec,1,size(mspec,2));
                    mxSpec = max(Spec,[],2);
                    Spec = Spec./repmat(mxSpec,1,size(mspec,2));
                    meanSpec = 20*log10(nanmean(10.^Spec./20,1));
                    meanMin = meanSpec-min(meanSpec);
                    spNorm = meanMin./max(meanMin);
                    specNorm = [specNorm;spNorm];
                    % calculate ICI dist
                    ici = histcounts(diff(dPARAMS.t(m(ia))*60*60*24),0:.01:1);
                    ICI = [ICI;ici];
                    % calculate mean waveform
                    wav = norm_wav(mean(dPARAMS.csnJ(m(ia),:),1));
                    meanWav = [meanWav;wav];
                end
                % plot all labeled & unlabeled clicks in this bin
                figure(99);clf
                subplot(1,3,1)
                hold on
                for q = 1:size(specNorm,1)
                    plot(dPARAMS.fmsp,specNorm(q,:));
                end
                hold off
                xlim([dPARAMS.fmsp(1) dPARAMS.fmsp(end)]);
                xlabel('Frequency (kHz)');
                ylabel('Normalized Amplitude');
                legend(string(labs),'Location','southeast');
                subplot(1,3,2)
                hold on
                for q = 1:size(ICI,1)
                    bar(0:.01:1,ICI(q,:),'FaceAlpha',0.75);
                end
                hold off
                grid on
                xticks([0 0.25 0.5 0.75 1]);
                xlabel('ICI (s)');
                ylabel('Counts');
                title(['Bin ',num2str(i),', Labels ',num2str(labs),': ',p.mySpID(labs).Name]);
                subplot(1,3,3)
                hold on
                for q = 1:size(meanWav,1)
                    plot(meanWav(q,:) + ones(size(meanWav,2),1)'.*rand(1,min(size(meanWav,2))));
                end
                hold off
                xlabel('Sample Number');
                ylabel('Normalized Amplitude');
                while isempty(b)
                    % check for false negative(s)
                    b = input(['Enter label # which should have been classified as ',dPARAMS.lab,', or 99 if none: ']);
                    if ~isempty(b) && ~ismember(b,labs) && b~=99
                        fprintf('WARNING: Entry not allowed\n');
                        b = [];
                    elseif ~isempty(b) && ismember(b,labs) %&& ~ismember(b,FN)
                        x = find(labs==b);
                        fpfnTD{1,dPARAMS.lab} = [fpfnTD{1,dPARAMS.lab};BT(i),0,1,conf(x),numClicks(x)];
                        %FN = [FN,b];
                        for k = 1:length(p.mySpID)
                            cMat{dPARAMS.lab,k}(1) = cMat{dPARAMS.lab,k}(1)+1;
                        end
                        if b==0
                            cMat{dPARAMS.lab,end}(2) = cMat{dPARAMS.lab,end}(2)+1;
                        else
                            cMat{dPARAMS.lab,b}(2) = cMat{dPARAMS.lab,b}(2)+1;
                        end
                    end
                end
            end  %%%%%%%%%%%%%
        else
            % calculate mean spec, ICI dist, & mean waveform for any
            % unlabeled clicks in this bin
            m = find(dPARAMS.t(:,1)>= BT(i) & dPARAMS.t(:,1) < BT(i)+p.binDur/(24*60));
            [noLab, ia] = setdiff(dPARAMS.t(m),tIDbin);
            if ~isempty(noLab)
                labs = [labs;0];
                conf = [conf;NA];
                numClicks = [numClicks;length(noLab)];
                
                % calculate unlabeled click params
                specs = dPARAMS.cspJ(m(ia),:);
                mspec = mean(specs,1);
                minSpec = min(mspec,[],2);
                Spec = mspec - repmat(minSpec,1,size(mspec,2));
                mxSpec = max(Spec,[],2);
                Spec = Spec./repmat(mxSpec,1,size(mspec,2));
                meanSpec = 20*log10(nanmean(10.^Spec./20,1));
                meanMin = meanSpec-min(meanSpec);
                spNorm = meanMin./max(meanMin);
                specNorm = [specNorm;spNorm];
                % calculate ICI dist
                ici = histcounts(diff(dPARAMS.t(m(ia))*60*60*24),0:.01:1);
                ICI = [ICI;ici];
                % calculate mean waveform
                wav = norm_wav(mean(dPARAMS.csnJ(m(ia),:),1));
                meanWav = [meanWav;wav];
            end
            % plot all labeled & unlabeled clicks in this bin
            figure(99);clf
            subplot(1,3,1)
            hold on
            for q = 1:size(specNorm,1)
                plot(dPARAMS.fmsp,specNorm(q,:));
            end
            hold off
            xlim([dPARAMS.fmsp(1) dPARAMS.fmsp(end)]);
            xlabel('Frequency (kHz)');
            ylabel('Normalized Amplitude');
            legend(string(labs),'Location','southeast');
            subplot(1,3,2)
            hold on
            for q = 1:size(ICI,1)
                bar(0:.01:1,ICI(q,:),'FaceAlpha',0.75);
            end
            hold off
            grid on
            xticks([0 0.25 0.5 0.75 1]);
            xlabel('ICI (s)');
            ylabel('Counts');
            title(['Bin ',num2str(i),', Labels ',num2str(labs),': ',p.mySpID(labs).Name]);
            subplot(1,3,3)
            hold on
            for q = 1:size(meanWav,1)
                plot(meanWav(q,:) + ones(size(meanWav,2),1)'.*rand(1,min(size(meanWav,2))));
            end
            hold off
            xlabel('Sample Number');
            ylabel('Normalized Amplitude');
            while isempty(b)
                % check for false negative(s)
                b = input(['Enter label # which should have been classified as ',dPARAMS.lab,', or 99 if none: ']);
                if ~isempty(b) && ~ismember(b,labs) && b~=99
                    fprintf('WARNING: Entry not allowed\n');
                    b = [];
                elseif ~isempty(b) && ismember(b,labs) %&& ~ismember(b,FN)
                    x = find(labs==b);
                    fpfnTD{1,dPARAMS.lab} = [fpfnTD{1,dPARAMS.lab};BT(i),0,1,conf(x),numClicks(x)];
                    %FN = [FN,b];
                    for k = 1:length(p.mySpID)
                        cMat{dPARAMS.lab,k}(1) = cMat{dPARAMS.lab,k}(1)+1;
                    end
                    if b==0
                        cMat{dPARAMS.lab,end}(2) = cMat{dPARAMS.lab,end}(2)+1;
                    else
                        cMat{dPARAMS.lab,b}(2) = cMat{dPARAMS.lab,b}(2)+1;
                    end
                end
            end
        end
        % remove bin indicators from time series plots
        xH201a.XData = [];xH201a.YData = [];
        xH201b.XData = [];xH201b.YData = [];
        xH201c.XData = [];xH201c.YData = [];
        xH201d.XData = [];xH201d.YData = [];
        
    end
    % advance to next bout with test bins for this label (if there is one)
    if dPARAMS.k < dPARAMS.nb
        fprintf(['No test bins for label ',num2str(p.mySpID(dPARAMS.lab).zID_Label),' in this session, moving to next session with test bins for this label\n']);
        next_tbidx = find(dPARAMS.binTimes(dPARAMS.ftb{1,dPARAMS.lab}) >= dPARAMS.eb(dPARAMS.k),1,'first');
        next_tbTime = dPARAMS.binTimes(dPARAMS.ftb{1,dPARAMS.lab}(next_tbidx));
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
    % advance to next bout with test bins for this label (if there is one)
    if dPARAMS.k < dPARAMS.nb
        fprintf(['No test bins for label ',num2str(p.mySpID(dPARAMS.lab).zID_Label),' in this session, moving to next session with test bins for this label\n']);
        next_tbidx = find(dPARAMS.binTimes(dPARAMS.ftb{1,dPARAMS.lab}) >= dPARAMS.eb(dPARAMS.k),1,'first');
        next_tbTime = dPARAMS.binTimes(dPARAMS.ftb{1,dPARAMS.lab}(next_tbidx));
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
end