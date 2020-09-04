function FDR_FOR_tests_RC
global dHANDLES dPARAMS p zTD zID

% dPARAMS.lab = input('Enter label to test: ');
% thisLabField = sprintf('Label_%d', dPARAMS.lab);

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
        FNunID = [];
        
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
        cidx = find(zID(:,1)>= BT(i) & zID(:,1) < BT(i)+p.binDur/(24*60));
        
        if ~isempty(cidx)
            % identify unique labels in this bin
            labs = unique(zID(cidx,2));
            
            for j = 1:length(labs) % for each label in this bin
                
                % add to running tally of FDR test bins for this label
                zTD{dPARAMS.k,labs(j)+1}(3) = zTD{dPARAMS.k,labs(j)+1}(3)+1;
                
                % which clicks have this label?
                lidx = find(zID(cidx,2)==labs(j));
                
                % find clicks with this label in session vars
                [~, kInd, ~] = intersect(dPARAMS.t,zID(cidx(lidx),1));
                
                % average clicks with this label and plot them
                % calculate mean spectrum
                specs = dPARAMS.cspJ(kInd,:);
                specMax = max(specs,[],2);
                specNorm = mean(specs./specMax,1);
                % calculate mean ICI dist
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
                title(['Label ',num2str(labs(j)),': ',p.mySpID(labs(j)).Name]);
                subplot(1,3,3)
                plot(meanWav);
                xlabel('Sample Number');
                ylabel('Normalized Amplitude');
                
                % ask for user input
                z = input('Are these clicks correctly labeled? Enter 1 for yes, 0 for no: ');
                
                if z==0 % record false positives for this label
                    zTD{dPARAMS.k,labs(j)+1}(4) = zTD{dPARAMS.k,labs(j)+1}(4)+1;
                    
                    % check if also a false negative for another label
                    a = input('Enter label # these click should have been assigned to, or enter 0 if uncertain: ')
                    if a > 0
                        zTD{dPARAMS.k,a+1}(5) = zTD{dPARAMS.k,a+1}(5)+1;
                        zTD{dPARAMS.k,a+1}(6) = zTD{dPARAMS.k,a+1}(6)+1;
                        FN = [FN,a];
                    end
                end
            end
        end
        
        % identify indices of unlabeled clicks in this bin
        m = find(dPARAMS.t(:,1)>= BT(i) & dPARAMS.t(:,1) < BT(i)+p.binDur/(24*60));
        [noLab, ia] = setdiff(dPARAMS.t(m),zID(cidx,1));
        
        if ~isempty(noLab)
            % average unlabeled clicks
            % calculate mean spectrum
            specs = dPARAMS.cspJ(m(ia),:);
            specMax = max(specs,[],2);
            specNorm = mean(specs./specMax,1);
            % calculate mean ICI dist
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
            set(gca,'fontSize',14);
            title('Unlabeled Clicks');
            subplot(1,3,3)
            plot(meanWav);
            xlabel('Sample Number');
            ylabel('Normalized Amplitude');
            
            % check if false negative for any label
            b = input('Enter label # these click should have been assigned to, or enter 0 if uncertain: ')
            if b > 0
                zTD{dPARAMS.k,b+1}(5) = zTD{dPARAMS.k,b+1}(5)+1;
                zTD{dPARAMS.k,b+1}(6) = zTD{dPARAMS.k,b+1}(6)+1;
                FN = [FN,b];
            end
        end
        
        if isempty(cidx) && isempty(noLab)
            disp('No clicks in this bin, moving on to next');
        end
        
        % record lack of false negatives for all other labels for this bin
        if ~isempty(labs)
            q = setdiff(1:size(zTD,2)-1,labs);
        else
            q = 1:size(zTD,2)-1
        end
        r = setdiff(q,FN); % don't double count bin if FN has already been identified
        
        for l = 1:length(r)
            zTD{dPARAMS.k,r(l)+1}(5) = zTD{dPARAMS.k,r(l)+1}(5)+1;
        end
        
        % remove bin indicators from time series plots
        xH201a.XData = [];xH201a.YData = [];
        xH201b.XData = [];xH201b.YData = [];
        xH201c.XData = [];xH201c.YData = [];
        xH201d.XData = [];xH201d.YData = [];
        
    end
    
else
    disp('No test bins in this session, moving to next session\n');
    dPARAMS.k = dPARAMS.k + 1;
end

end