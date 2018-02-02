function meanSPEC = norm_spec_simple(spectraSet,fimint,fimaxt)
% normalize spectrum given cutoff frequencies

% focus on the freq range of interest
prunedSpectra = spectraSet(:,fimint:fimaxt); 
mnSPEC = min(prunedSpectra,[],2); % before 121 instead of flowt
SPEC = prunedSpectra - repmat(mnSPEC,1,size(prunedSpectra,2));  % make low freq part = 0
mxSPEC = max(SPEC,[],2);
SPEC = SPEC ./ repmat(mxSPEC,1,size(prunedSpectra,2));  % make high freq part = 1

meanSPEC = 20*log10(nanmean(10.^SPEC./20,1));
meanMin = meanSPEC-min(meanSPEC);
meanSPEC = meanMin./max(meanMin);
