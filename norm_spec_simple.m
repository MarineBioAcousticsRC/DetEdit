function meanSPEC = norm_spec_simple(spectraSet, flowt,fimint,fimaxt)
% normalize spectrum given cutoff frequencies
mnSPEC = min(spectraSet(:,1:121),[],2);
SPEC = spectraSet - repmat(mnSPEC,1,size(spectraSet,2));  % make low freq part = 0
mxSPEC = max(SPEC(:,1:121),[],2);
SPEC = SPEC ./ repmat(mxSPEC,1,size(spectraSet,2));  % make high freq part = 1

meanSPEC = 20*log10(nanmean(10.^SPEC./20,1));
meanMin = meanSPEC-min(meanSPEC(1:121));
meanSPEC = meanMin./max(meanMin(1:121));
