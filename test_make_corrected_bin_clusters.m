% test make corrected bin clusters
%
%
IDcorrectedFname = "J:\GOM_Y4A_01\GOM_Y4A_01_TPWS\GOM_Y4A_01_disk01a_Delphin_ID1.mat";
IDorigFname = "J:\GOM_Y4A_01\GOM_Y4A_01_binLabels_ID\GOM_Y4A_01_disk01a_Delphin_ID1.mat";
TPWSName = "J:\GOM_Y4A_01\GOM_Y4A_01_TPWS\GOM_Y4A_01_disk01a_Delphin_TPWS1.mat";
netOutDir = "J:\updatedNetsTest";

IDcorrected = load(IDcorrectedFname);
IDorig = load(IDorigFname);
zIDcorrected = IDcorrected.zID;
zIDorig = IDorig.zID;

[newSpec, newWaveEnv, newICI, newID] = ...
    make_corrected_bin_clusters(zIDcorrected, zIDorig,TPWSName);

% Now load training set, add these examples to set
load("D:\Shared drives\MBARC_GOM\GOM x SIO x NOAA\Code\Feb2023\GOM_deepSites_Jan2023_bin_train.mat")
% load network, run some iterations
load("D:\Shared drives\MBARC_GOM\GOM x SIO x NOAA\Code\Feb2023\TrainedNet\GOM_deepSites_Jan2023_trainedNetwork_bin.mat")
newExamples_nostd = [newSpec, newICI, newWaveEnv];
% ICIstart = trainTestSetInfo.setSpecHDim+1;
% trainTestSetInfo.iciStd = mean(max(trainDataAll(:,ICIstart:(ICIstart+trainTestSetInfo.setICIHDim-1)),[],2));
%
% normSpec1 = newSpec-trainTestSetInfo.specStd(1);
% normSpec = normSpec1./(trainTestSetInfo.specStd(2)-trainTestSetInfo.specStd(1));
%
% normICI = newICI./max(newICI,[],2);
% normICI(isnan(normICI))=0;
% normWave = newWaveEnv/trainTestSetInfo.waveStd;
% newTrainingData = [normSpec,normICI,normWave];
newTrainingData = nn_fn_standardize_data(trainTestSetInfo,newExamples_nostd)
oldNet = net;
% [~, trainPrefs] = nn_build_network(trainTestSetInfo);
dlnet = dag2dlnetwork(oldNet);
trainPrefs.BatchNormalizationStatistics = "auto";
trainPrefs.MiniBatchSize = 25;
trainPrefs.InitialLearnRate = 0.0003;
trainPrefs.MaxEpochs = 5;
trainPrefs.ValidationData = [];
trainPrefs.ValidationPatience = inf;
getRandInt = size(trainDataAll,1)./length(typeNames)/10;

trainDataAllstd = nn_fn_standardize_data(trainTestSetInfo,trainDataAll);
randSet = round(rand(1)*getRandInt):getRandInt:(size(trainDataAllstd,1));
newTrainingDataAugmented = [newTrainingData;trainDataAllstd(randSet,:)];
newIDAugmented = [newID;trainLabelsAll(randSet,:)];

% sanity check plots

figure(30)
subplot(1,2,1)
imagesc(trainDataAllstd')
set(gca,'ydir','normal')
colorbar

subplot(1,2,2)
imagesc(newTrainingData')
set(gca,'ydir','normal')
colorbar
mtit('Do these look similar? If not, stop here')
%%
newTrainingDataR = reshape(newTrainingDataAugmented,[1,size(newTrainingDataAugmented,2),1,...
    size(newTrainingDataAugmented,1)]);
newnet = trainnet(newTrainingDataR,categorical(newIDAugmented),dlnet,"crossentropy",trainPrefs);
netSaveFile = fullfile(netOutDir,sprintf("updatedNet_%s.mat",datestr(now,'yymmddHHMMSS')));
net = newnet;
save(netSaveFile,...
    'net','netTrainingInfo','trainTestSetInfo','typeNames','trainPrefs',...
    'newTrainingDataR','newIDAugmented')

% reclassify bin files
global REMORA
REMORA.nn.classify.searchSubDirsTF = 0;
REMORA.nn.classify.wildcard = 'GOM_Y4A_01_disk01a';
%dirTemp = dirList{iDir};
REMORA.nn.classify.inDir = 'J:\GOM_Y4A_01\GOM_Y4A_01_ClusterBins_125_noNorm'; %fullfile(dirTemp.folder,dirTemp.name);
REMORA.nn.classify.saveDir = 'J:\GOM_Y4A_01\New_label_test';
excludeLabels = [];
if ~isfolder(REMORA.nn.classify.saveDir)
    mkdir(REMORA.nn.classify.saveDir)
end
REMORA.nn.classify.networkPath= netSaveFile;
fprintf('\n DET: Beginning classification of %s\n', REMORA.nn.classify.inDir)
nn_fn_classify_bins(excludeLabels)
fprintf('DETS: Successfully classified %s\n',REMORA.nn.classify.inDir)

REMORA.nn.exportzID.inDir = REMORA.nn.classify.saveDir;
[~,myFolderStr,~] = fileparts(REMORA.nn.classify.saveDir);
REMORA.nn.exportzID.wildcard = [];
%REMORA.nn.exportzID.saveDir = fullfile(REMORA.nn.classify.saveDir,strrep(myFolderStr,'_noNorm','_ID_noNorm'));
REMORA.nn.exportzID.saveDir = [REMORA.nn.classify.saveDir,'_ID'];
nn_fn_exportzID

% regenerate labels for data coming after.
