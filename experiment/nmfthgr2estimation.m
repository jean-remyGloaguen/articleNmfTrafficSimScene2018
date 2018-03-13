function [config, store, obs] = nmfthgr2estimation(config, setting, data)            
% nmfth2estimation ESTIMATION step of the expLanes experiment NMFThreshold         
%    [config, store, obs] = nmfth2estimation(config, setting, data)                
%      - config : expLanes configuration state                                     
%      - setting   : set of factors to be evaluated                                
%      - data   : processing data stored during the previous step                  
%      -- store  : processing data to be saved for the other steps                 
%      -- obs    : observations to be saved for analysis                           
                                                                                   
% Copyright: <userName>                                                            
% Date: 01-Dec-2017                                                                
                                                                                   
% Set behavior for debug mode                                                      
if nargin==0, NMFThresholdGrafic('do', 2, 'mask', {}); return; else store=[]; obs=[]; end

dataset = setting.dataset;
type = setting.type;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DICTIONARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(type,'nmf')
    dictionary.W = data.dictionary.W;
    dictionary.frequency = data.dictionary.frequency;
    dictionary.indTraffic = data.dictionary.indTraffic;
else
    dictionary = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BASELINE GLOBALE ERROR + FILTRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

creationSceneDir = strcat(config.inputPath, dataset, filesep, setting.gType, filesep);

files = dir(strcat(creationSceneDir,filesep,'*.wav'));
globalName = cell(1,length(files)/2);
ind = 1;

for ii = 1:length(files)
    if strfind(files(ii).name,'traffic')
        globalName{ind} = files(ii).name(1:end-12);
        ind = ind+1;
    end
end

numberScene = length(globalName);
nmf = cell(numberScene,1);

if strcmp(type,'nmf')
    nmf{1}.W0 = dictionary.W;
end

if isempty(config.sequentialData)
    sequentialData = cell(numberScene,1);
else
    sequentialData = config.sequentialData;
end

parfor ii = 1:numberScene
    
    fileTraffic = audioread(strcat(creationSceneDir,globalName{ii},'_traffic.wav'));
    fileRest = audioread(strcat(creationSceneDir,globalName{ii},'_interfering.wav'));
    
    fileTraffic(fileTraffic == 0) = eps;
    fileTot = fileRest+fileTraffic;
    
    [Vtraffic] = audio2SpectrogramEXP(fileTraffic',setting);
    [Lp,Leq] = estimationLpEXP(Vtraffic,setting);
    nmf{ii}.LpTraffic =  Lp{1};
    nmf{ii}.LeqTraffic = Leq(1);
    
    [V,Vlinear] = audio2SpectrogramEXP(fileTot',setting);
    [Lp,Leq] = estimationLpEXP(Vlinear,setting);
    nmf{ii}.LpGlobal = Lp{1};
    nmf{ii}.LeqGlobal = Leq(1);
    
    indice = indiceEstimationEXP(fileTot);
    
    switch type
        case 'filter'
            [LeqFiltre,LpFiltre] = filtrePasseBasEXP(Vlinear,setting);
            nmf{ii}.LeqTrafficEstimate = LeqFiltre;
            nmf{ii}.LpTrafficEstimate = LpFiltre;
            sequentialData{ii} = 0;
            nmf{ii}.cost = 0;
            nmf{ii}.indice = indice;
            
        case 'nmf'
            [NMF,sequentialData{ii}] =...
                NMFestimationEXP(V,dictionary,setting,sequentialData{ii});
            
            indice = indiceEstimationEXP(fileTot);
            
            if strcmp(setting.nmfType,'thresholded')
                nmf{ii}.H = NMF.H;
                nmf{ii}.W = NMF.W;
            end
            
            nmf{ii}.cost = NMF.cost(end);
            nmf{ii}.indice = indice;
    end 
end

config.sequentialData = sequentialData;
store.nmf = nmf;
