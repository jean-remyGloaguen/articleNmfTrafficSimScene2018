function [config, store, obs] = nmfthgr1dictionary(config, setting, data)            
% nmfth1dictionary DICTIONARY step of the expLanes experiment NMFThreshold         
%    [config, store, obs] = nmfth1dictionary(config, setting, data)                
%      - config : expLanes configuration state                                     
%      - setting   : set of factors to be evaluated                                
%      - data   : processing data stored during the previous step                  
%      -- store  : processing data to be saved for the other steps                 
%      -- obs    : observations to be saved for analysis                           
                                                                                   
% Copyright: <userName>                                                            
% Date: 01-Dec-2017                                                                
                                                                                   
% Set behavior for debug mode                                                      
if nargin==0, NMFThresholdGrafic('do', 1, 'mask',{0 0 0 0 1}); return; else store=[]; obs=[]; end
                                                                                   
if strcmp(setting.type,'nmf')
    K = setting.numberElement;
    inputPath = config.inputPath;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% NAME & NUMBER SAMPLE EXTRACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cityCarSamples = dir(strcat(inputPath,'dictionary',filesep,'*cityCar*'));
    roadCarSamples = dir(strcat(inputPath,'dictionary',filesep,'*roadCar*'));
    stopCarSamples = dir(strcat(inputPath,'dictionary',filesep,'*stopCar*'));
    sample = [cityCarSamples; roadCarSamples; stopCarSamples];
    indiceClass = ones(1,length(sample));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% AUDIO SAMPLE EXTRACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    W = cell(1,length(sample));
    freqAxis = cell(1,length(sample));
    indTraffic = W;
    
    parfor ind = 1:length(sample)
        file = audioread(strcat(inputPath,'dictionary',filesep,sample(ind).name));
        
        [W{ind},freqAxis{ind},indTraffic{ind}] = creationDictionaryEXP(file,indiceClass(ind),setting);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLUSTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    W = cell2mat(W);
    indTraffic = cell2mat(indTraffic);
    W(W<eps) = 0;
    
    if size(W,2) > K
        idx = zeros(K,1);
        switch setting.reduceSizeW
            case 'rand'
                idx = randperm(size(W,2),K);
                W = W(:,idx);
                indTraffic = indTraffic(idx);
                
            case 'kmeans'
                [~,C] = kmeans(W',K,'MaxIter',200);
                C = (C./repmat(sum(C,2),1,size(C,2)))';
                for ii = 1:K
                    [~,idx(ii)] = min(0.5*sum((repmat(C(:,ii),1,size(W,2))-W).^2,1));
                end
                indTraffic = indTraffic(idx);
                W = C;
                
        end
    end
    
    dictionary.className = setting.dictionaryClass;
    dictionary.indTraffic = indTraffic;
    dictionary.W = W;
    dictionary.frequency = freqAxis{1};
    dictionary = orderfields(dictionary);
else
    
    dictionary = [];
end

store.dictionary = dictionary;