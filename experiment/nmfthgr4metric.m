function [config, store, obs] = nmfthgr4metric(config, setting, data)                
% nmfth4metric METRIC step of the expLanes experiment NMFThreshold                 
%    [config, store, obs] = nmfth4metric(config, setting, data)                    
%      - config : expLanes configuration state                                     
%      - setting   : set of factors to be evaluated                                
%      - data   : processing data stored during the previous step                  
%      -- store  : processing data to be saved for the other steps                 
%      -- obs    : obsrmservations to be saved for analysis                           
                                                                                   
% Copyright: <userName>                                                            
% Date: 01-Dec-2017                                                                
                                                                                   
% Set behavior for debug mode                                                      
if nargin==0, NMFThresholdGrafic('do', 4, 'mask', {2 0 0 0 1, ...
	3 0 6 0 3, 0 3 11 1 1, 1 1 10 1 3, 1 0 42:47 20:40}); return; else store=[]; obs=[]; end                                                                                   

levels = data.levels;
pref = setting.p0;

LeqTrafficEstimate = zeros(length(levels),1);
LeqTrafficEstimatedB = zeros(length(levels),1);
LeqTrafficExact = zeros(length(levels),1);
LeqTrafficExactdB = zeros(length(levels),1);
LeqGlobal = zeros(length(levels),1);
rmsedB = zeros(length(levels),1);
cost = zeros(length(levels),1);
rmse = zeros(length(levels),1);
nrmse = zeros(length(levels),1);

%% EXTRACTION LEVELS
for ii = 1:length(levels)
    LpTrafficExact = levels{ii}.LpTraffic;
    LpTrafficExact(LpTrafficExact<2e-5) = pref;
    
    LpTrafficEstimate = levels{ii}.LpTrafficEstimate{1};
    LpTrafficEstimate(LpTrafficEstimate==0) = pref;
    
    LpGlobal = levels{ii}.LpGlobal;
    LpGlobal(LpGlobal==0) = pref;
       
    LeqTrafficEstimate(ii) = levels{ii}.LeqTrafficEstimate(1);
    LeqTrafficExact(ii) = levels{ii}.LeqTraffic;
    LeqGlobal(ii) = levels{ii}.LeqGlobal;
    
    if size(LpGlobal,2)>size(LpTrafficEstimate,2)
        LpGlobal = LpGlobal(:,1:size(LpTrafficEstimate,2));
    elseif size(LpGlobal,2)<size(LpTrafficEstimate,2)
        LpTrafficEstimate = LpTrafficEstimate(:,1:size(LpGlobal,2));
    end
    
    
    %% linear sound level
    rmse(ii) = rmseEXP(LpTrafficExact,LpTrafficEstimate);
    nrmse(ii) = rmseEXP(LpTrafficExact./LpGlobal,LpTrafficEstimate./LpGlobal);
    
    %% sound level in dB
    LpTrafficExactdB{ii} = 20*log10(LpTrafficExact/pref);
    LpTrafficEstimatedB{ii} = 20*log10(LpTrafficEstimate/pref);
    LeqTrafficEstimatedB(ii) = 20*log10(LeqTrafficEstimate(ii)/pref);
    LeqTrafficExactdB(ii) = 20*log10(LeqTrafficExact(ii)/pref);
    
    rmse_temp = rmseEXP(LpTrafficExactdB{ii},LpTrafficEstimatedB{ii});
    rmsedB(ii) = mean(rmse_temp);
    cost(ii) = levels{ii}.cost(end);
    
   if isfield(levels{ii},'indice')
        obs.LAeq(ii) = levels{ii}.indice.LAeq;
        obs.Leq(ii) = levels{ii}.indice.Leq;
        obs.L10(ii) = levels{ii}.indice.L10;
        obs.L50(ii) = levels{ii}.indice.L50;
        obs.L90(ii) = levels{ii}.indice.L90;
        obs.L05(ii) = levels{ii}.indice.L05;
        obs.L1k(ii) = levels{ii}.indice.L1k;
        obs.L2k(ii) = levels{ii}.indice.L2k;
        obs.L5k(ii) = levels{ii}.indice.L5k;
   else
       pause
    end
end
mae = sum(abs(LeqTrafficEstimatedB-LeqTrafficExactdB))/length(LeqTrafficEstimatedB);
rmseLeq = rmseEXP(LeqTrafficEstimate,LeqTrafficExact);
rmseLeqdB = rmseEXP(LeqTrafficEstimatedB,LeqTrafficExactdB);

obs.LeqGlobal = 20*log10(LeqGlobal./pref);
obs.LeqTrafficExa = LeqTrafficExactdB;
obs.LeqTrafficEst = LeqTrafficEstimatedB;
obs.rmsedB = rmsedB;
obs.cost = cost;
obs.rmseLeq = rmseLeq;
obs.rmseLeqdB = rmseLeqdB;
obs.mae = mae;
obs.rmse = rmse;
obs.nrmse = nrmse;

