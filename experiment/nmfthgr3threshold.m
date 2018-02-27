function [config, store, obs] = nmfthgr3threshold(config, setting, data)
% nmfth3threshold THRESHOLD step of the expLanes experiment NMFThreshold
%    [config, store, obs] = nmfth3threshold(config, setting, data)
%      - config : expLanes configuration state
%      - setting   : set of factors to be evaluated
%      - data   : processing data stored during the previous step
%      -- store  : processing data to be saved for the other steps
%      -- obs    : observations to be saved for analysis

% Copyright: <userName>
% Date: 01-Dec-2017

% Set behavior for debug mode
if nargin==0, NMFThresholdGrafic('do', 3, 'mask',{}); return; else store=[]; obs=[]; end

nmf =  data.nmf;
levels = cell(1,length(data.nmf));

if strcmp(setting.type,'nmf') && strcmp(setting.nmfType,'thresholded')
    distanceMethod = setting.distanceMethod;
    methodThreshold = setting.methodThreshold;
    
    valueThresholdCheck = 1;
    switch  methodThreshold
        case 'firm'
            threshold(1) = setting.thresholdFirmHigh;
            threshold(2) = setting.thresholdFirmLow;
            if threshold(1) < threshold(2)
                valueThresholdCheck = 0;
            end
        otherwise
            threshold = setting.threshold;
    end
    
    [soundMix] = cutSpectrogramEXP(nmf{1}.W0,setting);
    W0 = soundMix.W(1:soundMix.ind,:);
    [F,K] = size(W0);
    
    if valueThresholdCheck == 1
        for ii = 1:length(nmf)
            W = nmf{ii}.W(1:F,:);
            H = nmf{ii}.H;
            
            switch distanceMethod
                case 'cosine'
                    dist = sum(W.*W0)./(sqrt(sum(W.^2,1)).*sqrt(sum(W0.^2,1)));
                    dist(isnan(dist)) = 0;
                    [~, order] = sort(dist,'descend');
                    
                case 'none'
                    order = 1:K;
            end
            Wn = W(:,order);
            Hn = H(order,:);
            
            switch distanceMethod
                case 'none'
                    Wtraffic = Wn;
                    Htraffic = Hn;
                    
                otherwise
                    dist = dist(order);
                    vec = zeros(1,K);
                    
                    switch methodThreshold
                        case 'hard'
                            vec(dist>threshold) = 1;
                            vec(dist<=threshold) = 0;
                            
                        case 'firm'
                            vec(dist>threshold(1)) = 1;
                            vec(dist<=threshold(1) & dist>threshold(2)) = 2;
                            vec(dist<=threshold(2)) = 3;
                            
                            vec(vec==1) = 1;
                            if ~isempty(dist(vec==2))
                                a = dist(vec==2);
                                vec(vec==2) = (a-a(end))/(a(1)-a(end));
                            end
                            vec(vec==3) = 0;
                    end
                    Wtraffic = Wn.*repmat(vec,F,1);
                    Htraffic = Hn;
            end

            [LpTrafficEstimate,LeqTrafficEstimate] = estimationLpEXP(Wtraffic*Htraffic,setting);
            
            levels{ii}.LpTrafficEstimate = LpTrafficEstimate;
            levels{ii}.LeqTrafficEstimate = LeqTrafficEstimate;
            levels{ii}.LpTraffic = nmf{ii}.LpTraffic;
            levels{ii}.LeqTraffic = nmf{ii}.LeqTraffic;
            levels{ii}.LpGlobal = nmf{ii}.LpGlobal;
            levels{ii}.LeqGlobal = nmf{ii}.LeqGlobal;
            levels{ii}.cost = nmf{ii}.cost;
            levels{ii}.indice = data.nmf{ii}.indice;
        end
    else
        for ii = 1:length(nmf)
            levels{ii}.LpTrafficEstimate{1} = nan;
            levels{ii}.LeqTrafficEstimate(1) = nan;
            levels{ii}.LpTraffic = nmf{ii}.LpTraffic;
            levels{ii}.LeqTraffic = nmf{ii}.LeqTraffic;
            levels{ii}.LpGlobal = nmf{ii}.LpGlobal;
            levels{ii}.LeqGlobal = nmf{ii}.LeqGlobal;
            levels{ii}.cost = nmf{ii}.cost;
            levels{ii}.indice = data.nmf{ii}.indice;
        end
    end
    store.levels = levels;
else
    for ii = 1:length(nmf)
        levels{ii}.LeqTrafficEstimate = nmf{ii}.LeqTrafficEstimate;
        levels{ii}.LpTrafficEstimate = nmf{ii}.LpTrafficEstimate;
        levels{ii}.LpTraffic = nmf{ii}.LpTraffic;
        levels{ii}.LeqTraffic = nmf{ii}.LeqTraffic;
        levels{ii}.LpGlobal = nmf{ii}.LpGlobal;
        levels{ii}.LeqGlobal = nmf{ii}.LeqGlobal;
        levels{ii}.cost = nmf{ii}.cost;
        levels{ii}.indice = data.nmf{ii}.indice;
    end
    store.levels = levels;
end
