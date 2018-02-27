function [NMF] = nmfSupervise_MMEXP(setting,sequentialData,soundMix,beta)
%% function [NMFOUT] = beta_nmfSupervise_MM(V, W, HInit, beta, NMF, setting, sample, ind)

W = soundMix.Wcut;
V = soundMix.V_scene;

binTemp = size(V,2);

if isempty(sequentialData)
    iteration = setting.iteration;
    rng(soundMix.seed)
    HInit = rand(size(W,2),binTemp);
else
    % continuing step of the sequential run
    iteration = setting.iteration-sequentialData.numberIteration;
    HInit = sequentialData.H;
end

NMF = algo_nmfSupervised(HInit,W,V,beta,iteration,soundMix,setting);
