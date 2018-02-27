function NMF = nmfSemiSupervisedEXP(setting,sequentialData,soundMix,beta)
%% [W, H, cost] = nmfSemiSupervised(V,W,H,NMF,soundMix,ind)
% W_(A x (B+C)) = [M V];
% On entre le dictionnaire normal calcul� auquel on ajoute C �l�ments
% random

Wrand = setting.SS_sizeWrand;
W = soundMix.Wcut;
V = soundMix.V_scene;

binTemp = size(V,2);

if isempty(sequentialData)
    iteration = setting.iteration;
    rng(soundMix.seed); YInit =  rand(size(W,1),Wrand);
    rng(soundMix.seed); HInit = rand(size(W,2),binTemp);
    rng(soundMix.seed); ZInit = rand(Wrand,binTemp);
else
    % continuing step of the sequential run
    iteration = setting.iteration-sequentialData.numberIteration;
    HInit = sequentialData.H;
    ZInit = sequentialData.Z;
    YInit = sequentialData.Y;
end

NMF = algo_nmfSemiSupervised(HInit,W,YInit,ZInit,V,beta,iteration,soundMix,setting);
