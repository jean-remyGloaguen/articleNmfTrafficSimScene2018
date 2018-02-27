function [LeqFiltre,LpFiltre] = filtrePasseBasEXP(X,setting)
%% Y = filtrePasseBas(X,soundMix,varargin)
% Fonction filtre pour un signal temporel x filtrer par un filtre de
% Butterworth. Les paramètres se trouvent dans le structure soundMix :
% - la fréquence de coupure, .cutOffFreq,
% - l'ordre du filtre, .order.
% Si ils ne sont pas présent : par défaut cutOffFreq = 1000 et order = 2.
% Si on renseigne un troisième argument d'une valeur quelconque, on trace
% la fonction filtre par la fonction fvtool.

sr = setting.sr;
nfft = setting.nfft;
fc = setting.cutOffFreq;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FONCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frequency = linspace(0,sr/2,nfft/2+1);
[~,ind] = min(abs(fc-frequency));
X{1} = X{1}(1:ind,:);

[LpFiltre,LeqFiltre] = estimationLpEXP(X,setting);

