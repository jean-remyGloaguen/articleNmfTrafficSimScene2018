function [NMF] = separationTrafficEXP(NMF,soundMix,setting)
%% [NMF] = separationActivateur(NMF,soundMix,ind)
% séparation des activateurs pour chaque classe du dictionnaire
%

Htraffic = NMF.H(soundMix.indTraffic==1,:);
Wtraffic = soundMix.W(:,soundMix.indTraffic==1);

Vtraffic_full = Wtraffic*Htraffic;

if strcmp(setting.domain,'mel')
    
    Vtraffic_full = mel2SpectreEXP(Vtraffic_full,setting.sr,setting.nfft);
    [Lp,Leq] = estimationLpEXP(Vtraffic_full,setting);
    
elseif strcmp(setting.domain,'thirdOctave')
    
    [Lp,Leq] = estimationLpEXP(Vtraffic_full,setting);
else
    
    Vtraffic_cut = Wtraffic(1:soundMix.ind,:)*Htraffic;
    [Lp,Leq] = estimationLpEXP([{Vtraffic_full} {Vtraffic_cut}],setting);
end

NMF.LeqTrafficEstimate = Leq;
NMF.LpTrafficEstimate = Lp;
