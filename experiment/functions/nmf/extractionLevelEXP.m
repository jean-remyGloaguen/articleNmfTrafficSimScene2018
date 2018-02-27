function NMF = extractionLevelEXP(NMF,soundMix,V,beta,setting,varargin)

sparsity = setting.sparsity; 
smoothness = setting.smoothness; 
nmfType = setting.nmfType; 

Vap = NMF.Vap;
H = NMF.H;

if strcmp(setting.domain,'mel')
    Y = mel2SpectreEXP([{V} {Vap}],setting.sr,setting.nfft,setting.cutOffFreq);
    V = Y{1};
    Vap = Y{2};
end

switch nmfType
    case 'supervised'
        costTemp = betadivEXP(V,Vap,beta,H,sparsity,smoothness);
        dictionary = soundMix.W;
        activator = H;
        
    case 'semi-supervised'
        Y = varargin{1};
        Z = varargin{2};
        size_Wrand = setting.SS_sizeWrand; 
        constraint = setting.SS_constraint; 

        YFull = [Y; zeros(size(soundMix.Wfull,1)-size(Y,1),size_Wrand)];
        dictionary = [soundMix.W YFull];
        activator = [H; Z];

        switch constraint
            case 0
                costTemp = betadivEXP(V,Vap,beta,H,sparsity,smoothness);

            otherwise
                betaM = setting.SS_betaM;
                lambda = setting.SS_lambda;

                div = betadivEXP(repmat(soundMix.Wcut,size(Y,2),1),Y(:)*ones(1,size(soundMix.Wcut,2)),betaM);
                Cp = exp((-1/lambda)*sum(div));
                costTemp = betadivEXP(V,Vap,beta,H,sparsity,smoothness)+Cp;
        end
        NMF.Y = Y;
        NMF.Z = Z;
end

Vap_full = dictionary*activator;

[Lp,Leq] = estimationLpEXP([{Vap} {Vap_full}],setting);

NMF.H = H;
NMF.cost = costTemp;
NMF.LeqGlobalEstimate = Leq;
NMF.LpGlobalEstimate = Lp;

NMF = separationTrafficEXP(NMF,soundMix,setting);
NMF = orderfields(NMF);