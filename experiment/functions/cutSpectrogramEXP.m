function [soundMix] = cutSpectrogramEXP(W,setting)

switch setting.domain
    case 'mel'
        sr = setting.sr;
        numberMel = setting.numberMel; 

        [W,~] = spectre2MelEXP(W,numberMel,setting.cutOffFreq,sr);
        ind = numberMel;
        
    case 'thirdOctave'
        frequency = linspace(0,setting.sr/2,size(W,1));
        [W,~,~,~] = NarrowToNthOctaveEXP(frequency,W,3);
        ind = size(W,1);
        
    otherwise
        frequency = linspace(0,setting.sr/2,size(W,1));
        [~,ind] = min(abs(setting.cutOffFreq-frequency));
end

soundMix.ind = ind;
soundMix.seed = rng;
soundMix.W = W;
