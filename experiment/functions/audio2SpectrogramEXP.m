function [X,Xlinear,time,frequency,XScale] = audio2SpectrogramEXP(x,setting)
%% [X,Xlinear,time,XScale,frequency] = audio2Spectrogram(x,setting)

sr = setting.sr;
nfft = setting.nfft;
window = setting.nfft;
noverlap = floor(window*setting.noverlap/100);

if size(x,1)>size(x,2)
    x = x';
end

X = cell(1,size(x,1));
XScale = X;
frequency = X;
time = X;
Xlinear = X;

for ii = 1:size(x,1)
    [spec,frequency{ii},time{ii}] = spectrogram(x(ii,:),window,noverlap,nfft,sr);
    spec = spec./nfft;
    Xlinear{ii} = abs(spec);
    
    if strcmp('thirdOctave',setting.domain)
        [X{ii},XScale{ii},~,~] = NarrowToNthOctaveEXP(frequency{ii},Xlinear{ii},3);
    elseif strcmp('mel',setting.domain)
        [X{ii},XScale{ii}] = spectre2MelEXP(Xlinear{ii},setting.numberMel,setting.cutOffFreq,setting.sr);
    else
        X{ii} = Xlinear{ii};
        XScale{ii} = frequency{ii};
    end
end
