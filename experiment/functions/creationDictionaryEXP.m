function [W,Xscale,indTraffic] = creationDictionaryEXP(file,indClass,setting)


if min(size(file)) == 2        % stereo to mono
    file = mean(file,2)';
else
    file = file';
end

spec = audio2SpectrogramEXP(file,setting);

switch setting.frameSize
    case 0 % W composed of fft of sounds get from the rms of the spectrogram
        X = rms(spec{1},2);
        
    otherwise   % W composed of several temporal window get by the rms of the selected window
        step = (setting.window*(setting.noverlap/100))/setting.sr;
        time = step:step:size(spec{1},2)*step;
        frameSize_temp = setting.frameSize*1e-3;    % size of the temporal window
        timeEnd = floor(time(end));
        if time(end) < 1
           timeEnd = frameSize_temp; 
        end
        [~,idxBeg] = min(abs(repmat((0:frameSize_temp:timeEnd-frameSize_temp)',1,length(time))-repmat(time,timeEnd/frameSize_temp,1)),[],2);
        [~,idxEnd] = min(abs(repmat((frameSize_temp:frameSize_temp:timeEnd)',1,length(time))-repmat(time,timeEnd/frameSize_temp,1)),[],2);

        if length(idxEnd) < floor(time(end)/frameSize_temp)
            idxBeg = [idxBeg; idxEnd(end)+1];
            idxEnd = [idxEnd; length(time)];
        end
        X = zeros(size(spec{1},1),length(idxEnd));
        
        for ii = 1:length(idxBeg)
            X(:,ii) = rms(spec{1}(:,idxBeg(ii):idxEnd(ii)-1),2);
        end
end

Xscale = linspace(0,setting.sr/2,size(X,1));
W = X./repmat(sum(X),size(X,1),1);

if indClass == 1
    indTraffic = ones(1,size(W,2));
else
    indTraffic = zeros(1,size(W,2));
end