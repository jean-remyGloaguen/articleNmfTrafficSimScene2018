function [Y,binMel,wts] = spectre2MelEXP(X,numberMel,cutOffFreq,sr)

if iscell(X)
    Y = cell(size(X,1),size(X,2));
    for ii = 1:size(X,1)
        for jj = 1:size(X,2)
            [Y{ii,jj},wts,binMel] = audspec(X{ii,jj}, sr, numberMel,'mel', 20, cutOffFreq, 0, 1);
        end
    end
else 
    [Y,wts,binMel] = audspec(X, sr, numberMel,'mel', 20, cutOffFreq, 0, 1);

end