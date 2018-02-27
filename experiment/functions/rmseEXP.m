function errorRMSE = rmseEXP(X,Y,varargin)
%% rmse(X,Y,DIMENSION)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(X,2)>size(Y,2)
    X = X(:,1:size(Y,2));
elseif size(X,2)<size(Y,2)
    Y = Y(:,1:size(X,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FONCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ~isfinite(X)
%     X(ii) = 0;
% elseif ~isfinite(Y)
%     Y(ii) = 0;
% end
if size(X,1)~=1
    errorRMSE = sqrt(mean(bsxfun(@minus,Y,X).^2));
else
    errorRMSE = sqrt(sum((X-Y).^2));
end
