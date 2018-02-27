function cost = betadivEXP(V,Vap,beta,varargin)

if nargin > 3
    H = varargin{1};
    sparsity = varargin{2};
    smoothness = varargin{3};
else
    H = [];
    sparsity = 0;
    smoothness = 0;
end

switch beta
    case 2
        cost = sum((V(:)-Vap(:)).^2)/2;
    case 1
        ind_0 = find(V(:)<=eps);
        ind_1 = 1:length(V(:));
        ind_1(ind_0) = [];
        
        cost = sum(V(ind_1).*log(V(ind_1)./Vap(ind_1)) - V(ind_1) + Vap(ind_1) ) + sum(Vap(ind_0));        
    case 0
        ind_0 = find(V(:)<=eps);
        ind_1 = 1:length(V(:));
        ind_1(ind_0) = [];
        
        cost = sum(V(ind_1)./Vap(ind_1)-log(V(ind_1)./Vap(ind_1)))-length(V(:))+sum(Vap(ind_0));

    otherwise
        cost = sum(V(:).^beta+(beta-1)*Vap(:).^beta-beta*V(:).*Vap(:).^(beta-1))/(beta*(beta-1));
end

if sparsity~=0
    cost = cost+sparsity*sum(H(:));
end

if smoothness~=0
    cost = cost+0.5*smoothness*sum(sum((H(:,1:end-1)-H(:,2:end)).^2));
end
    

