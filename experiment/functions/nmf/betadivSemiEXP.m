function cost = betadivSemiEXP(V,NMF,setting,varargin)


beta = setting.beta;
sparsity = setting.sparsity;
smoothness = setting.SM_weight;
betaM = setting.SS_betaM;
mu = setting.SS_weight;
lambda = setting.SS_sensitivity;

H = NMF.H;
W = NMF.W;
Y = NMF.Y;
Vtot = NMF.Vtot;

term = zeros(1,size(W,2));
if mu == 0
    cost = betadivEXP(V,Vtot,beta,H,sparsity,smoothness);
else

    div = betadivEXP(Y(:)*ones(1,size(W,2)),repmat(W,size(Y,2),1),betaM);
	Cm = exp((-1/lambda)*sum(div));

    switch beta
        case 2
            cost = sum((V(:)-Vtot(:)).^2)/2+mu*sum(term);
        case 1
            ind_0 = find(V(:)<=eps);
            ind_1 = 1:length(V(:));
            ind_1(ind_0) = [];
            cost = sum(V(ind_1).*log(V(ind_1)./Vtot(ind_1))-V(ind_1)+Vtot(ind_1))+sum(Vtot(ind_0))...
                +mu*Cm;
        case 0
            ind_0 = find(V(:)<=eps);
            ind_1 = 1:length(V(:));
            ind_1(ind_0) = [];

            cost = sum(V(ind_1)./Vtot(ind_1)-log(V(ind_1)./Vtot(ind_1)))-length(V(:))+sum(Vtot(ind_0))...
                +mu*Cm;

        otherwise
            cost = sum(V(:).^beta + (beta-1)*Vtot(:).^beta - beta*V(:).*Vtot(:).^(beta-1) )/(beta*(beta-1))...
                +mu*Cm;
    end
end
if sparsity~=0
    cost = cost+sparsity*sum(H(:));
end

if smoothness~=0
    if strcmp(setting.smoothnessForm,'all')
        basis = [H;NMF.Z];
        cost = cost+0.5*smoothness*sum((sum((basis(:,1:end-1)-basis(:,2:end)).^2)));
    else
        cost = cost+0.5*smoothness*sum((sum((H(:,1:end-1)-H(:,2:end)).^2)));
    end
end