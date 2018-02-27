function NMF = algo_nmfSemiSupervisedEXP(H,W,Y,Z,V,iteration,soundMix,setting)

beta = setting.beta;
sparsity = setting.sparsity;
smoothness = setting.SM_weight;
J = setting.SS_sizeWrand;
K = size(W,2);

if setting.SS_weight ~= 0
    mu = setting.SS_weight;
    betaM = setting.SS_betaM;
    lambdaSS = setting.SS_lambda;
else
    mu = 0;
end

if ~isequal(beta,0) || ~isequal(beta,1) || ~isequal(beta,2)
    if beta < 1
        gamma = 1/(2-beta);
    elseif beta >= 1 && beta <= 2
        gamma = 1;
    else
        gamma = 1/(beta-1);
    end
end

if smoothness~=0
    smoothForm = setting.smoothnessForm;
    
    switch smoothForm
        case 'all'
            lambdaSM = [sum(W,1) sum(Y,1)];
            smoothnessTot(1:K+J,1) = smoothness;
        case 'traffic'
            lambdaSM = sum(W,1);
            smoothnessTot = zeros(K,1);
            smoothnessTot(soundMix.indTraffic==1) = smoothness;
    end
else
    smoothnessTot = 0;
    lambdaSM = NaN;
end

Vap = W*H;
Vtot = Vap+Y*Z;

NMF = [];

switch beta
    case 2
        for iter = 1:iteration
            if mu == 0
                Y = Y.*(((V*Z')./(Vtot*Z')).^gamma);
            else
                div = betadivEXP(Y(:)*ones(1,size(W,2)),repmat(W,size(Y,2),1),betaM);
                Cm = exp((-1/lambdaSS)*sum(div));

                num = lambdaSS*V*Z'+mu*Y.^(betaM-1)*Cm;
                den = (lambdaSS*Vtot)*Z'+mu*Y.^(betaM-2).*Cm.*sum(W,2);
                Y = Y.*(num./den);
                Y(Y<2e-5) = 2e-5;
            end
            Y = Y./repmat(sum(Y),size(Y,1),1);
            Vtot = Vap+Y*Z;
            
            if smoothness == 0 || strcmp(smoothForm,'traffic')
                Z = Z.*(((Y'*V)./(Y'*Vtot)).^gamma);
            elseif smoothness == 1 && strcmp(smoothForm,'all')
                Z = updateSmoothEssidEXP(Z,Y,V,Vtot,beta,smoothnessTot,sparsity,lambdaSM(K+1:end));
            end
            Vtot = Vap+Y*Z;
            
            if smoothness == 0
                H = H.*(((W'*V-sparsity)./(W'*Vtot)).^gamma);
                H(H<0)=0;
            else
                H = updateSmoothEssidEXP(H,W,V,Vtot,beta,smoothnessTot,sparsity,lambdaSM(1:K));
            end
            
            Vap = W*H;
            Vtot = Vap+Y*Z;
        end
        
    case 1
        F = size(V,1);
        N = size(H,2);
        
        for iter = 1:iteration
            if mu == 0
                Y = Y.*((((V.*(Vtot).^(-1))*Z')./repmat(sum(Z,2)',F,1)).^gamma);
            else
                div = betadivEXP(Y(:)*ones(1,size(W,2)),repmat(W,size(Y,2),1),betaM);
                Cm = exp((-1/lambdaSS)*sum(div));

                num = lambdaSS*V.*(Vtot).^(-1)*Z'+mu*Y.^(betaM-1)*Cm;
                den = (lambdaSS*repmat(sum(Z,2)',F,1))+mu*Y.^(betaM-2).*Cm.*sum(W,2);
                Y = Y.*(num./den);
                Y(Y<2e-5) = 2e-5;
            end
            
            Y = Y./repmat(sum(Y), size(Y, 1), 1);
            Vtot = Vap+Y*Z;
            
            Z = Z.*(((Y'*(V.*(Vtot).^(-1)))./repmat(sum(Y,1)',1,N)).^gamma);
            Vtot = Vap+Y*Z;
            
            if smoothness == 0
                H = H.*(((W.'*(V.*(Vtot).^(-1)))./(sparsity+repmat(sum(W,1)',1,N))).^gamma);
            else
                H = updateSmoothEssidEXP(H,W,V,Vtot,beta,smoothnessTot,sparsity,lambdaSM);
            end
            
            Vap = W*H;
            Vtot = Vap+Y*Z;
        end
    case 0
        for iter = 1:iteration
            if mu == 0
                Y = Y.*((((V.*Vtot.^(-2))*Z')./(Vtot.^(-1)*Z')).^gamma);
            else
                div = betadivEXP(Y(:)*ones(1,K),repmat(W,J,1),betaM);
                Cm = exp((-1/lambdaSS)*(sum(div)));
                
                num = lambdaSS*V.*(Vtot).^(-2)*Z'+mu*Y.^(betaM-1)*Cm;
                den = lambdaSS*Vtot.^(-1)*Z'+mu*Y.^(betaM-2).*Cm.*repmat(sum(W,2),1,J);
                Y = Y.*((num./den).^gamma);
                Y(Y<eps) = eps;
            end
            
            Y = Y./repmat(sum(Y), size(Y, 1), 1);
            Vtot = Vap+Y*Z;
            
            Z = Z.*(((Y'*(Vtot.^(-2).*V))./(Y'*Vtot.^(-1))).^gamma);
            Vtot = Vap+Y*Z;
            
            if smoothness == 0
                H = H.*(((W'*(Vtot.^(-2).*V))./(sparsity+W'*Vtot.^(-1))).^gamma);
            else
                H = updateSmoothEssidEXP(H,W,V,Vtot,beta,smoothnessTot,sparsity,lambdaSM);
            end
            
            Vap = W*H;
            Vtot = Vap+Y*Z;
        end
    otherwise
        for iter = 1:iteration
            if mu == 0
                Y = Y.*((((V.*Vtot.^(beta-2))*Z')./(Vtot.^(beta-1)*Z')).^gamma);
            else
                div = betadivEXP(Y(:)*ones(1,K),repmat(W,J,1),betaM);
                Cm = exp((-1/lambdaSS)*(sum(div)));
                
                num = V.*(Vtot).^(beta-2)*Z'+mu*Y.^(betaM-1)*Cm;
                den = Vtot.^(beta-1)*Z'+mu.*Y.^(-1).*repmat(K,1,J);
                Y = Y.*((num./den).^gamma);
            end
            
            Y = Y./repmat(sum(Y), size(Y, 1), 1);
            Vtot = Vap+Y*Z;
            
            Z = Z.*(((Y'*(Vtot.^(beta-2).*V))./(Y'*Vtot.^(beta-1))).^gamma);
            Vtot = Vap+Y*Z;
            
            if smoothness == 0
                if beta < 2
                    H = H.*(((W'*(Vtot.^(beta-2).*V))./(sparsity+W'*Vtot.^(beta-1))).^gamma);
                else
                    H = H.*(((W'*(Vtot.^(beta-2).*V)-sparsity)./(W'*Vtot.^(beta-1))).^gamma);
                end
            else
                H = updateSmoothEssidEXP(H,W,V,Vtot,beta,smoothnessTot,sparsity,lambdaSM);
            end
            Vap = W*H;
            Vtot = Vap+Y*Z;
            
        end
end

NMF.H = H;
NMF.W = W;
NMF.Z = Z;
NMF.Y = Y;
NMF.Vtot = Vtot;