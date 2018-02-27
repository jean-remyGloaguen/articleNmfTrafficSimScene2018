function NMF = algo_nmfSupervisedEXP(H,W,V,iteration,soundMix,setting)

sparsity = setting.sparsity; 
smoothness = setting.SM_weight; 
beta = setting.beta;
Vap = W*H;

NMF = [];

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
    lambda = sum(W,1);
    smoothForm = setting.smoothnessForm;
    
    switch smoothForm
        case 'all'
            smoothnessTot(1:size(W,2),1) = smoothness;
        case 'traffic'
            smoothnessTot = zeros(size(W,2),1);
            smoothnessTot(soundMix.indTraffic==1) = smoothness;
    end
else
    smoothnessTot = 0;
    lambda = NaN;
end

switch beta
    case 2
        for iter = 1:iteration(end)
            if smoothness==0
                H = H.*(W'*V-sparsity)./(W'*Vap);
                H(H<0)=0;
            else
                H = updateSmoothEssidEXP(H,W,V,Vap,beta,smoothnessTot,sparsity,lambda);
            end
            H(isnan(H)) = 0;
            Vap = W*H;
        end
    case 1
        for iter = 1:iteration
            
            if smoothness==0
                H = H.*((W'*(V.*Vap.^(-1)))./(sparsity+W'*Vap.^(0)));
            else
                H = updateSmoothEssidEXP(H,W,V,Vap,beta,smoothnessTot,sparsity,lambda);
            end
            H(isnan(H)) = 0;
            Vap = W*H;
            
        end
    case 0
        for iter = 1:iteration
            
            if smoothness==0
                H = H.*((W'*(V.*Vap.^(-2)))./(sparsity+W'*Vap.^(-1))).^(1/2);
            else
                H = updateSmoothEssidEXP(H,W,V,Vap,beta,smoothnessTot,sparsity,lambda);
            end
            H(isnan(H)) = 0;
            Vap = W*H;

        end
    otherwise
        if beta < 1
            for iter = 1:iteration
                if smoothness==0
                    H = H.*((W'*(V.*Vap.^(beta-2)))./(sparsity+W'*Vap.^(beta-1))).^(gamma);
                else
                    H = updateSmoothEssidEXP(H,W,V,Vap,beta,smoothnessTot,sparsity,lambda);
                end
                H(isnan(H)) = 0;
                Vap = W*H;
                
            end
        else
            for iter = 1:iteration
                
                if smoothness==0
                    H = H.*((W'*(V.*Vap.^(beta-2))-sparsity)./(W'*Vap.^(beta-1))).^(gamma);
                    H(H<0)=0;
                else
                    H = updateSmoothEssidEXP(H,W,V,Vap,beta,smoothness,sparsity,lambda);
                end
                H(isnan(H)) = 0;
                Vap = W*H;

            end
        end
end
NMF.H = H;
NMF.Vap = Vap;