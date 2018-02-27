function [NMF] = algo_nmfUnsupervisedEXP(H,W,W0,V,iteration,setting)

sparsity = setting.sparsity;
beta = setting.beta;
Vap = W*H;
K = size(W,2); 

if beta == 0
    smoothness = setting.smoothnessBeta0;
elseif beta == 1
    smoothness = setting.smoothnessBeta1;
elseif beta == 2
    smoothness = setting.smoothnessBeta2;
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

cost = zeros(1,iteration);
switch beta
    case 2
        for iter = 1:iteration
            
            W = W.*(V*H')./(Vap*H');
            W = (1-setting.weightW)*W0+setting.weightW*W;
            W(isnan(W)) = 0;
            W = W./repmat(sum(W),size(W,1),1);
            Vap = W*H;
            
            if smoothness == 0
                H = H.*(W'*V-sparsity)./(W'*Vap);
            else
                H_1 = [zeros(K,1) H(:,1:end-1)];
                H_2 = [H(:,2:end) zeros(K,1)];
                H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                H = H.*(W'*V-sparsity+2*smoothness.*(H_1+H_2))./(W'*Vap+2*smoothness.*(H+H_12));
            end
            
            H(isnan(H)) = 0;

            Vap = W*H;
            cost(iter) = betadivEXP(V,Vap,beta);
        end
    case 1
        for iter = 1:iteration
            W = W .* ((V./Vap)*H')./(repmat(sum(H,2)',size(V,1),1));
            W(isnan(W)) = 0;
            W = W./repmat(sum(W),size(W,1),1);
            Vap = W*H;

            if smoothness == 0
                H = H.*((W'*(V./Vap))./(sparsity+W'*Vap.^(0)));
            else
                H_1 = [zeros(K,1) H(:,1:end-1)];
                H_2 = [H(:,2:end) zeros(K,1)];
                H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                H = H.*((W'*(V.*Vap.^(-1))+2*smoothness.*(H_1+H_2))./(sparsity+W'*Vap.^(0)+2*smoothness.*(H+H_12)));
            end
            H(isnan(H)) = 0;
            
            Vap = W*H;
            cost(iter) = betadivEXP(V,Vap,beta);
        end
    case 0
        for iter = 1:iteration
            
            W = W .* (((V.*Vap.^(-2))*H')./(Vap.^(-1)*H')).^gamma;
            W = W./repmat(sum(W),size(W,1),1);
            Vap = W*H;
            
            if smoothness == 0
                H = H.*((W'*(V.*Vap.^(-2)))./(sparsity+W'*Vap.^(-1))).^(1/2);
            else
                H_1 = [zeros(K,1) H(:,1:end-1)];
                H_2 = [H(:,2:end) zeros(K,1)];
                H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                H = H.*((W'*(V.*Vap.^(-2))+2*smoothness.*(H_1+H_2))./(sparsity+W'*Vap.^(-1))+2*smoothness.*(H+H_12)).^(1);
            end
            H(isnan(H)) = 0;
            
            Vap = W*H;
            cost(iter) = betadivEXP(V,Vap,beta);
        end
    otherwise
        if beta < 1
            for iter = 1:iteration
                if smoothness==0
                    H = H.*((W'*(V.*Vap.^(beta-2)))./(sparsity+W'*Vap.^(beta-1))).^(gamma);
                else
                    H = updateSmoothEssid(H,W,V,Vap,beta,smoothness,sparsity,lambda);
                    H(isnan(H)) = 0;
                end
                
                Vap = W*H;
                cost(iter) = betadivEXP(V,Vap,setting.beta);
            end
        else
            for iter = 1:iteration
                
                if smoothness==0
                    H = H.*((W'*(V.*Vap.^(beta-2))-sparsity)./(W'*Vap.^(beta-1))).^(gamma);
                    H(H<0)=0;
                else
                    H = updateSmoothEssid(H,W,V,Vap,beta,smoothness,sparsity,lambda);
                    H(isnan(H)) = 0;
                end
                
                Vap = W*H;
                cost(iter) = betadivEXP(V,Vap,setting.beta);
            end
        end
end

NMF.Vap = Vap;
NMF.H = H;
NMF.W = W;
NMF.cost = cost;