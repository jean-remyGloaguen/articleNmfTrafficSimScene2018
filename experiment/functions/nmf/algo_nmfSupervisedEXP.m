function NMF = algo_nmfSupervisedEXP(H,W,V,iteration,soundMix,setting)

sparsity = setting.sparsity; 
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

if setting.SM_weight~=0
    switch setting.smoothnessForm
        case 'all'
            smoothness(1:size(W,2),1) = setting.SM_weight;
        case 'traffic'
            smoothness = zeros(size(W,2),1);
            smoothness(soundMix.indTraffic==1) = setting.SM_weight;
    end
else
    smoothness = 0;
end

switch beta
    case 2
        for iter = 1:iteration(end)
            if smoothness==0
                H = H.*(W'*V-sparsity)./(W'*Vap);
                H(H<0)=0;
            else
                H_1 = [zeros(K,1) H(:,1:end-1)];
                H_2 = [H(:,2:end) zeros(K,1)];
                H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                H = H.*(W'*V-sparsity+2*smoothness.*(H_1+H_2))./(W'*Vap+2*smoothness.*(H+H_12));
            end
            H(isnan(H)) = 0;
            Vap = W*H;
        end
    case 1
        for iter = 1:iteration
            
            if smoothness==0
                H = H.*((W'*(V.*Vap.^(-1)))./(sparsity+W'*Vap.^(0)));
            else
                H_1 = [zeros(K,1) H(:,1:end-1)];
                H_2 = [H(:,2:end) zeros(K,1)];
                H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                H = H.*((W'*(V.*Vap.^(-1))+2*smoothness.*(H_1+H_2))./(sparsity+W'*Vap.^(0)+2*smoothness.*(H+H_12)));
            end
            H(isnan(H)) = 0;
            Vap = W*H;
            
        end
    case 0
        for iter = 1:iteration
            
            if smoothness==0
                H = H.*((W'*(V.*Vap.^(-2)))./(sparsity+W'*Vap.^(-1))).^(1/2);
            else
                H_1 = [zeros(K,1) H(:,1:end-1)];
                H_2 = [H(:,2:end) zeros(K,1)];
                H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                H = H.*((W'*(V.*Vap.^(-2))+2*smoothness.*(H_1+H_2))./(sparsity+W'*Vap.^(-1))+2*smoothness.*(H+H_12)).^(1);
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
                    H_1 = [zeros(K,1) H(:,1:end-1)];
                    H_2 = [H(:,2:end) zeros(K,1)];
                    H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                    H = H.*((W'*(V.*Vap.^(beta-2))+2*smoothness.*(H_1+H_2))./(sparsity+W'*Vap.^(beta-1))+2*smoothness.*(H+H_12)).^(gamma);
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
                    H_1 = [zeros(K,1) H(:,1:end-1)];
                    H_2 = [H(:,2:end) zeros(K,1)];
                    H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                    H = H.*((W'*(V.*Vap.^(beta-2))-sparsity+2*smoothness.*(H_1+H_2))./(W'*Vap.^(beta-1))+2*smoothness.*(H+H_12)).^(gamma);
                end
                H(isnan(H)) = 0;
                Vap = W*H;

            end
        end
end
NMF.H = H;
NMF.Vap = Vap;