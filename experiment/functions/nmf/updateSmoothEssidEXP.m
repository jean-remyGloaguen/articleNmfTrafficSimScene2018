function [H] = updateSmoothEssidEXP(H,W,V,V_ap,beta,smoothness,sparsity,lambda)

N = size(H,2);

switch beta
    case 0
        Ht = H;
        a = 2*smoothness.*(lambda'.^2);
        c = 0;
        d = -W'*(V./(V_ap.^2)).*(Ht.^2);
        
        % n = 1
        b = W'*(V_ap(:,1).^(-1))-smoothness.*(lambda'.^2).*(H(:,2));
        H(:,1) = resolutionThirdOrder(a,b,c,d(:,1),smoothness);
        
        % n = 2:N-1
        for n = 2:N-1
            b = W'*(V_ap(:,n).^(-1))-smoothness.*(lambda'.^2).*(H(:,n-1)+H(:,n+1));
            H(:,n) = resolutionThirdOrder(a,b,c,d(:,n),smoothness);
        end
        
        % n = N
        b = W'*(V_ap(:,N).^(-1))-smoothness.*(lambda'.^2).*(H(:,N-1));
        H(:,N) = resolutionThirdOrder(a,b,c,d(:,N),smoothness);
        
        
    case 1
        if any(smoothness==0)
            smoothness(smoothness==0)=1e-12;
        end
        
        Ht = H;
        a_n = 2*smoothness.*(lambda'.^2);
        a_1 = a_n/2;
        a_N = a_n/2;
        bt = sparsity+lambda';
        c = Ht.*(W'*(V./V_ap));
        
        % n = 1
        b(:,1) = bt-smoothness.*(lambda'.^2).*H(:,2);
        H(:,1) = (sqrt(b(:,1).^2+4*a_1.*c(:,1))-b(:,1))./a_1;
        
        % n = 2:N-1
        b(:,2:N-1) = repmat(bt,1,size(H,2)-2)-repmat(a_n,1,size(H,2)-2).*(H(:,1:N-2)+H(:,3:N));
        for n = 2:N-1
            H(:,n) = (sqrt(b(:,n).^2+4*a_n.*c(:,n))-b(:,n))./(2*a_n);
        end
        
        % n = N
        b(:,N) = bt-smoothness.*(lambda'.^2).*(H(:,N-1));
        H(:,N) = (sqrt(b(:,N).^2+4*a_N.*c(:,N))-b(:,N))./a_N;
        
        H(H<0)=0;
        
    case 2

        if any(smoothness==0)
            smoothness(smoothness==0)=1e-12;
        end
        
        Ht = H;
        lambda = sum(W)';
        a = (W'*V_ap)./Ht+2*repmat(smoothness.*(lambda.^2),1,size(H,2));
        
        % n = 1
        b = W'*V(:,1)+smoothness.*(lambda.^2).*H(:,2)-sparsity;
        H(:,1) = b./a(:,1);
        
        %         n = 2:N-1
        for n = 2:N-1
            b = W'*V(:,n)+smoothness.*(lambda.^2).*(H(:,n+1)+H(:,n-1))-sparsity;
            H(:,n) = b./a(:,n);
        end
        
        %         n = N
        b = W'*V(:,N)+smoothness.*(lambda.^2).*H(:,N-1)-sparsity;
        H(:,N) = b./a(:,N);
        
        H(H<0)=0;
end
