function X = resolutionThirdOrder(a,b,c,d,alphaCoeff)

K = length(a);

if length(alphaCoeff)==1
    alphaCoeff(1:K,1) = alphaCoeff;
end

alphaNonNull = find(alphaCoeff~=0);
alphaNull = find(alphaCoeff==0);

%% alphaNonNull ~= 0
if ~isempty(alphaNonNull)
    p = c./a(alphaNonNull)-(b(alphaNonNull).^2)./(3*(a(alphaNonNull).^2));
    q = (2*(b(alphaNonNull).^3))./(27*(a(alphaNonNull).^3))+...
        d(alphaNonNull)./a(alphaNonNull)-b(alphaNonNull)*c./(3*a(alphaNonNull).^2);

    delta = q.^2+4*((p./3).^3);

    pos = find(delta > 0);
    nul = find(delta == 0);
    neg = find(delta <0);

    X = zeros(K,1);
    if ~isempty(pos)
        C1 = nthroot((-q(pos)+sqrt(delta(pos)))/2,3);
        C2 = nthroot((-q(pos)-sqrt(delta(pos)))/2,3);

        X(pos) = C1+C2-(b(pos)./(3*a(pos)));
    end
    if ~isempty(nul)

        Xnul(:,1) = (3*q(nul))./p(nul)-(b(nul)./(3*a(nul)));
        Xnul(:,2) = -(3*q(nul))./(2*p(nul))-(b(nul)./(3*a(nul)));
        Xnul(:,3) = Xnul(:,2);
        X(nul) = min(abs(Xnul),[],2);
    end
    if~isempty(neg)
        deltaLoop = abs(delta(neg));
        jj = -0.5+1j*(sqrt(3)/2);
        c1 = (-q(neg)+1j*sqrt(deltaLoop))/2;
        c2 = (-q(neg)-1j*sqrt(deltaLoop))/2;
        C1 = nthroot(abs(c1),3).*exp(1j*angle(c1)/3);
        C2 = nthroot(abs(c2),3).*exp(1j*angle(c2)/3);

        Xneg(:,1) = C1+C2-(b(neg)./(3*a(neg)));
        Xneg(:,2) = jj*C1+(jj)^2*C2-(b(neg)./(3*a(neg)));
        Xneg(:,3) = (jj)^2*C1+jj*C2-(b(neg)./(3*a(neg)));
        X(neg) = min(abs(Xneg),[],2);
    end
end
%% alphaNull == 0
if ~isempty(alphaNull)
    X(alphaNull) = sqrt(-d(alphaNull)./b(alphaNull));
end


