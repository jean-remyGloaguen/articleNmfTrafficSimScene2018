function [Lp,Leq] = estimationLpEXP(Xinput,setting)

if ~iscell(Xinput)
    Xin{1} = Xinput;
else
    Xin = Xinput;
end

sr = setting.sr;
window = setting.window;
noverlap = setting.noverlap;
step = (window*(noverlap/100))/sr;

Lp = cell(1,size(Xin,2));
Leq = zeros(1,size(Xin,2));

for ii = 1:size(Xin,2)
    X = Xin{ii};
    X(2:end,:) = X(2:end,:)*2;
    time =  step:step:size(X,2)*step;
    
    B = sqrt(abs(X).^2);
    B(2:end,:) = sqrt(2)*B(2:end,:);
    Bw = B/sqrt(1.5);           %hanning window effect "scaling effect" see https://www.ap.com/technical-library/fft-scaling-for-noise/

    a = sqrt(sum(Bw.^2));
    Leq(ii) = rms(a);
    
    T_tau = floor(time(end)/setting.temporalInterval); % nombre interval
    N_tau = floor(length(time)/T_tau); % nombre echantilon par interval
    ind_beg = linspace(0,(T_tau-1)*N_tau,T_tau);
    ind_end = linspace(N_tau,T_tau*N_tau,T_tau);
    
    Lp{ii} = zeros(1,T_tau);
    for jj = 1:T_tau
        Lp{ii}(:,jj) = rms(a(ind_beg(jj)+1:ind_end(jj)));
    end
    Lp{ii}(isnan(Lp{ii}))=0;
    Lp{ii}(isinf(Lp{ii}))=0;
end
