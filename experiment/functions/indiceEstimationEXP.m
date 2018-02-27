function indice = indiceEstimationEXP(x)


indice.Leq = 20*log10(rms(x)/2e-5);

%% LAeq %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = 2/sqrt(2)*abs(fft(x)/length(x));
frequency = linspace(0,44100,length(X));

X = X(frequency>20);
frequency = frequency(frequency>20);
X = X(frequency<20000);
frequency = frequency(frequency<20000);

[X_third,~] = NarrowToNthOctaveEXP(frequency,X,3);

A_weight = [0.0012 0.0021 0.0037 0.0062 0.0098 0.0150 0.0222 0.0313 0.0428...
    0.0570 0.0743 0.0935 0.1151 0.1384 0.1607 0.1824 0.2000 0.2143 0.2244 0.2296...
    0.2323 0.2296 0.2244 0.2119 0.1977 0.1762 0.1500 0.1219 0.0935 0.0686]'.*1e-4;

indice.LAeq = 20*log10(sqrt(sum((X_third(2:end)+A_weight).^2))/2e-5);

%% L50, L10, L90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0:1/44100:length(x)/44100-1/44100;
dt = 0.125;

T_tau = floor(t(end)/dt); % nombre interval
N_tau = floor(length(t)/T_tau); % nombre echantilon par interval
ind_beg = linspace(0,(T_tau-1)*N_tau,T_tau);
ind_end = linspace(N_tau,T_tau*N_tau,T_tau);

tt = linspace(0,t(end),T_tau)+(dt/2);
Lp = zeros(1,T_tau);
for jj = 1:T_tau
    Lp(:,jj) = rms(x(ind_beg(jj)+1:ind_end(jj)));
end

Lp_sort = sort(20*log10(Lp/2e-5),'ascend');
t_10 = (1-0.1)*tt(end);
t_50 = 0.5*tt(end);
t_90 = (1-0.9)*tt(end);

[~,idx_90] = min(abs(tt-t_90));
[~,idx_50] = min(abs(tt-t_50));
[~,idx_10] = min(abs(tt-t_10));

indice.L10 = Lp_sort(idx_10);
indice.L50 = Lp_sort(idx_50);
indice.L90 = Lp_sort(idx_90);

%% Band level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indice.L05 = 20*log10(X_third(15)/2e-5);
indice.L1k = 20*log10(X_third(18)/2e-5);
indice.L2k = 20*log10(X_third(21)/2e-5);
indice.L5k = 20*log10(X_third(25)/2e-5);