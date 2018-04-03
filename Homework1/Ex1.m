clc; close all; clear global; clearvars;

%% LOAD 1 REALIZATION OF THE PROCESS

load('inputsignal2.mat', 'x');

Nsamples=length(x);

%% SPECTRAL ANALYSIS

% Autocorrelation (unbiased estimate)
[rx_full]=autocorrelation_Unb(x);
L=floor(Nsamples/2);%L should be lower than the length of the r.p. because of the high variance when n approaches K
rx=rx_full(1:L);

% Blackman-Tukey correlogram 
% The length of 2*L+1 is because of page 86 note 24
w_bt=window(@hamming,2*L+1);    
Pbt1=correlogram(x, w_bt, rx, L);
%Pbt2=correlogram(x, w_rect, rx, L);

% Periodogram
X=fft(x);
Pper=(1/Nsamples)*(abs(X)).^2;

% Welch periodogram
S=60;   %overlap
D=120;   %window length
%w_welch=window(@hamming,D);
w_welch=window(@rectwin,D);
%w_welch=kaiser(D,5);
[Welch_P, Ns] = welchPSD(x, w_welch, S);
var_Welch=Welch_P.^2/Ns;


% Analytical PSD: compute the transform of rx(n) on paper and plot it
% according to the requirements
%sigmaw=0.1
b = zeros(1,800);
for i=1:length(b)
    %change to 0.1 for es2
    b(i) = 10*log10(2);
end

% 30 is random choice just to see the plot
b(ceil(0.17*800)) = 10*log10(Nsamples);
b(ceil(0.78*800)) = 0.8*10*log10(Nsamples);

%% Choice of N
N = 25;

[copt, Jmin]=predictor(rx, N);
t=40;
Jvect=zeros(t,1);

for i=1:length(Jvect)
    [c_it, J_it]=predictor(rx, i);
    Jvect(i)=J_it;
end

figure('Name', 'J over N');
plot(1:t,10*log10(Jvect));
title('J_{min} over N');
xlim([1 t]);
xlabel('N'); ylabel('J_{min} [dB]');
% coeff=[1; copt];
% A = tf([1 copt.'], 1,1);
% figure, pzmap(A)

%% AR model
% Coefficients of Wiener filter
[a, s_white, d]=findAR(N, rx);
[H_w, omega] = freqz(1, [1; a], Nsamples, 'whole');

%% Final spectral plot
figure('Name', 'Spectral Analysis');
hold on;
plot((1:Nsamples)/Nsamples, 10*log10(Welch_P), 'r')
plot((1:Nsamples)/Nsamples, 10*log10(abs(Pbt1)), 'b')
plot((1:Nsamples)/Nsamples, 10*log10(Pper), 'g:')
plot(omega/(2*pi), 10*log10(s_white*(abs(H_w)).^2), 'Color', 'm');
plot((1:Nsamples)/Nsamples, b, 'k:');
title('Spectral analysis');
legend('Welch', 'Correlogram', 'Periodogram', ['AR(' int2str(N) ')'], 'Actual value', 'Location', 'SouthWest');
hold off;
xlabel('Normalized frequency');
ylabel('Estimated PSD (dB)');
ylim([-15 30]);
