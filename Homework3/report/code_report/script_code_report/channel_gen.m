clc; close all; clear global; clearvars;

% MUST BE USED ONCE TO GENERATE AN INPUT FOR
% ALL RECEIVERS

set(0,'defaultTextInterpreter','latex')    % latex format
T = 1;              % Symbol period
Tc = T/4;           % upsampling period
Q = T/Tc;           
snr_db = 10;
snr_lin = 10^(snr_db/10);
sigma_a = 2;

% Filter setup
alpha = 0.67;
beta = 0.7424;
qc_num = [0 0 0 0 0 beta];
qc_denom = [1 -alpha];
qc = impz(qc_num, qc_denom);
q_c = impz(qc_num, qc_denom);
q_c = [0; 0; 0; 0; 0; q_c( q_c >= max(q_c)*10^(-2) )]; 
% remove all components under max(qc/100)
aaaaa = find(qc>(max(qc)/100));
qc = qc(1:(aaaaa(end))+1);
% stem(0:length(qc)-1,qc);
E_qc = sum(qc.^2);              % Energy of the impulse response

% Input signal: from a PN sequence generate QPSK
L = 1023;
x = PNSeq(L);
in_bits = bitmap(x(1:length(x)-1));

% Upsampling
a_prime = upsample(in_bits,Q);

% Noise
sigma_w = sigma_a * E_qc / snr_lin;     % N0

% Output
s_c = filter(qc_num, qc_denom, a_prime);
wc = wgn(length(s_c),1,sigma_w,'complex');
r_c = s_c + wc;

%% FIGURES
abs_qc = abs(qc);

figure()
stem(abs_qc);
xlabel('$m\frac{T}{4}$');
ylabel('$q_c$')
xlim([1 length(abs_qc)]);
grid on

[Qc, f] = freqz(qc,1,1024,'whole');
% f = f/(2*pi);
f = linspace(0,4,length(f));

figure()
plot(f,10*log10(abs(Qc)));
xlim([0 2]);
xlabel('f'), grid on
ylabel('$| Q_c(f) |$ $[dB]$')

grid on