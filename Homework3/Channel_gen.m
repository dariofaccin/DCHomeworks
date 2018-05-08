clc; close all; clear global; clearvars;

% MUST BE USED ONCE TO GENERATE AN INPUT FOR
% ALL RECEIVERS

T = 1;              % Symbol period
Tc = T/4;           % upsampling period
Q = T/Tc;           
snr = 10;
sigma_a = 2;

% Filter setup
alpha = 0.67;
beta = 0.7424;
qc_num = [0 0 0 0 0 beta];
qc_denom = [1 -alpha];
qc = impz(qc_num, qc_denom);
% remove al components under max(qc/100)
qc = qc([1:17]);
stem([0:length(qc)-1],qc);

E_qc = sum(qc.^2);

% Input signal: from a PN sequence generate QPSK
L = 1023;
x = PNSeq(L);
in_bits = bitmap(x(1:length(x)-1));

% Upsampling
a_prime = upsample(in_bits,Q);

% Noise
sigma_w = sigma_a * E_qc / (4*snr);     % N0

% Output
s_c = filter(qc_num, qc_denom, a_prime);
wc = wgn(length(s_c),1,sigma_w);
r_c = s_c+wc;

save('rec_input.mat', 'in_bits', 'r_c', 'qc', 'E_qc', 'wc','sigma_w');


%% FIGURES
abs_qc = abs(qc);

figure()
stem(abs_qc);
xlabel('n T/4');
xlim([1 length(abs_qc)]);

[Qc, f] = freqz(qc,1,4096,'whole');

figure()
plot(f,10*log10(abs(Qc)/length(Qc)));
xlim([0 2/T]);

