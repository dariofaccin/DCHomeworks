clc; close all; clear global; clearvars;

T = 1;              % Symbol period
Tc = T/4;           % upsampling period
Q = T/Tc;           

alpha = 0.67;       % filter coefficients
beta = 0.7424;

qc_num = [0 0 0 0 0 beta];       % qc filter
qc_denom = [1 -alpha];
qc = impz(qc_num, qc_denom);
E_qc = sum(qc.^2);

x = PNSeq(1023);     % generate a strem of 0 and 1 

% x_i = x(1:ceil(length(x)/2)-1);
% x_q = x(ceil(length(x)/2)+1:end);

x_i = x(1:2:end);
x_q = x(2:2:end-1);

x_i = x_i(1:length(x_q));

a_prime_i = upsample(x_i, Q);
a_prime_q = upsample(x_q, Q);

a_prime = a_prime_i + 1i*a_prime_q;

s_c = filter(qc_num, qc_denom, a_prime);

noise = wgn(length(s_c),1,10);

r_c = s_c+noise;

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

save('channel_output.mat', 'r_c');