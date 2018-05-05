clc; close all; clear global; clearvars;

T = 1;
Tc = T/4;

Q = T/Tc;

alpha = 0.67;
beta = 0.7424;

qc_num = [0 0 0 0 0 beta];
qc_denom = [1 -alpha];
qc = impz(qc_num, qc_denom);

E_qc = sum(qc.^2);

x = PNSeq(1023);

x_i = x(1:ceil(length(x)/2)-1);
x_q = x(ceil(length(x)/2)+1:end);

a_prime_i = upsample(x_i, Q);
a_prime_q = upsample(x_q, Q);

a_prime = a_prime_i + 1i*a_prime_q;

s_c = filter(qc_num, qc_denom, a_prime);

noise = wgn(length(s_c),1,10);

r_c = s_c + noise;

%% FIGURES

abs_qc = abs(qc);

figure()
plot(abs_qc);
xlabel('n T/4');
xlim([1 length(abs_qc)]);

Qc = fftshift(fft(qc));

figure()
plot(10*log10(abs(Qc)/length(Qc)));

save('channel_output.mat', 'r_c');