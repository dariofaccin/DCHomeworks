function [output, sigma_w, qc] = channel_sim(x, snr, sigma_a)

T = 1;
Tc = T/4;
Q = T/Tc;

alpha = 0.67;
beta = 0.7424;

qc_num = [0 0 0 0 0 beta];
qc_denom = [1 -alpha];
qc = impz(qc_num, qc_denom);

E_qc = sum(qc.^2);

sigma_w = sigma_a * E_qc / (4*snr);

x_i = x(1:2:end);
x_q = x(2:2:end-1);
x_i = x_i(1:length(x_q));

a_prime_i = upsample(x_i, Q);
a_prime_q = upsample(x_q, Q);

a_prime = a_prime_i + 1i*a_prime_q;

s_c = filter(qc_num, qc_denom, a_prime);

noise = wgn(length(s_c),1,sigma_w);

r_c = s_c + noise;

output = r_c;