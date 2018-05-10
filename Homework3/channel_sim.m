function [output, sigma_w, qc] = channel_sim(x, snr, sigma_a)

T = 1;
Tc = T/4;
Q = T/Tc;

alpha = 0.67;
beta = 0.7424;

snr_db = snr;
snr_lin = 10^(snr_db/10);

qc_num = [0 0 0 0 0 beta];
qc_denom = [1 -alpha];
qc = impz(qc_num, qc_denom);
aaaaa = find(qc>(max(qc)/100));
qc = qc(1:(aaaaa(end))+1);
E_qc = sum(qc.^2);

sigma_w = sigma_a * E_qc / snr_lin;

a_prime = upsample(x,Q);

s_c = filter(qc_num, qc_denom, a_prime);

% noise = wgn(length(s_c),1,sigma_w,'complex');

r_c = s_c;
output = r_c;