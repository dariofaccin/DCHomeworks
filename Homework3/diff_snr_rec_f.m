clc; close all; clear global; clearvars;

% Input
load('Useful.mat');
SNR_vect = 8:14;
Pe_FBA = zeros(length(SNR_vect),1);
errors = zeros(length(SNR_vect),1);
sigma_a = 2;
gm = conj(qc(end:-1:1));
h = conv(qc,gm);
h = h(h>max(h)/100);
h = h(3:end-2);
h_T = downsample(h,4);
t0_bar = length(gm);
r_gm = xcorr(gm,gm);
M1 = 5;
N2 = floor(length(h_T)/2);
D = 4;
M2 = N2 + M1 - 1 - D;
L1 = 0;
L2 = 4;
	
for i=1:length(SNR_vect)
    snr_db = SNR_vect(i);
    snr_lin = 10^(snr_db/10);
    [r_c, sigma_w, qc] = channel_sim(in_bits, snr_db, sigma_a);
    r_c = r_c + w(:,i);
    r_c_prime = filter(gm,1,r_c);    
    r_c_prime = r_c_prime(t0_bar:end);
    x = downsample(r_c_prime,4);
    rw_tilde = sigma_w/4 .* downsample(r_gm, 4);
	[c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);
    psi = conv(c_opt, h_T);
    psi = psi/max(psi);
    
    y = conv(x, c_opt);
    y = y/max(psi);
	psi_pad = [psi; 0; 0];
	indexD = find(psi_pad == max(psi_pad));
	detected = FBA(y, psi_pad(indexD:end), L1, L2);
	detected = detected(5:end);
	nerr = length(find(in_bits(1:length(detected))~=detected(1:end)));
	Pe_FBA(i) = nerr/length(detected);
end

figure, semilogy(SNR_vect, Pe_FBA, 'r');

save('Pe_FBA_ATTENZIONE.mat','Pe_FBA');