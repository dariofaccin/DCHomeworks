clc; close all; clear global; clearvars;

% Input
load('Useful.mat');
SNR_vect = 8:14;
Pe_VA = zeros(length(SNR_vect),1);
errors = zeros(length(SNR_vect),1);

sigma_a = 2;
M = 4;

gm = conj(qc(end:-1:1));
h = conv(qc,gm);
t0_bar = find(h == max(h));
h = h(h>max(h)/100);
h = h(3:end-2);
h_T = downsample(h,4);
r_gm = xcorr(gm,gm);

in_bits_2  =  in_bits(1+4-0 : end-4+4-2);

for i=1:length(SNR_vect)
    snr_db = SNR_vect(i);
    snr_lin = 10^(snr_db/10);
    [r_c, sigma_w, qc] = channel_sim(in_bits, snr_db, sigma_a);
    r_c = r_c + w(:,i);
    
    % Filtering received signal
    r_c_prime = filter(gm,1,r_c);
    
    r_c_prime = r_c_prime(t0_bar:end);
    x = downsample(r_c_prime,4);
    rw_tilde = sigma_w/4 .* downsample(r_gm, 4);
    M1 = 5;
	N2 = 2;
	D = 2;
	M2 = N2 + M1 - 1 - D;
	[c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);

	psi = conv(c_opt, h_T);
	psi = psi/max(psi);
	
	y = conv(x, c_opt);
	y = y/max(psi);
	
	detected = VBA(y, psi, 0, 4, 4, 4);
	detected = detected(D+1:end);

	[Pe_VA(i),errors(i)] = SER(in_bits_2(1:length(detected)), detected);
end


figure();
semilogy(SNR_vect, Pe_VA, 'r--');
hold on; grid on;

save('Pe_VA.mat','Pe_VA'),