clc; close all; clear global; clearvars;

load('Useful.mat');

SNR_vect = [8 11 14];

sigma_a = 2;
M = 4;

gm = conj(qc(end:-1:1));
h = conv(qc,gm);
h = h(h>max(h)/100);
h = h(3:end-2);
h_T = downsample(h,4);
t0_bar = length(gm);
r_gm = xcorr(gm,gm);

N2 = floor(length(h_T)/2);
D = 2;
Pe_DFE = zeros(10, length(SNR_vect));

for M1=1:10
	for snr=1:length(SNR_vect)
		channel_snr_db = SNR_vect(snr);
		[r_c, sigma_w, ~] = channel_sim(in_bits, channel_snr_db, sigma_a);
		r_c = r_c + w(:,snr);
		M2 = N2 + M1 - 1 - D;
		r_w = sigma_w/4 .* downsample(r_gm, 4);
		
		r_c_prime = filter(gm,1,r_c);
    
		r_c_prime = r_c_prime(t0_bar:end);
		x = downsample(r_c_prime,4);
		[c_opt, JminDFE] = Adaptive_DFE(h_T, r_w, sigma_a, M1, M2, D);
		
		psi = conv(c_opt, h_T);
		psi = psi/max(psi);

		b = - psi(end - M2 + 1:end);
		detected = equalization_DFE(x, c_opt, b, M1, M2, D);

		errors = length(find(in_bits(1:length(detected))~=detected));
		Pe_DFE(M1,snr) = errors/length(in_bits(1:length(detected)));
	end
end

[min, idx] = min(Pe_DFE(:));

[idx_snr, idx_m1] = ind2sub(size(Pe_DFE), idx);

% save('JminDFE.mat', 'JminDFE', 'idx_m1', 'idx_d', 'M2');