clc; close all; clear global; clearvars;

load('Useful.mat', 'qc');
load('Input_symbols.mat');

symbols = bitmap(x).';

SNR_vect = 0:14;
sigma_a = 2;	% Input variance
M = 4;			% Constellation cardinality

gm = conj(qc(end:-1:1));		% Matched filter: complex conjugate of qc
h = conv(qc,gm);				% Impulse response

t0_bar = find(h == max(h));		% Timing phase: peak of h
h = h(h>max(h)/100);
h = h(3:end-2);
h_T = downsample(h,4);

r_gm = xcorr(gm,gm);			% Matched filter autocorrelation
Pbit_DFE_uncode = zeros(length(SNR_vect),1);
N2 = floor(length(h_T)/2);
N1 = N2;
M1 = 5;
D = 4;
M2 = N2 + M1 - 1 - D;
tic
for i=1:length(SNR_vect)
	snr_db = SNR_vect(i);
	snr_lin = 10^(snr_db/10);
	% Single carrier channel simulation
	[r_c, sigma_w, qc] = channel_sim(symbols, snr_db, sigma_a);
	% additive complex-Gaussian noise
	w = wgn(length(r_c),1, 10*log10(sigma_w), 'complex');
	r_c = r_c + w;
	r_c_prime = filter(gm,1,r_c);
	r_c_prime = r_c_prime(t0_bar:end);
	% signal at the input of the DFE
	x_aa = downsample(r_c_prime,4);
	rw_tilde = sigma_w/4 .* downsample(r_gm, 4);
	% DFE
	[c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);
	psi = conv(c_opt, h_T);
	b = - psi(end - M2 + 1:end);
	y_hat = x_aa/max(psi);
	% equalized signal
	detected = equalization_DFE(y_hat, c_opt, b, M1, M2, D);
	hard_d = zeros(length(detected),1);
	for k=1:length(hard_d)
		hard_d(k) = QPSK_detector(detected(k));
	end
	hard_bits = ibmap(hard_d);
	Pbit_DFE_uncode(i) = length(find(x(1:length(hard_bits))~=hard_bits))/length(hard_bits);
end
toc

figure();
semilogy(SNR_vect, Pbit_DFE_uncode, 'b', 'Marker', '^');
grid on;
ylim([10^-5 10^-1]); xlim([0 14]);
legend('Uncoded QPSK + DFE');

save('Pbit_DFE_uncoded.mat', 'Pbit_DFE_uncode');