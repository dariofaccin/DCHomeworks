clc; close all; clear global; clearvars;

load('Useful.mat', 'qc');
load('Input_symbols.mat');

SNR_vect = 1.4:0.05:1.8;
sigma_a = 2;				% Input variance
M = 4;					% Constellation cardinality

gm = conj(qc(end:-1:1));		% Matched filter: complex conjugate of qc
h = conv(qc,gm);			% Impulse response

t0_bar = find(h == max(h));		% Timing phase: peak of h
h = h(h>max(h)/100);
h = h(3:end-2);
h_T = downsample(h,4);

r_gm = xcorr(gm,gm);			% Matched filter autocorrelation
Pbit_DFE_code = zeros(length(SNR_vect),1);
N2 = floor(length(h_T)/2);
N1 = N2;
M1 = 5;
D = 4;
M2 = N2 + M1 - 1 - D;
tic
parfor i=1:length(SNR_vect)
	snr_db = SNR_vect(i);
	snr_lin = 10^(snr_db/10);
	% Single carrier channel simulation
	[r_c, sigma_w, qc] = channel_sim(symbols_ak, snr_db, sigma_a);
	% Additive complex-Gaussian noise
	w = wgn(length(r_c),1, 10*log10(sigma_w), 'complex');
	r_c = r_c + w;
	r_c_prime = filter(gm,1,r_c);
	r_c_prime = r_c_prime(t0_bar:end);
	% Signal at the input of the DFE
	x_aa = downsample(r_c_prime,4);
	rw_tilde = sigma_w/4 .* downsample(r_gm, 4);
	% DFE
	[c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);
	psi = conv(c_opt, h_T);
	b = - psi(end - M2 + 1:end);
	y_hat = x_aa/max(psi);
	% Equalized signal
	detected = equalization_DFE(y_hat, c_opt, b, M1, M2, D);
	% Decoder as LLR
	llr = zeros(2*length(detected),1);
	Jmin_lin = 10^(Jmin/10);
	noise_var = (Jmin_lin-sigma_a*abs(1-max(psi))^2)/(abs(max(psi))^2);
	llr(1:2:end) = -2*real(detected)/(noise_var/2);
	llr(2:2:end) = -2*imag(detected)/(noise_var/2);
	% Deinterleave and decode bits
	llr_deinter = deinterleaver(llr);
	decoded = LDPC_decoder(llr_deinter).';
	Pbit_DFE_code(i) = length(find(x(1:length(decoded))~=decoded))/length(decoded);
end
toc