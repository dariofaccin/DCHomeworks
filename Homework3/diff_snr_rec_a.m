clc; close all; clear global; clearvars;

% Input
load('Useful.mat');
SNR_vect = 8:14;
Pe_LE = zeros(length(SNR_vect),1);
errors = zeros(length(SNR_vect),1);
awgn_bound = zeros(length(SNR_vect),1);

sigma_a = 2;
M = 4;

for i=1:length(SNR_vect)
    snr_db = SNR_vect(i);
    snr_lin = 10^(snr_db/10);
    [r_c, sigma_w, qc] = channel_sim(in_bits, snr_db, sigma_a);
    r_c = r_c + w(:,i);
    gm = conj(qc(end:-1:1));
    h = conv(qc,gm);
    h = h(h>max(h)/100);
    h = h(3:end-2);
    h_T = downsample(h,4);

    % Filtering received signal
    r_c_prime = filter(gm,1,r_c);
    
    t0_bar = length(gm);
    r_c_prime = r_c_prime(t0_bar:end);
    x = downsample(r_c_prime,4);

    r_gm = xcorr(gm,gm);
    rw_tilde = sigma_w/4 .* downsample(r_gm, 4);

    M1 = 5;
    M2 = 0;
    D = 2;
    [c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);

    detected = equalization_LE(x, c_opt, M1, D, max(conv(c_opt, h_T)));

    [Pe_LE(i), errors(i)] = SER(in_bits(1:length(detected)), detected);
    awgn_bound(i) = 4*(1-1/sqrt(M))*qfunc(sqrt(snr_lin/(sigma_a/2)));
end

save('Pe_LE.mat','Pe_LE', 'awgn_bound');