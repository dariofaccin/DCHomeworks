clc; close all; clear global; clearvars;

load('Useful.mat');

sigma_a = 2;
sigma_w = sigma_a / 10;
M = 4;

gm = conj(qc(end:-1:1));
h = conv(qc,gm);
h = h(h>max(h)/100);
h = h(3:end-2);
h_T = downsample(h,4);
t0_bar = length(gm);
r_gm = xcorr(gm,gm);
rw_tilde = sigma_w/4 .* downsample(r_gm, 4);
N2 = floor(length(h_T)/2);
N1 = N2;

M1_span = 2:20;
D_span = 2:20;

Jvec = zeros(19);
for k=1:length(M1_span)
    for l=1:length(D_span)
        M1 = M1_span(k);
        D = D_span(l);
        M2 = N2 + M1 - 1 - D;
        [c, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);
        Jvec(k,l) = Jmin;
    end
end

figure, mesh(2:20, 2:20, reshape((Jvec(:, :)), size(Jvec(:, :), 2), size(Jvec(:, :), 2)))
title('Jmin for DFE, SNR = 10 (dB)'); view(160,20);
xlim([2 20]); ylim([2 20]);
xlabel('D'), ylabel('M1'), zlabel('Jmin (dB)')

[min, idx] = min(Jvec(:));

[idx_d, idx_m1] = ind2sub(size(Jvec), idx);