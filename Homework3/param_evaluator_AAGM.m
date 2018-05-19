clc; close all; clear global; clearvars;

load('Useful.mat');

sigma_a = 2;
sigma_w = sigma_a / 10;
M = 4;

Fpass = 0.2;
Fstop = 0.3;
Dpass = 0.057501127785;
Dstop = 0.01;
dens  = 20;
[N, Fo, Ao, W] = firpmord([Fpass, Fstop], [1 0], [Dpass, Dstop]);
g_AA  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(g_AA);

qg_up = conv(qc, g_AA);
qg_up = qg_up.';
t0_bar = find(qg_up == max(qg_up));

qg = downsample(qg_up(1:end), 2);
g_m = conj(flipud(qg));

h = conv(qg, g_m);
h = h(h ~= 0);
N0 = (sigma_a * 1) / (4 * 10);

r_g = xcorr(conv(g_AA, g_m));
r_w = N0 * downsample(r_g, 2);

N2 = floor(length(h)/2);
N1 = N2;

M1_span = 2:20;
D_span = 2:20;

Jvec = zeros(19);
for k=1:length(M1_span)
    for l=1:length(D_span)
        M1 = M1_span(k);
        D = D_span(l);
        M2 = N2 + M1 - 1 - D;
        [c, Jmin] = WienerC_frac(h, r_w, sigma_a, M1, M2, D, N1, N2);
        Jvec(k,l) = Jmin;
    end
end

figure, mesh(2:20, 2:20, reshape((Jvec(:, :)), size(Jvec(:, :), 2), size(Jvec(:, :), 2)))
title('Jmin for AA+GM');
xlabel('D'), ylabel('M1'), zlabel('Jmin [dB]')

[min, idx] = min(Jvec(:));

[idx_d, idx_m1] = ind2sub(size(Jvec), idx);