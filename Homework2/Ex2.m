close all; clear global; clearvars; clc;

Tc = 1;

% Rice factor
K_db = 2;
K = 10^(K_db/10);
C=sqrt(K/(K+1));

% Doppler spread
fd = 40 * 10^(-5);
Tp = 0.1*(1/fd);

% Denominator and Numerator (pg 317)
a_dopp=[1, -4.4153, 8.6283, -9.4592, 6.1051, -1.3542, -3.3622, 7.2390, -7.9361, 5.1221, -1.8401, 2.8706e-1];
b_dopp=[1.3651e-4, 8.1905e-4, 2.0476e-3, 2.7302e-3, 2.0476e-3, 9.0939e-4, 6.7852e-4, 1.3550e-3, 1.8076e-3, 1.3550e-3, 5.3726e-4, 6.1818e-5, -7.1294e-5, ...
    -9.5058e-5, -7.1294e-5, -2.5505e-5, 1.3321e-5, 4.5186e-5, 6.0248e-5, 4.5186e-5, 1.8074e-5, 3.0124e-6];

[h_dopp] = impz(b_dopp, a_dopp);

% Energy normalization
hds_nrg = sum(h_dopp.^2);
b_dopp = b_dopp/sqrt(hds_nrg);
Md=1-C^2;

% Number of coefficients (given)
N_h = 1;
num_samples = 7500;

% Transient is determined by the pole closest to the unit circle
poles = abs(roots(a_dopp));
most_imp = max(poles);
transient = 5*Tp*ceil(-1/log(most_imp));

% Required samples
h_samples_needed = num_samples + transient;
w_samples_needed = ceil(h_samples_needed / Tp);

% h_mat = zeros(N_h, h_samples_needed - transient);

% Noise: Complex Gaussian, zero mean, unit variance
w = wgn(w_samples_needed,1,0,'complex');

hprime = filter(b_dopp, a_dopp, w);

% Interpolation
t = 1:length(hprime);
t_fine = Tc/Tp:Tc/Tp:length(hprime);
h_fine = interp1(t, hprime, t_fine, 'spline');
h_fine = h_fine*sqrt(Md);

% Drop the transient and add C because of LOS components
h = h_fine(transient+1:end)+C;

% Plot 7500 samples
figure('Name', 'Coefficient h0');
plot(abs(h));
xlabel('nT_c'), ylabel('|h_0(nT_c)|');
title('|h_0(nT_C)|');
xlim([0 num_samples]);

%% PDF of h0_bar

h_bar=h/sqrt(C^2+Md);

% mean and variance of the realization of h
h_sd=std(abs(h_bar));
h_mean=mean(abs(h_bar));

% Magnitude
abs_h=abs(h_bar);

% Histogram to see the distribution behavior
Nbins=100;
figure('Name','Distribution');
hist(abs_h, Nbins);
ylabel('Number of samples');
xlabel('Value');
title('Distribution h_0');