clc; close all; clear global; clearvars;
load('Input_symbols.mat');
set(0,'defaultTextInterpreter','latex')    % latex format

snr = 1;
snr_lin = 10^(snr/10);
M = 512;                % Block size
Npx = 5;
sigma_a = 2;

[r_c, sigma_w, g_srrc, g] = channel_OFDM(symbols_ak, snr, sigma_a);

a = r_c;

a_pad = [a; ones(M - mod(length(a), M), 1) * (-1-1i)];
a_matrix = reshape(a_pad, M, []); % it should mantain columnwise order

A_matrix = ifft(a_matrix);
A_matrix = [A_matrix(M-Npx+1:M, :); A_matrix]; % very powerful operationt
r = reshape(A_matrix, [], 1);


s_after_grrc = filter(g_srrc, 1, r);


t0 = 23;			% JUST A RANDOM NUMBER
g = g(1+t0 : end);
G = fft(g, 512);
G = G(:);

% r = r(1+t0 : end - mod(length(r), M+Npx) + t0);

r = r(1+t0:end);

r = downsample(r,4);
% perform the DFT
r = r(length(r)-2025*(M+Npx)+1:end);
r_matrix = reshape(r, M+Npx, []);
r_matrix = r_matrix(Npx + 1:end, :);

G_inv = G.^(-1);
x_matrix = fft(r_matrix);

y_matrix = bsxfun(@times, x_matrix, G_inv);

% Detect and compute BER
% sigma_i after the DFT and the scaling by G_i of each branch
sigma_i = 0.5*sigma_w*M*abs(G_inv).^2;
% Compute Log Likelihood Ratio
% It is different for each branch
llr_real = -2*bsxfun(@times, real(y_matrix), sigma_i.^(-1));
llr_imag = -2*bsxfun(@times, imag(y_matrix), sigma_i.^(-1));
llr_real_ar = reshape(llr_real, [], 1);
llr_imag_ar = reshape(llr_imag, [], 1);
llr = zeros(numel(llr_real) + numel(llr_imag), 1);
llr(1:2:end) = llr_real_ar;
llr(2:2:end) = llr_imag_ar;
% Drop the zero padding
% llr = llr(1:length(2*symbols_ak));
% Decode the bits
llr = deinterleaver(llr); % Deinterleave the loglikelihood ratio first
tic
dec_bits = LDPC_decoder(llr).';
toc

% %% PLOT
% % square-root raised cosine in T
% stem(-(length(g_srrc)-1)/2:(length(g_srrc)-1)/2,g_srrc)
% xlim([-(length(g_srrc)-1)/2 (length(g_srrc)-1)/2 ]), grid on
% 
% % square-root raised cosine in f
% [H, f] = freqz(g_srrc,1,1024,'whole');
% G_srrc = fftshift(H);
% % f = linspace(-(length(G_srrc)-1)/2,(length(G_srrc)-1)/2,length(f));
% figure, plot(f,abs(G_srrc)), grid on, xlim([0 2*pi])
% 
% % qc from homework 3
% figure, stem(0:length(qc)-1,qc)
% xlim([0 length(qc)-1]), grid on
% 
% % qr
% figure, stem(q_r), grid on