clc; close all; clear global; clearvars;
load('Input_symbols.mat');
set(0,'defaultTextInterpreter','latex')    % latex format

snr = 1;
snr_lin = 10^(snr/10);
M = 512;                % Block size
Npx = 16;
sigma_a = 2;

[r_c, sigma_w, g_srrc, g, t0] = channel_OFDM(symbols_ak, snr, sigma_a);
G = fft(g,512).';

% % a_pad = [a; ones(M - mod(length(a), M), 1) * (-1-1i)];
a_matrix = reshape(r_c(1:end-mod(length(r_c),M+Npx)), M+Npx, []); % it should mantain columnwise order
rn = a_matrix(Npx+1:end,:);
x_k = fft(rn);
K_i = 1./G;
y_matrix = x_k.*K_i;

% Detect and compute BER
% sigma_i after the DFT and the scaling by G_i of each branch
sigma_i = 0.5*sigma_w*M*abs(K_i).^2;
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
% figure, stem(g);
% f = linspace(0,1,length(G));
% plot(f,abs(G)), xlim([0 0.5]);
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