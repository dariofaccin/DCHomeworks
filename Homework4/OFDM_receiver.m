clc; close all; clear global; clearvars;
load('Input_symbols.mat');
set(0,'defaultTextInterpreter','latex')    % latex format

snr = 2.5;
snr_lin = 10^(snr/10);
% Sub-channels
M = 512;   
% cycle prefix length
sigma_a = 2;

Npx = 11;
% OFDM symulation
[r_c, sigma_w, g_srrc, g, t0, q_r] = channel_OFDM(symbols_ak, snr, sigma_a, Npx);
G = fft(g,512).';
a_matrix = reshape(r_c(1:end-mod(length(r_c),M+Npx)), M+Npx, []);
rn = a_matrix(Npx+1:end,:);
x_k = fft(rn);
K_i = 1./G;
y_matrix = x_k.*K_i;

% Detect and compute BER
% sigma_i after the DFT and the scaling by G_i of each branch
sigma_i = 0.5*sigma_w*M*abs(K_i).^2;
% Compute Log Likelihood Ratio
% It is different for each branch
llr_real = -2*real(y_matrix).*sigma_i.^(-1);
llr_imag = -2*imag(y_matrix).*sigma_i.^(-1);
llr_real_ar = reshape(llr_real, [], 1);
llr_imag_ar = reshape(llr_imag, [], 1);
llr = zeros(numel(llr_real) + numel(llr_imag), 1);
llr(1:2:end) = llr_real_ar;
llr(2:2:end) = llr_imag_ar;
% Drop the zero padding
% llr = llr(1:2*length(x));
% Decode the bits
llr = deinterleaver(llr); % Deinterleave the loglikelihood ratio first
dec_bits = LDPC_decoder(llr).';
nerr = length(find(x(1:length(dec_bits))~=dec_bits));
Pbit = nerr/length(dec_bits);