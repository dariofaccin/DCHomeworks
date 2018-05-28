clc; close all; clear global; clearvars;
load('Input_symbols.mat');
set(0,'defaultTextInterpreter','latex')    % latex format

snr = 2;
snr_lin = 10^(snr/10);
% Sub-channels
M = 512;   
% cycle prefix length
Npx = 8;
sigma_a = 2;
% OFDM symulationd
[r_c, sigma_w, g_srrc, g, t0, q_r] = channel_OFDM(symbols_ak, snr, sigma_a);
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
tic
dec_bits = LDPC_decoder(llr).';
toc

nerr = length(find(x(1:length(dec_bits))~=dec_bits))
Pbit = nerr/length(dec_bits)



% %% PLOT

% square-root raised cosine in T
stem(-(length(g_srrc)-1)/2:(length(g_srrc)-1)/2,g_srrc)
xlim([-(length(g_srrc)-1)/2 (length(g_srrc)-1)/2 ]), grid on
xlabel('$nT_c$')
ylabel('$g_{\sqrt{rcos}}(nT_c)$')

% square-root raised cosine in f
[G_srrc, f] = freqz(g_srrc,1,1024,'whole');
f = f/(2*pi);
figure, plot(f,10*log10(abs(G_srrc))), grid on, xlim([0 0.5])
ylim([-15 5])

% qc from homework 3
figure, stem(0:length(qc)-1,qc)
xlabel('$nT_c$')
ylabel('$q_c$')
xlim([0 length(qc)-1]), grid on

% qr
figure, stem(q_r), grid on
xlabel('$nT_c$')
ylabel('$q_r(nT_c)$')
xlim([0 length(q_r)-1]);

% Qr
[Q_r, f] = freqz(q_r,1,1024,'whole');
f = f/(2*pi);
plot(f,10*log10(abs(Q_r)))
grid on
xlim([0 0.5]);
ylim([-15 5]);



