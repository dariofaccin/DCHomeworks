clc; close all; clear global; clearvars;
load('Input_symbols.mat');
load('qc')
set(0,'defaultTextInterpreter','latex')    % latex format

snr = 1;
M = 512;                % Block size
Npx = 7;
sigma_a = 2/M;


a = symbols_ak;
a_pad = [a; ones(M - mod(length(a), M), 1) * (-1-1i)];
a_matrix = reshape(a_pad, M, []); % it should mantain columnwise order

A_matrix = ifft(a_matrix);
A_matrix = [A_matrix(M-Npx+1:M, :); A_matrix]; % very powerful operationt
s = reshape(A_matrix, [], 1);

snr_lin = 10^(snr/10);

% square-root raised cosine
N    = 26;         % Order
Fc   = 0.46;       % Cutoff Frequency
TM   = 'Rolloff';  % Transition Mode
R    = 0.0625;     % Rolloff
DT   = 'sqrt';     % Design Type
Beta = 0.5;        % Window Parameter
win = kaiser(N+1, Beta);
g_srrc  = firrcos(N, Fc, R, 2, TM, DT, [], win);

q_r = conv(g_srrc,qc);
q_r = conv(q_r,g_srrc);

%% PLOT
% square-root raised cosine in T
stem(-(length(g_srrc)-1)/2:(length(g_srrc)-1)/2,g_srrc)
xlim([-(length(g_srrc)-1)/2 (length(g_srrc)-1)/2 ]), grid on

% square-root raised cosine in f
[H, f] = freqz(g_srrc,1,1024,'whole');
G_srrc = fftshift(H);
% f = linspace(-(length(G_srrc)-1)/2,(length(G_srrc)-1)/2,length(f));
figure, plot(f,abs(G_srrc)), grid on, xlim([0 2*pi])

% qc from homework 3
figure, stem(0:length(qc)-1,qc)
xlim([0 length(qc)-1]), grid on

% qr
figure, stem(q_r), grid on