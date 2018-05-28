clc; close all; clear global; clearvars;
set(0,'defaultTextInterpreter','latex')    % latex format

%% Uncoded OFDM + DFE + AWGN
SNR_vect_uncoded = 0:14;
load('Pbit_DFE_uncoded.mat');
load('Pbit_OFDM_uncoded.mat');
load('Pbit_AWGN_uncoded.mat');
figure();
semilogy(SNR_vect_uncoded, Pbit_DFE_uncode,'g');
hold on; grid on;
semilogy(SNR_vect_uncoded, Pbit_OFDM_uncode, 'b');
semilogy(SNR_vect_uncoded, Pbit_AWGN_uncode, 'k');
ylim([10^-5 10^-1]); xlim([4 14]);
xlabel('SNR'); ylabel('$P_{bit}$');
legend('Uncoded DFE','Uncoded OFDM','Uncoded AWGN');
set(legend,'Interpreter','latex');

%% Coded OFDM + DFE + AWGN
SNR_vect_coded = 2:0.1:3;
figure();
load('Pbit_DFE_coded.mat');
load('Pbit_OFDM_coded.mat');
load('Pbit_AWGN_coded.mat');
semilogy(SNR_vect_coded, Pbit_DFE_code,'g');
hold on; grid on;
semilogy(SNR_vect_coded, Pbit_OFDM_code, 'b');
semilogy(SNR_vect_coded, Pbit_AWGN_code, 'k');
ylim([10^-5 10^-1]); xlim([2 3]);
xlabel('SNR'); ylabel('$P_{bit}$');
legend('Coded DFE','Coded OFDM', 'Coded AWGN');
set(legend,'Interpreter','latex');