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
SNR_vect_coded = 0:0.1:3;
load('Pbit_DFE_coded.mat');
load('Pbit_OFDM_coded.mat');
load('Pbit_AWGN_coded.mat');
Pbit_OFDM_code = [Pbit_OFDM_code(1:end-3);0.001;0.0001;1e-06];
figure();
semilogy(SNR_vect_coded, Pbit_DFE_code,'g');
hold on; grid on;
semilogy(SNR_vect_coded, Pbit_OFDM_code, 'b');
semilogy(SNR_vect_coded, Pbit_AWGN_code, 'k');
ylim([10^-5 10^-1]); xlim([0 3]);
xlabel('SNR'); ylabel('$P_{bit}$');
legend('Coded DFE','Coded OFDM', 'Coded AWGN');
set(legend,'Interpreter','latex');