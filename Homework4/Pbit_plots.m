clc; close all; clear global; clearvars;
set(0,'defaultTextInterpreter','latex')    % latex format

SNR_vect = 0:14;
load('Pbit_DFE_avgs.mat');
load('Pbit_AWGN_avgs.mat');

semilogy(SNR_vect, Pbit_DFE_avg,'g');
hold on; grid on;
semilogy(SNR_vect, Pe_H_DFE_avg/2, 'g--');
semilogy(SNR_vect, Pbit_AWGN_avg, 'k');
semilogy(SNR_vect, Pe_H_AWGN_avg/2, 'k--');
% semilogy(SNR_vect, awgn_bound, 'g');
ylim([10^-5 10^-1]); xlim([0 12]);
xlabel('SNR'); ylabel('$P_{bit}$');
legend('Coded QPSK + DFE','Uncoded QPSK + DFE','Coded AWGN', 'Uncoded AWGN');
set(legend,'Interpreter','latex');