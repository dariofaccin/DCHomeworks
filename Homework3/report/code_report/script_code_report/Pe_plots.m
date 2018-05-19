clc; close all; clear global; clearvars;
set(0,'defaultTextInterpreter','latex')    % latex format

SNR_vect = 8:14;
load('Pe_LE_avgs.mat');
load('PE_DFE_avgs.mat');
load('Pe_AA_GM_avgs.mat');
load('Pe_AA_NOGM_avgs.mat');
load('Pe_VA_avgs.mat');
load('Pe_FBA.mat');
load('Pe_AWGN_SIM_avgs.mat');


figure()
semilogy(SNR_vect, Pe_LE_avg, 'b--');
hold on; grid on;
semilogy(SNR_vect, Pe_DFE_avg,'b');
semilogy(SNR_vect, Pe_AA_GM_avg, 'k--');
semilogy(SNR_vect, Pe_AA_NOGM_avg, 'k');
semilogy(SNR_vect, Pe_VA_avg, 'r--');
semilogy(SNR_vect, Pe_FBA, 'r');
semilogy(SNR_vect, Pe_AWGN_SIM_avg, 'g--');
semilogy(SNR_vect, awgn_bound, 'g');
ylim([10^-4 10^-1]); xlim([8 14]);
xlabel('SNR'); ylabel('$P_e$');
legend('MF+LE@T','MF+DFE@T','AAF+MF+DFE@$\frac{T}{2}$','AAF+DFE@$\frac{T}{2}$',...
	'VA','FBA','MF b-S','MF b-T');
set(legend,'Interpreter','latex');