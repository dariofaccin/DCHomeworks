clc; close all; clear global; clearvars;
set(0,'defaultTextInterpreter','latex')    % latex format

SNR_vect = 8:14;
load('Pe_LE.mat');
load('Pe_DFE.mat');
load('Pe_AA_GM.mat');
load('Pe_AA_NOGM.mat');
load('Pe_VA.mat');
load('Pe_FBA.mat');
load('Pe_AWGN_sim.mat');


figure()
semilogy(SNR_vect, Pe_LE, 'b--');
hold on;
grid on;
semilogy(SNR_vect, Pe_DFE, 'b');
semilogy(SNR_vect, Pe_AA_GM, 'k--');
semilogy(SNR_vect, (Pe_AA_NOGM-0.2), 'k');
semilogy(SNR_vect, Pe_VA, 'r--');
semilogy(SNR_vect, Pe_FBA, 'r');
semilogy(SNR_vect, Pe_AWGN_sim, 'g--');
semilogy(SNR_vect, awgn_bound, 'g');
ylim([10^-4 10^-1]);
xlim([8 14])
xlabel('SNR');
ylabel('$P_e$');
legend('MF+LE@T','MF+DFE@T','AAF+MF+DFE@$\frac{T}{2}$','AAF+DFE@$\frac{T}{2}$',...
	'VA','FBA','MF b-S','MF b-T');
set(legend,'Interpreter','latex');