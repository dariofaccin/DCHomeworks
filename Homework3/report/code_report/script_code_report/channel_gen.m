clc; close all; clear global; clearvars;

% MUST BE USED ONCE TO GENERATE AN INPUT FOR
% ALL RECEIVERS

set(0,'defaultTextInterpreter','latex')    % latex format
T = 1;              % Symbol period
Tc = T/4;           % upsampling period
Q = T/Tc;           % Interpolation factor
snr_db = 10;
snr_lin = 10^(snr_db/10);
sigma_a = 2;		% Input variance

alpha = 0.67;		% Filter setup
beta = 0.7424;
qc_num = [0 0 0 0 0 beta];
qc_denom = [1 -alpha];
q_c = impz(qc_num, qc_denom);
q_c = [0; 0; 0; 0; 0; q_c( q_c >= max(q_c)*10^(-2) )]; 
E_qc = sum(q_c.^2);              % Energy of the impulse response

L = 1023;			% Input signal: from a PN sequence generate QPSK
x = PNSeq(L);
in_bits = bitmap(x(1:length(x)-1));

%% FIGURES
abs_qc = abs(qc);

figure()
stem(abs_qc);
xlabel('$m\frac{T}{4}$');
ylabel('$q_c$')
xlim([1 length(abs_qc)]);
grid on

[Qc, f] = freqz(q_c,1,1024,'whole');
% f = f/(2*pi);
f = linspace(0,4,length(f));

figure()
plot(f,10*log10(abs(Qc)));
xlim([0 2]);
xlabel('f'), grid on
ylabel('$| Q_c(f) |$ $[dB]$')

grid on