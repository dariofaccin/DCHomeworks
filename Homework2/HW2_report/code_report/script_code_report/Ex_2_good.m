clc
clearvars
close all
set(0,'defaultTextInterpreter','latex')    % latex format

% Given Parameters
Tc = 1;
fd = (40*10^-5)/Tc;             % Doppler spread given
Tp = 1/10*(1/fd); 
N_h0 = 7500;                    % samples first plot
N_t = 80000;                    % samples of second plot

K_dB = 2;                       % Rice Factor in dB
K = 10^(K_dB/10);               % Rice Factor in linear units   
C = sqrt(K/(K+1));              

[a_ds, b_ds] = ClassicalDS();      % Parameters of the IIR filter which implement
                                   % the classical Doppler Spectrum (page 317) 
h_dopp = impz(b_ds, a_ds);         % Impulse response 
E_d = sum(h_dopp.^2);              % Energy of the impulse response
b_ds = b_ds/sqrt(E_d);             % Normalization
Md = 1-C^2;                        % normalization of the statistical power

% Doppler spectrum
[H_dopp,w]=freqz(h_dopp,1,1024,'whole',1/Tp);  
DS = abs(H_dopp).^2;

% figure                                  
subplot(121), plot(h_dopp,'r'), ylabel('$|h_{ds}|$'), hold on 
stem(1:length(h_dopp),real(h_dopp));
axis([0 Tp -0.15 0.25]), grid on
legend('continuous |h_{ds}|','sampled |h_{ds}|');
title('Impulse response of the IIR filter');
subplot(122), plot(w,10*log10(DS)), ylabel('$|D(f)|$'), grid on;
hold on, plot([fd fd], [-60 20], 'r--'), text(4.1e-4, 10, '$f_d$'); 
xlim([0 3*fd]), xlabel('f');
ylim([-60 20]);
legend('Doppler Spectrum');
title('Doppler Spectrum')

poles = abs(roots(a_ds));             % poles' magnitude 
most_imp = max(poles);  
tr = 5*Tp*ceil(-1/log(most_imp));     % transient as 5*Neq*Tp

h_samples_needed = N_t+tr;            % total length including the transient
w_samples_needed = ceil(h_samples_needed/Tp);
w = wgn(w_samples_needed,1,0,'complex');      % w ~ CN(0,1)
hprime_Tp = filter(b_ds, a_ds, w);

t = 1:length(hprime_Tp);                 % interpolation to Tq
Tq = Tc/Tp;
t_fine = Tq:Tq:length(hprime_Tp);
h_prime_Tq = interp1(t, hprime_Tp, t_fine, 'spline');
sigma = sqrt(Md);
h_prime_Tq = h_prime_Tq*sigma;            % impose the desired power delay profile
h0 = h_prime_Tq(tr+1:end)+C;              % remove the transient and add C 

figure, plot(abs(h0(1:N_h0)))             % first asked plot     
xlabel('$nT_c$')
ylabel('$|h_0(nT_c)|$')
xlim([1 N_h0]), grid on
title('Impulse response of the channel')

%% ESTIMATE OF THE PDF OF H_p=|h0|/sqrt(M)
h_p = h0/sqrt(C^2+Md);        
abs_h = abs(h_p);             % amplitude response of h_p
a=linspace(0,10,3000);
% Rice distribution
th_pdf = 2*(1+K).*a.*exp(-K-(1+K).*a.^2).*besseli(0,2.*a*sqrt(K*(1+K)));
% Estimate of the pdf
[y,t] = hist(abs_h,30);          
est_pdf = y/max(y);

figure                        % second asked plot
bar(t,est_pdf,'g'), hold on,plot(a,th_pdf,'r--','LineWidth',2); 
ylabel('Number of samples')
xlabel('Value')
title('Histogram of $\frac{|h_0|}{\sqrt{M_{|h_0|}}}$')
legend('estimate pdf','theoretical pdf');
axis([0 3 0 1.3]);
grid on

%% SPECTRUM ESTIMATION
% Theoretical PSD
Np = N_t*1/Tp;
[H_dopp,w] = freqz(h_dopp,1,Np,'whole');
H_dopp = (1/Np)*abs(H_dopp).^2;
DS = fftshift(H_dopp);
f2=[-Np/2+1:Np/2];

% Welch estimator
D = 40000;
S = D/2;                   % overlap
w_welch=window(@bartlett,D);
[Welch_P, N] = welchPSD(h0', w_welch, S);  
Welch_P = Welch_P/N;
Welch_cent=fftshift(Welch_P);
f1=[-N/2+1:N/2];

C_comp = 10*log10(C^2);             % Deterministic component
PSD_theo = 10*log10(Md*DS);
PSD_theo(length(PSD_theo)/2) = C_comp;

figure,
plot(f1, 10*log10(Welch_cent)) , hold on, plot(f2,PSD_theo,'r')
ylim([-40 0])
xlim([-5*N*fd  5*N*fd]);
xticks([-5*N*fd -4*N*fd -3*N*fd -2*N*fd -1*N*fd 0 1*N*fd 2*N*fd 3*N*fd 4*N*fd 5*N*fd])
xticklabels({'-5f_d','-4f_d','-3f_d','-2f_d','-f_d','0','f_d','2f_d','3f_d','4f_d','5f_d'});
ylabel('H(f) [dB]')
xlabel('f')
legend('Welch Periodogram','Theoretical PSD')
title('Spectrum Estimate')
grid on










