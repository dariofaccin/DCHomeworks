clc; close all; clear global; clearvars;

%% SETUP OF THE GIVEN PARAMETERS
Tc=1;
fd=(40*10^-5)/Tc;
Tp=1/10*(1/fd); 
Nsamples=80000;
hplot_samples=7500;

KdB=2;
K=10^(KdB/10); 
C=sqrt(K/(K+1)); 

%setup of the doppler filter coefficients (page 317)
a_ds=[1, -4.4153, 8.6283, -9.4592, 6.1051, -1.3542, -3.3622, 7.2390, -7.9361, 5.1221, -1.8401, 2.8706e-1];
b_ds=[1.3651e-4, 8.1905e-4, 2.0476e-3, 2.7302e-3, 2.0476e-3, 9.0939e-4, 6.7852e-4, 1.3550e-3, 1.8076e-3, 1.3550e-3, 5.3726e-4, 6.1818e-5, -7.1294e-5, ...
    -9.5058e-5, -7.1294e-5, -2.5505e-5, 1.3321e-5, 4.5186e-5, 6.0248e-5, 4.5186e-5, 1.8074e-5, 3.0124e-6];
% The energy of the doppler filter need to be normalized to 1
h_ds=impz(b_ds, a_ds);
E_hds=sum(h_ds.^2);
b_ds=b_ds/sqrt(E_hds);
%normalize the power of h0
Md=1-C^2;

%% SIMULATION OF THE CHANNEL IMPULSE RESPONSE
%compute the transient as 5*Neq*Tp because the final sampling time is Tq 
pvec=abs(roots(a_ds));
pmax=max(pvec);
Neq=ceil(-1/log(pmax));
transient=5*Tp*Neq;

h_nsamples=Nsamples+transient;
w_samples=ceil(h_nsamples/Tp);
% Complex-valued Gaussian white noise with zero mean and unit variance
w=wgn(w_samples,1,0,'complex');

h_full=filter(b_ds, a_ds, w);
figure, plot(impz(b_ds,a_ds)), ylabel('|h_{ds}|')

%interpolation to Tq
t = 1:length(h_full);
t_int = Tc/Tp:Tc/Tp:length(h_full);
h_int = interp1(t, h_full, t_int, 'spline');

%multiply by sqrt(M_h0) to give the desired power
h_int=h_int*sqrt(Md);

%drop the transient and add C because of LOS component
h=h_int(transient+1:end)+C;


%plot of 7500 samples of h
figure, 
plot(abs(h(1:hplot_samples)))
xlabel('nT_c')
ylabel('|h_0(nT_c)|')
xlim([1 hplot_samples])
title('Impulse response of the channel')

%% ESTIMATE OF THE PDF OF H0bar

h_bar=h/sqrt(C^2+Md); % the normalization here does not make much sense
% as M_h0=1-C^2, but it's to keep the formulas as in the book

% mean and variance of the realization of h
h_sd=std(abs(h_bar));
h_mean=mean(abs(h_bar));
% magnitude
mag_h=abs(h_bar);
% fit the data with a Rice ditribution to derive the parameters 
[v, s]=ricefit(mag_h);
% parameters v ans s are related to the rice factor K (wikipedia rician
% distribution)
K_est=v^2/(2*s^2); % estimated value of K 

% plot an histogram of h to get an idea of the shape of the distribution
Nbins=100;
figure,
hist(mag_h, Nbins)
ylabel('Number of samples')
xlabel('Value')
title('Histogram of h_0')

% compute the theoretical distribution
x=linspace(0,3,1000);
v_th=sqrt(K/(K+1));
s_th=sqrt(1/(2*(K+1)));
th=ricepdf(x,v_th,s_th );


%% SPECTRUM ESTIMATION

% Welch estimator
S=600;   %overlap
D=1200;   %window length

w_welch=window(@hamming,D);
%w_welch=kaiser(D,5);

[Welch_P, Ns] = welchPSD(h', w_welch, S);

f=1/Tc:1/Tc:Nsamples;
Welch_centered=fftshift(Welch_P);
figure,
plot(10*log10(Welch_centered))
%ylim([-10 30])
%xlim([Nsamples/2-5*Nsamples*fd Nsamples/2+5*Nsamples*fd])
xticks([39850 39900 39950 40000 40050 40100 40150])
xticklabels({'-150','-100','-50','0','50','100','150'});
ylabel('H(f) [dB]')
xlabel('f')
