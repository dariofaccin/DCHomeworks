clc; close all; clear global; clearvars;

L=127;               % length of PN sequence
Nh=20;             % Bound on the length of h
% Additive noise
sigdB = -8;
sigmaw = 10^(sigdB/10);
w = wgn(4*L,1,sigdB);
w_0 = w([1:2:end-1]);
w_1 = w([2:2:end]);

PN = PNSeq(L);      % ML sequence repeated once
x=[PN ; PN];

%% POLYPHASE REALIZATION
a1 = -0.9635;
a2 = 0.4642;
h = impz(1, [1 a1 a2]);
h = h(1:Nh);
[h_even,h_odd] = polyphase(h,Nh);

%% ESTIMATE OF h0 WITH THE CORRELATION METHOD
z_0=filter(h_even, 1, x);       
d_0 = z_0 + w_0;
h0_cor=corr_method(x, d_0);     % correlation method
var_h0_cor = sigmaw/(L/2);          % variance of the estimate  

%% ESTIMATE OF h1 WITH THE CORRELATION METHOD
z_1=filter(h_odd, 1, x);
d_1 = z_1 + w_1;
h1_cor=corr_method(x, d_1);         % correlation method
var_h1_cor = sigmaw/(L/2);          % variance of the estimate

% Estimated impulse response
h_cor=zeros(Nh,1);
h0_spaced=zeros(Nh,1);
h1_spaced=zeros(Nh,1);
for i=1:Nh
    if (i<=L)
    h_cor(2*i-1)=h0_cor(i);
    h_cor(2*i)=h1_cor(i);
    h0_spaced(2*i-1) = h0_cor(i);
    h1_spaced(2*i) = h1_cor(i);
    end
end

h_cor=h_cor(1:Nh);
var_cor = var_h0_cor + var_h1_cor;
h0_spaced=h0_spaced(1:Nh);
h1_spaced=h1_spaced(1:Nh);

% Plot 
% figure, 
% stem(0:Nh-1,h0_spaced,'ro'), hold on
% stem(0:Nh-1,h1_spaced,'ko'), hold on
% plot([0:Nh-1],h,'bx','LineWidth',1.5);
% legend('h_0 component','h_1 component','Analytic impulse response');
% xlabel('n'), ylim([-0.4 1.2]);
% title('Polyphase components with the Correlation method')
% grid on

%% ESTIMATE WITH THE LS METHOD
h0_ls=LS(x, d_0, L);
h1_ls=LS(x, d_1, L);

h_ls=zeros(Nh,1);
h0_spaced_ls=zeros(Nh,1);
h1_spaced_ls=zeros(Nh,1);
for i=1:Nh
    if (i<=L)
    h_ls(2*i-1)=h0_ls(i);
    h_ls(2*i)=h1_ls(i);
    h0_spaced_ls(2*i-1) = h0_ls(i);     
    h1_spaced_ls(2*i) = h1_ls(i);
    end
end
h_ls=h_ls(1:Nh);
h0_spaced_ls=h0_spaced_ls(1:Nh);
h1_spaced_ls=h1_spaced_ls(1:Nh);

% Plot 
% figure, 
% stem(0:Nh-1,h0_spaced_ls,'ro'), hold on
% stem(0:Nh-1,h1_spaced_ls,'ko'), hold on
% plot([0:Nh-1],h,'bx','LineWidth',1.5);
% legend('h_0 component','h_1 component','Analytic impulse response');
% xlabel('n'), ylim([-0.4 1.2]);
% title('Polyphase components with the LS method')
% grid on

%% comparison cor vs ls
% figure, 
% stem(0:Nh-1,h_cor,'bo')
% hold on
% stem(0:Nh-1,h_ls,'ro')
% grid on
% plot([0:Nh-1],h,'bx','LineWidth',1.5);
% legend('h_cor','h_ls','Analytic impulse response');

%% COST FUNCTION
d0_hat = filter(h0_cor,1,x);
d1_hat = filter(h1_cor,1,x);
d_hat_cor = PS(d0_hat, d1_hat); 
error_cor = d - d_hat_cor;
E_cor = sum(error_cor(L:2*L-1).^2);
E_L_cor = 10*log10(E_cor/L);

d0_hat = filter(h0_ls,1,x);
d1_hat = filter(h1_ls,1,x);
d_hat_ls = PS(d0_hat, d1_hat);
error_ls = d - d_hat_ls;
E_ls = sum(error_ls(L:2*L).^2);
E_L_ls(1,n) = 10*log10(E_ls/L);



















