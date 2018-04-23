clc; close all; clear global; clearvars;

N = [1:1:20];
L = [31 63 127 255 511 1023];
sigdB = -8;
sigmaw = 10^(sigdB/10);
a1 = -0.9635;
a2 = 0.4642;
noise = wgn(4*max(L),1,sigdB);
% save('good_noise3','noise')
load good_noise3

%%
index=0;
for l=1:length(L)
    index = index+1;
for n=1:length(N)
    w = noise(1:4*L(l));
    w_0 = w([1:2:end]);
    w_1 = w([2:2:end]);
    PN = PNSeq(L(l));         % ML sequence repeated once
    x=[PN ; PN];
    
    h = impz(1, [1 a1 a2]);   % Analytical h
   [h_even,h_odd] = polyphase(h,length(h));
    
    % scheme pag 239
    z_0=filter(h_even, 1, x);
    z_1=filter(h_odd, 1, x);
    d_0 = z_0 + w_0;
    d_1 = z_1 + w_1;
    d = PS(d_0, d_1);
    
    % Correlation method
    h0_cor=corr_method(x, d_0);
    h1_cor=corr_method(x, d_1);
    if N(n)<L(l)
        h0_cor=h0_cor(1:ceil(N(n)/2));
        h1_cor=h1_cor(1:floor(N(n)/2));
    end
    h_cor = PS(h0_cor, h1_cor);
        
    % Cost function
    d0_hat = filter(h0_cor,1,x);
    d1_hat = filter(h1_cor,1,x);
    d_hat_cor = PS(d0_hat, d1_hat); 
    error_cor = d - d_hat_cor;
    E_cor = sum(error_cor(L(l):2*L(l)).^2);
    E_L_cor(1,n) = 10*log10(E_cor/L(l));
    
    % ls method
    h0_ls=LS(x, d_0, L(l));
    h1_ls=LS(x, d_1, L(l));
    if N(n)<L(l)
        h0_ls=h0_ls(1:ceil(N(n)/2));
        h1_ls=h1_ls(1:floor(N(n)/2));
    end
    h_ls = PS(h0_ls, h1_ls);
    
    % Cost Function
    d0_hat = filter(h0_ls,1,x);
    d1_hat = filter(h1_ls,1,x);
    d_hat_ls = PS(d0_hat, d1_hat);
    error_ls = d - d_hat_ls;
    E_ls = sum(error_ls(L(l):2*L(l)).^2);
    E_L_ls(1,n) = 10*log10(E_ls/L(l));
    
end
Cost_cor(:,index) = E_L_cor;
Cost_ls(:,index) = E_L_ls;

end

plot_est(Cost_cor,Cost_ls,sigdB);

% figure
% grid on
% plot(N,Cost_cor,'o-')
% hold on
% plot(N, Cost_ls,'--')
% a = sigdB*ones(1,max(L));
% hold on
% plot(a,'b--','LineWidth',2)
% xlim([1 20]);  



    
% figure
% grid on
% plot(N,Cost_ls,'-')
% hold on
% a = sigdB*ones(1,20);
% plot(a,'b--')
% xlim([1 20]);
% title('ls')
    
    
    
    
    
    
    
    
    
    
    
    
    
    