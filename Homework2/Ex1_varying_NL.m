clc; close all; clear global; clearvars;

N = [1:1:20];
L = [3 7 15 31 63 127];

sigdB = -8;
sigmaw = 10^(sigdB/10);
a1 = -0.9635;
a2 = 0.4642;
index = 0;
%%
for l=1:length(L)
    index = index+1;
for n=1:length(N)
    w = wgn(4*L(l),1,sigdB);
    w_0 = w([1:2:end-1]);
    w_1 = w([2:2:end]);
    PN = PNSeq(L(l));         % ML sequence repeated once
    x=[PN ; PN];
    
    h = impz(1, [1 a1 a2]);
    h = h(1:N(n));            % Analytical h
    x_up = upsample(x,2);
    d = filter(1,[1 a1 a2],x_up)+w;            % ideal d(k)
    
    [h_even,h_odd] = polyphase(h,N(n));
 
    z_0=filter(h_even, 1, x);
    z_1=filter(h_odd, 1, x);
    d_0 = z_0 + w_0;
    d_1 = z_1 + w_1;
    
    % Correlation method
    h0_cor=corr_method(x, d_0);
    h1_cor=corr_method(x, d_1);
    if N(n)<L(l)
        h0_cor=h0_cor(1:ceil(N(n)/2));
        h1_cor=h1_cor(1:floor(N(n)/2));
    end
    temp = length(h0_cor)+length(h1_cor);
    h_cor=zeros(temp,1);
    if mod(temp,2)==0
        for i=1:temp/2
            h_cor(2*i-1)=h0_cor(i);
            h_cor(2*i)=h1_cor(i);
        end
    elseif  i==1
        h_cor(i)=h0_cor(i);
    else
        for i=1:length(h0_cor)-1
            h_cor(2*i-1)=h0_cor(i);
            h_cor(2*i)=h1_cor(i);
        end
        h_cor(2*(i+1)-1) = h0_cor(i+1);
    end
    % Cost function
    d_hat_cor = filter(h_cor,1,x_up)+w;
    error_cor = d - d_hat_cor;
    E_cor = sum(abs(error_cor(L(l)-1:2*L(l)-2)).^2);
    E_L_cor(1,n) = 10*log10(E_cor/L(l));
    
    % ls method
    h0_ls=LS(x, d_0, L(l));
    h1_ls=LS(x, d_1, L(l));
    h_ls=zeros(N(n),1);
    for i=1:N(n)
        if (i<=L)
        h_ls(2*i-1)=h0_ls(i);
        h_ls(2*i)=h1_ls(i);
        end
    end
    h_ls=h_ls(1:N(n));  
    d_hat_ls = filter(h_ls,1,x_up)+w;
    error_ls = d - d_hat_ls;
    E_ls = sum(abs(error_ls(L(l)-1:2*L(l)-2)).^2);
    E_L_ls(1,n) = 10*log10(E_ls/L(l));
    
    figure
    stem(h,'r*'), hold on
    stem(h_cor,'bo'), hold on
    stem(h_ls,'go'), hold on
    
end
Cost_cor(:,index) = E_L_cor;
Cost_ls(:,index) = E_L_ls;

end

%%
figure
grid on
plot(N,Cost_cor,'--')
hold on
a = sigdB*ones(1,20);
plot(a,'b--')
xlim([1 20]);
title('Cor')
    
figure
grid on
plot(N,Cost_ls,'-')
hold on
a = sigdB*ones(1,20);
plot(a,'b--')
xlim([1 20]);
title('ls')
    
    
    
    
    
    
    
    
    
    
    
    
    
    