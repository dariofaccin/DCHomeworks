function [h_ls]=LS(x,d,L)
% INPUT
% x the input seq
% d filter output
% L half length PN seq
% OUTPUT
% h_ls the least squares estimate of h

%build matrix Phi and theta by I and o (page 246)
I=zeros(L);
for k=1:L
    I(:,k)=x(L-k+1:(2*L-k));
end
o=d(L:2*L-1);
Phi=I'*I;
theta=I'*o;

h_ls=inv(Phi)*theta;
%h_ls=h_ls(1:N);
end