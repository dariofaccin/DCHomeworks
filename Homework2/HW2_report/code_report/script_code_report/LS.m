function [h_ls]=LS(x,d,L)
% Compute the solution of the ls problem (pag. 246)
% x the input seq
% d filter output
% L half length PN seq
% h_ls the least squares estimate of h

I = zeros(L);
for k=1:L
    I(:,k)=x(L-k+1:(2*L-k));
end
o = d(L:2*L-1);
Phi = I'*I;
theta = I'*o;

h_ls = inv(Phi)*theta;

end