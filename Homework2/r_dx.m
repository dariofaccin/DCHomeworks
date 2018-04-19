function [rx]=r_dx(x,d)
% INPUT
% x the input sequence of length 2*L
% r the output of the filter
% OUTPUT
% rx the cross correlation of d and x of length L=length(x)/2

L=length(x)/2;
rx=zeros(L, 1);
for m=1:L-1
    rtemp=zeros(L,1);
    for k=1:L
        %starts using the samples of d after a transient of ength L-1
        rtemp(k)=d(L-2+k)*conj(x(L-1+k-m));
    end
    rx(m)=sum(rtemp)/L;
end
end