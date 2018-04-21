function [h_even h_odd] = polyphase(h,Nlim)

% Even samples
h_even=zeros(ceil(Nlim/2),1);
for k=1:(Nlim/2)
    h_even(k)=h(2*k-1);
end

% Odd samples
h_odd=zeros(floor(Nlim/2),1);
for k=1:Nlim/2
    h_odd(k)=h(2*k);
end

end