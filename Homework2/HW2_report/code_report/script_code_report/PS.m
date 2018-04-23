function [h] = PS(h0, h1)
% Compute Parallel-to-series from the two polyphase components

temp = length(h0)+length(h1);
h=zeros(temp,1);
if mod(temp,2)==0
    for i=1:temp/2
        h(2*i-1)=h0(i);
        h(2*i)=h1(i);
    end
elseif  length(h0)==1
        h(1)=h0(1);
else
    for i=1:length(h0)-1
        h(2*i-1)=h0(i);
        h(2*i)=h1(i);
    end
    h(2*(i+1)-1) = h0(i+1);
end

end