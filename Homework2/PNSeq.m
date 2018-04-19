function [pn] = PNSeq(L)

r = log2(L+1);
pn = zeros(L,1);

% Initial conditions (set to one, arbitrary)
% Must not be ALL zeros
pn(1:r) = ones(1,r).';

for l=r+1:L
    switch r
        case 1
            pn(l) = pn(l-1);
        case 2
            pn(l) = xor(pn(l-1), pn(l-2));
        case 3
            pn(l) = xor(pn(l-2), pn(l-3));
        case 4
            pn(l) = xor(pn(l-3), pn(l-4));
        case 5
            pn(l) = xor(pn(l-3), pn(l-5));
        case 6
            pn(l) = xor(pn(l-5), pn(l-6));
        case 7
            pn(l) = xor(pn(l-6), pn(l-7));
    end
end

% Bits are {-1, 1}
for i=1:L
	if pn(i)==0
        pn(i)=-1;
	end
end

end