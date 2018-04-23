function [pn] = PNSeq(L)
% Maximal length sequense of period L (pag. 233)

r = log2(L+1);
pn = zeros(L,1); 
pn(1:r) = ones(1,r).';    % Initial conditions

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
        case 8 
            pn(l) = xor(xor(pn(l-2),pn(l-3)),xor(pn(l-4),pn(l-8)));
        case 9
            pn(l) = xor(pn(l-5), pn(l-9));
        case 10
            pn(l) = xor(pn(l-7), pn(l-10));
            
            
            
    end
end

pn = 2*pn -1;      % pam modulation

end