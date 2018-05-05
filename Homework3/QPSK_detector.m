function [outsym] = QPSK_detector(insym)

if (real(insym)>0)
    if (imag(insym)>0)
        outsym = 1+1i;
    else
        outsym = 1-1i;
    end
else
    if (iamg(insym)>0)
        outsym = -1+1i;
    else
        outsym = -1-1i;
    end
end

end