function [c_opt, Jmin] = WienerC_frac(h, r_w, sigma_a, M1, M2, D, N1, N2)

    padding = 60;
    hpad = padarray(h, padding);

    % Padding the noise correlation 
    r_w_pad = padarray(r_w, padding);   
    p  = zeros(M1 ,1);

    for i = 0 : M1-1
        p(i + 1) = sigma_a * conj(hpad(N1 + padding + 1 + 2*D - i));
    end

    R = zeros(M1);
    for row = 0:(M1-1)
        for col = 0:(M1-1)
            f=zeros(length(h),1);
            for n=0:length(h)-1
                f(n+1) = hpad(N1 + padding + 1 + 2 * n - col)*conj(hpad(N1 + padding + 1 + 2 * n - row));
            end
            fsum = sum(f);
            if M2==0
                ssum=0;
            else
                s=zeros(M2,1);
                for j=1:M2
                    s(j) = hpad(N1 + padding + 1 + 2*(j+D) -col)*conj(hpad(N1 + padding + 1+2*(j+D) -row ));

                end

                ssum = sum(s);

            end

            

            R(row + 1, col + 1) = sigma_a * (fsum - ssum) + r_w_pad(padding + 1 ...
                + row - col + (floor(length(r_w) / 2 )));

        end

    end

    

    c_opt = R \ p;



    temp2 = zeros(M1, 1);

    for l = 0:M1-1

        temp2(l + 1) = c_opt(l + 1) * hpad(N1 + padding + 1 +2*D-l); 

    end

    

    Jmin = 10*log10(sigma_a * (1 - sum(temp2)));

end