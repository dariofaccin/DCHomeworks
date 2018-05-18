function [detected] = FBA(y, psiD, L1, L2)
% y: the data after filter c (hopefully no effect of the precursors)
% psiD: overall system impulse response after the cancellation 
% of precursors by filter c
% L1: number of considered precursors for each symbol
% L2: number of considered precursors for each symbol

M = 4;              % cardinality of the constellation
Ns = M^(L1+L2);     % # of states
K = length(y);      
%QPSK symbols
symb = [1+1i, 1-1i, -1+1i, -1-1i]; 

tStart = tic;   

% initialize matrix with the input data U
states_symbols = zeros(Ns, M);
statelength = L1 + L2; 
statevec = zeros(1, statelength);
U = zeros(Ns, M);
for state = 1:Ns
    for j = 1:M
        lastsymbols = [symb(statevec + 1), symb(j)]; 
        U(state, j) = lastsymbols * flipud(psiD);
    end    
    states_symbols(state,:) = lastsymbols(1:M); 
    % update statevec
    statevec(statelength) = statevec(statelength) + 1;
    i = statelength;
    while (statevec(i) >= M && i > 1)
        statevec(i) = 0;
        i = i-1;
        statevec(i) = statevec(i) + 1;
    end
end

%computation of matrix C (3D)
c = zeros(M, Ns, K+1);

for k = 1:K
    c(:, :, k) = (-abs(y(k) - U).^2).';
end
c(:,:,K+1) = 0;


% backward metric 
b = zeros(Ns, K+1);   %matrix

% the index has to go backwards
for k = K:-1:1      
    for i = 1:Ns
        % the index of the state is computed with the following
        % expression
        possible_state = mod(i-1, M^(L1 + L2 - 1))*M + 1;
        % the value of b is computed from b(k+1)
        b(i, k) = max(b(possible_state:possible_state+M-1, k+1) ...
            + c(:, i, k+1));
    end
end



% forward metric, state metric, log-likelihood function 
% f_old is set to -1
f_old = zeros(Ns, 1);   
f_new = zeros(Ns, 1); 
%will contain the symbols from which we select the one with the highest
%likelihood
likely = zeros(M, 1);
detected = zeros(K, 1);
row_step = (0:M-1)*M^(L1+L2-1);
for k = 1:K   
    for j = 1:Ns       
        in_vec = ceil(j/M) + row_step;
        f_new(j) = max(f_old(in_vec) + c(mod(j-1, 4)+1, in_vec, k).');
    end
    v = f_new + b(:, k);
    for beta = 1:M
        ind = states_symbols(:,M) == symb(beta);
        likely(beta) = max(v(ind));
    end
    [~, maxind] = max(likely);
    detected(k) = symb(maxind); 
    f_old = f_new;
end

toc(tStart)

end