function [detected] = VBA(r_c, hi, L1, L2, N1, N2)

M = 4;
symb = [1+1i, 1-1i, -1+1i, -1-1i]; % All possible transmitted symbols (QPSK)
Kd = 28; 
Ns = M ^ (L1+L2); % Number of states
r_c  =  r_c(1+N1-L1 : end-N2+L2);   % Discard initial and final samples of r
hi = hi(1+N1-L1 : end-N2+L2);   % Discard initial and final samples of hi

tStart = tic; 

survSeq = zeros(Ns, Kd);
detectedSymb = zeros(1, length(r_c));
cost = zeros(Ns, 1); 

statelength = L1 + L2; 
statevec = zeros(1, statelength); 
%matrix with the input values
U = zeros(Ns, M);
for state = 1:Ns
    for j = 1:M
        lastsymbols = [symb(statevec + 1), symb(j)]; 
        U(state, j) = lastsymbols * flipud(hi);
    end
    statevec(statelength) = statevec(statelength) + 1;
    i = statelength;
    while (statevec(i) >= M && i > 1)
        statevec(i) = 0;
        i = i-1;
        statevec(i) = statevec(i) + 1;
    end
end

for k = 1 : length(r_c)

    nextcost = - ones(Ns, 1);
    pred = zeros(Ns, 1);
    nextstate = 0;
    
    for state = 1 : Ns
        
        for j = 1 : M   
            nextstate = nextstate + 1;
            if nextstate > Ns, nextstate = 1; end
            u = U(state, j);
            newstate_cost = cost(state) + abs(r_c(k) - u)^2;
            if nextcost(nextstate) == -1 ...    
                    || nextcost(nextstate) > newstate_cost 
                nextcost(nextstate) = newstate_cost;
                pred(nextstate) = state;
            end
        end
    end
    
    temp = zeros(size(survSeq));
    for nextstate = 1:Ns
        temp(nextstate, 1:Kd) = ...
            [survSeq(pred(nextstate), 2:Kd), ... 
            symb(mod(nextstate-1, M)+1)];       
    end
    [~, decided_index] = min(nextcost);   
    detectedSymb(1+k) = survSeq(decided_index, 1); 
    survSeq = temp;
    
    % Update the cost to be used as cost at time k-1 in the next iteration
    cost = nextcost;
end

toc(tStart)

detectedSymb(length(r_c)+2 : length(r_c)+Kd) = survSeq(decided_index, 1:Kd-1);

detectedSymb = detectedSymb(Kd+1 : end);
detected = detectedSymb;
detected = detected(2:end); 
end