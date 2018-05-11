function [detected] = viterbi(rho, psi, N1, N2, L1, L2)

M = 4;
syms = [1+1i, 1-1i, -1+1i, -1-1i];
Ns = M^(L1+L2);
Kd = 30;
rho = rho(1+N1-L1:end-N2+L2);
psi = psi(1+N1-L1:end-N2+L2);

survSeq = zeros(Ns, Kd);
detectedSymb = zeros(1, length(rho));
cost = zeros(Ns, 1);

statelength = L1 + L2; % number of digits of the state, i.e. its length in base M=4
statevec = zeros(1, statelength); % symbol index, from the oldest to the newest
u_mat = zeros(Ns, M);
for state = 1:Ns
    
    % Set value of the current element of u_mat
    for j = 1:M
        lastsymbols = [syms(statevec + 1), syms(j)]; % symbols, from the oldest to the newest
        %         u_mat(state, j) = lastsymbols * flipud(psi);
    end
    
    % Update statevec
    statevec(statelength) = statevec(statelength) + 1;
    i = statelength;
    while (statevec(i) >= M && i > 1)
        statevec(i) = 0;
        i = i-1;
        statevec(i) = statevec(i) + 1;
    end
end

for k = 1 : length(rho)
    
    % Initialize the costs of the new states to -1
    costnew = - ones(Ns, 1);
    
    % Vector of the predecessors: the i-th element is the predecessor at
    % time k-1 of the i-th state at time k.
    pred = zeros(Ns, 1);
    
    % Counter for the new state. It is determined iteratively even though a
    % closed form expression exists: (mod(state-1, M^(L1+L2-1)) * M + j).
    newstate = 0;
    
    
    for state = 1 : Ns  % Cycle through all states, at time k-1
        
        for j = 1 : M   % M possibilities for the new symbol
            
            % Index of the new state: it's mod(state-1, M^(L1+L2-1)) * M + j
            newstate = newstate + 1;
            if newstate > Ns, newstate = 1; end
            
            % Desired signal u assuming the input sequence is the one given by the current
            % state "state", followed by a new assumed symbol given by j.
            u = u_mat(state, j);
            
            % Compute the cost of the new state assuming this input sequence, then update
            % the cost of the new state, and overwrite the predecessor, if this transition
            % has a lower cost than before.
            newstate_cost = cost(state) + abs(rho(k) - u)^2;
            if costnew(newstate) == -1 ...     % not assigned yet, or...
                    || costnew(newstate) > newstate_cost  % ...found path with lower cost
                costnew(newstate) = newstate_cost;
                pred(newstate) = state;
            end
        end
    end
    
    
    % Update the survivor sequence by shifting the time horizon of the matrix by one, and
    % rewrite the matrix with the new survival sequences sorted by current state.
    % Meanwhile, decide the oldest sample (based on minimum cost) and get rid of it to
    % keep only Kd columns in the matrix.
    temp = zeros(size(survSeq));
    for newstate = 1:Ns
        temp(newstate, 1:Kd) = ...    % In the new matrix except the last col
            [survSeq(pred(newstate), 2:Kd), ... % we write the data we had except the oldest one,
            syms(mod(newstate-1, M)+1)];        % and then the new symbol we just supposed to have received.
    end
    [~, decided_index] = min(costnew);      % Find the oldest symbol that yields the min cost
    detectedSymb(1+k) = survSeq(decided_index, 1); % and store it (decide it for good).
    survSeq = temp;
    
    % Update the cost to be used as cost at time k-1 in the next iteration
    cost = costnew;
end

detectedSymb(length(rho)+2 : length(rho)+Kd) = survSeq(decided_index, 1:Kd-1);
% Decided using the min cost from the last iteration
detectedSymb = detectedSymb(Kd+1 : end);
detected = detectedSymb;


end