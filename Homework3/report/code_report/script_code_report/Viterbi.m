function [detected] = Viterbi(r_c, hi, L1, L2, N1, N2)

if (L1 > N1) || (L2 > N2)
    disp('Check your input')
    return
end

M = 4;
symb = [1+1i, 1-1i, -1+1i, -1-1i]; % QPSK constellation
Kd = 28; % Size of the Trellis diagram (and of the matrix)
Ns = M ^ (L1+L2); % States
r_c  =  r_c(1+N1-L1 : end-N2+L2);   % Discard initial and final samples of r
hi = hi(1+N1-L1 : end-N2+L2);   % Discard initial and final samples of hi

survSeq = zeros(Ns, Kd);
detectedSymb = zeros(1, length(r_c));
cost = zeros(Ns, 1); % Define Gamma(-1) for each state (cost)

statelength = L1 + L2; % state length
statevec = zeros(1, statelength); % symb idx: old --> new
u_mat = zeros(Ns, M);
for state = 1:Ns
    for j = 1:M
        lastsymbols = [symb(statevec + 1), symb(j)]; % symbols: old --> new
        u_mat(state, j) = lastsymbols * flipud(hi);
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
    % Initialize the costs of the new states to -1
    costnew = - ones(Ns, 1);
    % Vector of the predecessors: the i-th element is the predecessor at
    % time k-1 of the i-th state at time k.
    pred = zeros(Ns, 1);
    % counter iteratively: (mod(state-1, M^(L1+L2-1)) * M + j).
    newstate = 0;
    for state = 1 : Ns  % All states
        for j = 1 : M   % M times
            % Index of the new state: it's mod(state-1, M^(L1+L2-1)) * M + j
            newstate = newstate + 1;
            if newstate > Ns, newstate = 1; end
            u = u_mat(state, j);
            % updatethe cost of the new state and overwrite the predecessor
			% if this transition has lower cost than before
            newstate_cost = cost(state) + abs(r_c(k) - u)^2;
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
        temp(newstate, 1:Kd) = ...
            [survSeq(pred(newstate), 2:Kd), ...
            symb(mod(newstate-1, M)+1)];
    end
    [~, decided_index] = min(costnew);	% Find the oldest symbol that yields the min cost
    detectedSymb(1+k) = survSeq(decided_index, 1); % and store it.
    survSeq = temp;
    cost = costnew;
end

detectedSymb(length(r_c)+2 : length(r_c)+Kd) = survSeq(decided_index, 1:Kd-1);
% Use the min cost from the last iteration
detectedSymb = detectedSymb(Kd+1 : end);
detected = detectedSymb;
detected = detected(2:end); % Discard first symbol (time k=-1)
end