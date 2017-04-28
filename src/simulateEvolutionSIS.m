function [nEvents, t, states, infectEdge, absorbed] = simulateEvolutionSIS(n, x0, Adj, beta, delta)

% Initialize variables
nEvents    = 0;                         %number of events
t          = 0;                         %list of times when an event occurs
states     = zeros(1, n);               %represents states of nodes  between two events (1 if infected, 0 otherwise)
infectEdge = zeros(1, 2);               %index of transmitter and receiver
x          = x0;                        %number of infected notes

% Select first contaminated edges
list = randperm(n);
states(1,list(1:x0)) = 1;


% While loop simulating events
i = 1;
while true
    eventPerformed = 0;
    
    % Sample a random times for each infected node
    infected = find(states(i,:)==1);
    infectionAttempts = t(i,1) + exprnd(1/beta, size(infected));
    remissionAttempts = t(i,1) + exprnd(1/delta, size(infected));
    
    while eventPerformed == 0
        
        % find minimum time
        [minInfectionTime, transmitterIndex]    = min(infectionAttempts);
        [minRemissionTime, remittedIndex] = min(remissionAttempts);
        
        transmitter = infected(transmitterIndex);
        remitted    = infected(remittedIndex);
        
        if minRemissionTime <= minInfectionTime
            % applying a remission
            nEvents                     = nEvents + 1;
            x                           = x-1;
            t(i+1,1)                    = minRemissionTime;
            states(i+1, :)              = states(i,:);
            states(i+1, remitted(1))    = 0;
            infectEdge(i,1)             = remitted(1);
            infectEdge(i,2)             = remitted(1);
            
            eventPerformed = 1;
            
        else
            % check if infection is successful
            neighbors     = getNeighbors(Adj, transmitter);
            saneNeighbors = getSaneNeighbors(Adj, transmitter, states(i,:));
            
            if(length(neighbors)>0)
                infectedIndex = randi([1, length(neighbors)], 1,1);
                if any(saneNeighbors==neighbors(infectedIndex))
                    
                    % Infection successfull
                    receiver              = neighbors(infectedIndex);
                    nEvents               = nEvents + 1;
                    x                     = x+1;
                    t(i+1,1)              = minInfectionTime;
                    states(i+1, :)        = states(i,:);
                    states(i+1, receiver) = 1;
                    infectEdge(i,1)       = transmitter;
                    infectEdge(i,2)       = receiver;
                    eventPerformed = 1;
                else
                    
                    % Infection failed, resample a new infection time for this
                    % node
                    
                    newInfectAttempt = exprnd(1/beta);
                    infectionAttempts(transmitterIndex) = infectionAttempts(transmitterIndex) + newInfectAttempt;
                    eventPerformed = 0;
                end
            else
                newInfectAttempt = exprnd(1/beta);
                infectionAttempts(transmitterIndex) = infectionAttempts(transmitterIndex) + newInfectAttempt;
                eventPerformed = 0;
            end
            
        end
    end
    
    %% Conditions for end of simulation
    if isequal(states(i+1, :),ones(1,n))
        absorbed = -1;
        break;
    elseif isequal(states(i+1, :), zeros(1,n))
        absorbed = 1;
        break;
    elseif nEvents>3*n
        absorbed = 0;
        break;
    end
    
    i=i+1;
    
end

end



function [list] = getSaneNeighbors(adj, transmitter, states)

% list of neighbors
neighbors = find(adj(transmitter, :) == 1);

% list of sane neighbors
list = neighbors(states(neighbors) == 0);
end


function [neighbors] = getNeighbors(adj, transmitter)

% list of neighbors
neighbors = find(adj(transmitter, :) == 1);

end