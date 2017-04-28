function [t, states, infectEdge] = simulateEvolutionSI(n, x0, Adj, lambda)

% Initialize variables
nEvents    = n-x0;                      %number of events
t          = zeros(nEvents+1, 1);       %list of times when an event occurs
states     = zeros(nEvents+1, n);       %represents states of nodes  between two events (1 if infected, 0 otherwise)
infectEdge = zeros(nEvents, 2);         %index of transmitter and receiver
x          = x0;                        %number of infected notes

% Select first contaminated edges
states(1,randi([1 n], x0, 1)) = 1;


% For loop simulating events
for i=1:nEvents
    
    % Sample a random time (exponential distribution)
    q      = lambda*x*(n-x)/(n-1);         %time rate
    t(i+1) = t(i)+ exprnd(1/q);                  %event sampling
    
    % Sample a random infected node (uniform distribution on infected nodes)
    infectedList        = find(states(i,:)==1);                      %retrieve list of infected nodes
    transmitterIndex    = randi(length(infectedList), 1,1);    %sample transmitter index
    transmitter         = infectedList(transmitterIndex);
    
    
    % Sample a random sane neighbor for the infected node
    saneNeighborsList   = getSaneNeighbors(Adj, transmitter, states(i,:));                       % retrieve list of sane neighbors
    receiverIndex       = randi(length(saneNeighborsList), 1, 1);
    receiver            = saneNeighborsList(receiverIndex);
    
    % Apply contamination
    x = x+1;
    states(i+1, :)        = states(i,:);
    states(i+1, receiver) = 1;
    infectEdge(i,1)       = transmitter;
    infectEdge(i,2)       = receiver;
    
end

end



function [list] = getSaneNeighbors(adj, transmitter, states)

% list of neighbors
neighbors = find(adj(transmitter, :) == 1);

% list of sane neighbors
list = neighbors(states(neighbors) == 0);
end