clc;
clear;
close all;


%% Constants definition
n       = 100;     %number of nodes
lambda  = 1;        %contamination intensity
x0      = 1;        %initial number of infected nodes
nEvents = n-x0;     %number of events that will occur (number of contaminations)

%% Building graph (complete topology)
Adj     = ones(n,n) - eye(n,n);
G       = graph(Adj);

%% Plot generated graph

figure
hold on;
axis off;
set(gca, 'DataAspectRatio', [1,1,1])
p1 = plot(G, '-.o','NodeLabel',{}, 'EdgeAlpha', 0.1,  'NodeColor', 'b', 'MarkerSize', 5);


%% Simulating behavior
[t, states, infectEdge] = simulateEvolutionSI(n, x0, Adj, lambda);

%% Movie
displayScenario(G, t, nEvents, states, infectEdge)


function [] = displayScenario(G, t, nEvents, states, infectEdge)

for i=1:nEvents
    clf
    % plot graph
    p = plot(G, '-.o','NodeLabel',{}, 'LineWidth', 0.0001,  'NodeColor', 'b', 'MarkerSize', 3);
    axis off
    % differentiate infected nodes
    infectedList = find(states(i,:)==1);
    highlight(p,infectedList,'NodeColor','red');
    
    % differentiate transmitter, receiver and edge
    highlight(p, infectEdge(i,:), 'NodeColor', [1 .5 0], 'EdgeColor',[1 1 0], 'LineWidth', 4, 'MarkerSize', 5)
    highlight(p, infectEdge(i,1), 'NodeColor','r', 'MarkerSize', 5)
    shg
    pause((t(i+1)-t(i))/100);
    
end
end
