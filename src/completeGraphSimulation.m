clc;
clear;
close all;


%% Constants definition
n       = 1000;      %number of nodes
lambda  = 1;      %contamination intensity
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
[t, states, infectEdge] = simulateEvolution(n, x0, Adj, lambda);

%% Movie
displayScenario(G, t, nEvents, states, infectEdge)


function [] = displayScenario(G, t, nEvents, states, infectEdge)

for i=1:nEvents
    clf
    % plot graph
    p = plot(G, '-.o','NodeLabel',{}, 'LineWidth', 0.0001,  'NodeColor', 'b', 'MarkerSize', 5);
    axis off
    % differentiate infected nodes
    infectedList = find(states(i,:)==1);
    highlight(p,infectedList,'NodeColor','red');
    
    % differentiate transmitter, receiver and edge
    highlight(p, infectEdge(i,:), 'NodeColor','g','EdgeColor','g', 'LineWidth', 2, 'MarkerSize', 10)
    shg
    pause(t(i+1)-t(i));
    
end
end
