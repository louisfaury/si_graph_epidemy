clc;
clear;
close all;


%% Constants definition
n       = 200;      %number of nodes
nCommunities = 4;   %number of communities
beta  = 1;        %contamination intensity
delta    = 0.01;     %remission intensity
x0      = 1;        %initial number of infected nodes

%% Building graph (sparse topology)
Adj     = generateSparseGraph(n, nCommunities);
G       = graph(Adj);

%% Plot generated graph

figure
hold on;
axis off;
set(gca, 'DataAspectRatio', [1,1,1])
p1 = plot(G, '-.o','NodeLabel',{}, 'EdgeAlpha', 0.1,  'NodeColor', 'b', 'MarkerSize', 3);


%% Simulating behavior
[nEvents, t, states, infectEdge] = simulateEvolutionSIS(n, x0, Adj, beta, delta);

%% Movie
displayScenario(G, t, nEvents, states, infectEdge)


function [] = displayScenario(G, t, nEvents, states, infectEdge)

for i=1:nEvents
    clf
    % plot graph
    p = plot(G, '-.o','NodeLabel',{}, 'LineWidth', 0.0001,  'NodeColor', 'b', 'MarkerSize', 2);
    axis off
    % differentiate infected nodes
    infectedList = find(states(i,:)==1);
    highlight(p,infectedList,'NodeColor','red');
    
    % differentiate transmitter, receiver and edge
    highlight(p, infectEdge(i,:), 'NodeColor','g','EdgeColor','g', 'LineWidth', 2, 'MarkerSize', 2)
    shg
    pause(t(i+1)-t(i));
    
end

end
