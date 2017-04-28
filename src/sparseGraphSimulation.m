clc;
clear;
close all;


%% Constants definition
n               = 500;          %number of nodes
nCommunities    = 5;            %number of communities
beta            = 0.25;         %contamination intensity
delta           = 0.5;          %remission intensity
x0              = 500;          %initial number of infected nodes
propEdge        = 0.3;          %proportion of edges inside clusters

%% Building graph (sparse topology)
Adj     = generateSparseGraph(n, nCommunities, propEdge);
G       = graph(Adj);
R       = max(abs(eig(double(Adj))));
disp(R);

%% Plot generated graph

figure
hold on;
axis off;
set(gca, 'DataAspectRatio', [1,1,1])
p1 = plot(G, '-.o','NodeLabel',{}, 'EdgeAlpha', 0.1,  'NodeColor', 'b', 'MarkerSize', 3);


%% Simulating behavior
[nEvents, t, states, infectEdge, absorbed] = simulateEvolutionSIS(n, x0, Adj, beta, delta);

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
    if infectEdge(i,1) ~= infectEdge(i,2)
        highlight(p, infectEdge(i,:), 'NodeColor', [1 .5 0], 'EdgeColor',[1 1 0], 'LineWidth', 4, 'MarkerSize', 5)
        highlight(p, infectEdge(i,1), 'NodeColor','r', 'MarkerSize', 5)
    else
         highlight(p, infectEdge(i,1), 'NodeColor','g', 'MarkerSize', 5)
    end
    
    shg
    pause((t(i+1)-t(i))/10);
    
end

end
