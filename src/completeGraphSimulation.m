clc;
clear;
close all;

%% Do you want to record ?
recordFlag = 1;

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
displayScenario(G, t, nEvents, states, infectEdge, recordFlag)


function [] = displayScenario(G, t, nEvents, states, infectEdge, recordFlag)

if recordFlag 
    video = VideoWriter('completeGraph.avi');
    %video.CompressionRatio = 100;
    video.FrameRate = 10;
    open(video);
end

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
    
    if recordFlag
          counter = round((t(i+1)-t(i))*video.FrameRate);
          for k=1:counter
            frame = getframe(gcf);
            writeVideo(video,frame);
          end
    end
    
    pause((t(i+1)-t(i))/10);
        
end

if recordFlag
    close(video)
end


    
end
