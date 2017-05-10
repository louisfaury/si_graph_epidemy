clc;
clear;
close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                       %%   
%%              Check effect of population size                          %%
%%                                                                       %%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Constants definition
fold            = 100;                   %number of iterations per config
n               = round(logspace(2,3,20));  % list of tested population sizes
nCommunities    = 8;                    %number of communities
beta            = 0.025;                 %contamination intensity
delta           = 0.5;                  %remission intensity
x0              = round(n/2);           %initial number of infected nodes
targetR         = 0.5;                    %target graph spectral radius       

%% Variables initialization
totalTime = zeros(fold ,size(n, 2));    %store times before total contamination


for i=1:length(n)
    for f=1:fold
        
        propEdge =1/ n(1,i)*nCommunities*targetR;
        %% Building graph (sparse topology)
        Adj     = generateSparseGraph(n(1,i), nCommunities, propEdge);
        R       = max(abs(eig(double(Adj))));
        disp(R);
        %% Simulating
        [t, ~, ~, ~, absorbed] = simulateEvolutionSIS(n(1,i),x0, Adj, beta, delta);
        disp(absorbed)
        totalTime(f, i) = t(end,1);
   
    end
    
end


%% Compute contamination time mean
meanTime = mean(totalTime,1);
varTime  = var(totalTime, 1);
stdTime = sqrt(varTime);

%% Plot results
figure 
a = plot(n, meanTime, 'b', 'LineWidth', 2);
b = jbfill(n,meanTime-stdTime,meanTime+stdTime,'r','r',1,0.3);
set(gca, 'FontSize', 14)
set(0,'defaulttextInterpreter','latex')
xlabel('n');
ylabel('T')
title('$E[T] = f(n)$')
grid minor
h = legend([a,b],{'Mean fold value','Std fold tube'},'Location','southeast');
set(h,'FontSize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                       %%   
%%              Check effect of spectral radius                          %%
%%                                                                       %%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

%% Constants definition
fold            = 100;                  %number of iterations per config
R               = 1:2:25;               % list of tested spectral radius
n               = 500; 
nCommunities    = 4;                    %number of communities
beta            = 0.05;                 %contamination intensity
delta           = 0.5;                  %remission intensity
x0              = round(n/2);           %initial number of infected nodes
       

%% Variables initialization
totalTime = zeros(fold ,size(n, 2));    %store times before total contamination


for i=1:length(R)
    actualR = 0;
    for f=1:fold
        
        propEdge =1/ n*nCommunities*R(i);
        %% Building graph (sparse topology)
        Adj           = generateSparseGraph(n, nCommunities, propEdge);
        actualR       =  actualR + max(abs(eig(double(Adj))))/fold;

        %% Simulating
        [t, ~, ~, ~, absorbed] = simulateEvolutionSIS(n,x0, Adj, beta, delta);
        disp(absorbed)
        disp( max(abs(eig(double(Adj)))))
        totalTime(f, i) = t(end,1);
   
    end
    meanR(i) = actualR;
    
end


%% Compute contamination time mean
meanTime = mean(totalTime,1);
varTime  = var(totalTime, 1);

%% Plot results
figure 
set(0,'defaulttextInterpreter','latex')
set(gca, 'FontSize', 14)
plot(R, meanTime, 'b', 'LineWidth', 2);
xlabel('$\rho$');
ylabel('T')
title('$E[T] = f(\rho)$')
grid minor
