clc;
clear;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                       %%   
%%              Check logarithmic bound on expectation                   %%
%%                                                                       %%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants definition
n       = round(logspace(1, 3, 100));      %tested population cardinals
fold    = 30;                              %number of tests per n value
lambda  = 1;                               %contamination intensity
x0      = 1;                               %initial number of infected nodes


%% Variables definitions
totalTime = zeros(fold ,size(n, 2));    %store times before total contamination

%% Launch all simulations
for i=1:size(n,2)
    disp(i)
    
    % Building graph (complete topology)
    Adj     = ones(n(1,i),n(1,i)) - eye(n(1,i),n(1,i));
    
    for f=1:fold
        [t, ~, ~] = simulateEvolutionSI(n(1,i), x0, Adj, lambda);
        totalTime(f, i) = t(end,1);
    end
    
end

%% Compute contamination time mean
meanTime = mean(totalTime,1);
varTime  = var(totalTime, 1);



%% Plot results
figure 

semilogx(n, meanTime, 'b', 'LineWidth', 2);
xlabel('population size');
ylabel('time')
grid minor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                       %%   
%%                     Check fluctuation                                 %%
%%                                                                       %%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

%% Constants definition
n       = 100;                             %tested population size
fold    = 10000;                           %number of tests
lambda  = 1;                               %contamination intensity
x0      = 1;                               %initial number of infected nodes


%% Variables definitions
totalTime = zeros(1 ,fold);    %store times before total contamination

% Building graph (complete topology)
Adj     = ones(n,n) - eye(n,n);


for f=1:fold
    disp(f)
    [t, ~, ~] = simulateEvolutionSI(n, x0, Adj, lambda);
    totalTime(1, f) = t(end,1);
    
end

%% Compute contamination time mean and fluctuations
meanTime = mean(totalTime,2);
fluct    = totalTime - meanTime;


%% Plot histogram
hist(abs(fluct), 100);
