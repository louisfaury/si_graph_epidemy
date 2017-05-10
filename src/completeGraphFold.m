clc;
clear;
close all;

set(0,'defaulttextInterpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                       %%   
%%              Check logarithmic bound on expectation                   %%
%%                                                                       %%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants definition
n       = round(logspace(1, 3, 10));      %tested population cardinals
fold    = 100;                              %number of tests per n value
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
stdTime  = sqrt(varTime);

%% Linear regression over lambda 
P = [log(n)',ones(size(log(n)'))];
lambda = 2 / ( [1 0]* ((P'*P)\(P'*meanTime'))) ; % <- normal equations 
disp(lambda);

%% Plot results
figure 
set(0,'defaulttextInterpreter','latex')
set(gca, 'FontSize', 14)
a = semilogx(n, meanTime, 'b', 'LineWidth', 2);
b = jbfill(n,meanTime-stdTime,meanTime+stdTime,'r','r',1,0.3);
xlabel('n','FontSize',14);
ylabel('E(T)','FontSize',14)
grid minor
title('$E[T] = f(log(n))$','FontSize',15,'FontWeight','bold','interpreter','latex')
h = legend([a,b],{'Mean fold value','Std fold tube'},'Location','southeast');
set(h,'FontSize',14);

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
figure
set(gca, 'FontSize', 14)
h = histogram(fluct, 100);
xlabel('Flucutation','FontSize',14)
title('Histogram representation of the fluctuation distribution ($n=100, \quad \lambda = 1$)','FontSize',14)


%% Compute cumulative distribution
cdf = zeros(1, size(h.BinCounts, 2));
abs = zeros(1, size(h.BinCounts, 2));
cdf(1,1) = h.BinCounts(1, 1);
abs(1,1) = h.BinEdges(1, 1)+0.5*h.BinWidth;
for i=2:size(cdf,2)
    cdf(i) = cdf(1, i-1)+h.BinCounts(1, i);
    abs(i) = h.BinEdges(1, i)+0.5*h.BinWidth;
end

% Normalization
cdf = cdf/cdf(1,end);

oneMinusCdf = ones(size(cdf))-cdf;

% Plot
figure
grid minor
set(gca, 'FontSize', 16)
plot(abs, oneMinusCdf, 'LineWidth', 2);
xlabel('t','FontSize',14)
ylabel('P($S_n \geq t$)','FontSize',14)
title('P($S_n \geq t$) = f(t)','FontSize',14)

