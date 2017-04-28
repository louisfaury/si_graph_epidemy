function [Adj] = generateSparseGraph(n, nCommunities)

Adj = zeros(n,n);

%% Create communities
sizeStdDev      = 0.05*n;
meanSize        = round(n/nCommunities);
communitySize   = zeros(1, nCommunities);

while true
    
    for i=1:nCommunities-1
        communitySize(1, i) = round(rand(1,1)*sizeStdDev + meanSize);
    end
    
    communitySize(1, nCommunities) = n - sum(communitySize(1,1:nCommunities-1));
    if communitySize(1, nCommunities)>0
        break
    end
    
end

%% Build edges inside communities
minIndex = 1;
for i=1:nCommunities
    
    maxIndex = minIndex + communitySize(1,i)-1;
    Adj(minIndex:maxIndex, minIndex:maxIndex) = randi([0,1], maxIndex-minIndex + 1,  maxIndex-minIndex + 1);
    minIndex = maxIndex+1;
    
end

%% Add a few edges between communities
sparseMat = (randi([0,100],n,n)>95);
Adj = Adj + double(sparseMat);
Adj = floor((Adj+Adj')/2);
Adj = (Adj>0);

%% Remove self links
Adj(logical(eye(size(Adj)))) = 0;

end