% function to calculate fingerprinting accuracy from functional
% connectivity matrices

% Input: 
% databaseEdgeWeights - cell array of size (1, N) where N is the number of
% subjects in database; each entry in the cell array should be a vector of
% edge weights
% targetEdgeWeights - vector containing edge weights of individual who is
% being fingerprinted
% correct_idx - index of subject in target array

% Output:
% accuracy - 1 if individual correctly identified, 0 if not

function accuracy = calculateFingerprintAccuracy(databaseEdgeWeights, targetEdgeWeights, correct_idx)
    nSubjects = size(databaseEdgeWeights, 1);   
    edgeWeightCorrelations = zeros(1, nSubjects);
    for i = 1:nSubjects
        currentDatabaseEdgeWeights = databaseEdgeWeights{i};
        edgeWeightCorrelations(i) = corr(targetEdgeWeights', currentDatabaseEdgeWeights');
    end
    
    [~, idx] = max(edgeWeightCorrelations);
    if idx == correct_idx
        accuracy = 1;
    else
        accuracy = 0;
    end
end