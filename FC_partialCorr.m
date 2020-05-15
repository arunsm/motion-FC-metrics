% Function to compute connectivity based on partial correlation
% INPUT: time series data as a Txp matrix, where T=number of time points,
% p=number of brain regions
% OUTPUT: Symmetric, weighted pxp adjacency matrix

function A = FC_partialCorr(t)
    A = partialcorr(t);
    % setting diagonal elements to zero
    for j = 1:size(A, 1)
        for k = j:size(A, 2)
            if j == k
                A(j, k) = 0;
            end
        end
    end
end