% Function to compute connectivity based on Tikhonov regularized partial correlation
% INPUT: time series data as a Txp matrix, where T=number of time points,
% p=number of brain regions; alpha = regularization parameter
% OUTPUT: Symmetric, weighted pxp adjacency matrix

function A = FC_tikhonovPartialCorr(t, alpha)

    % pausing until statistics toolbox is available
    while (~license('checkout', 'Statistics_Toolbox'))
        pause(30);
    end

    nBrainRegions = size(t, 2);
    C = cov(t); % covariance matrix
    P = -inv(C + alpha*eye(nBrainRegions)); % precision matrix after regularization
    A = (P ./ repmat(sqrt(abs(diag(P))), 1, nBrainRegions)) ./ repmat(sqrt(abs(diag(P)))', nBrainRegions, 1); % normalization

    % setting diagonal elements to zero
    for j = 1:size(A, 1)
        for k = j:size(A, 2)
            if j == k
                A(j, k) = 0;
            end
        end
    end
end