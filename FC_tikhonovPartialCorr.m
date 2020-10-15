% Function to compute connectivity based on Tikhonov regularized partial correlation
% INPUT: time series data as a Txp matrix, where T=number of time points,
% p=number of brain regions; alpha = regularization parameter
% OUTPUT: Symmetric, weighted pxp adjacency matrix

function A = FC_tikhonovPartialCorr(t, alpha)

    % pausing until statistics toolbox is available
    while (~license('checkout', 'Statistics_Toolbox'))
        pause(30);
    end

    nParcels = size(t, 2);
    C = cov(t); % covariance matrix
    P = -inv(C + alpha*eye(nParcels)); % precision matrix after regularization
    A = (P ./ repmat(sqrt(abs(diag(P))), 1, nParcels)) ./ repmat(sqrt(abs(diag(P)))', nParcels, 1); % normalization

    % setting diagonal elements to zero
    A = A - (diag(diag(A)));
end