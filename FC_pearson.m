% Function to compute connectivity based on Pearson correlation
% INPUT: time series data as a Txp matrix, where T=number of time points,
% p=number of brain regions
% OUTPUT: Symmetric, weighted pxp adjacency matrix

function A = FC_pearson(t)

    % pausing until statistics toolbox is available
    while (~license('checkout', 'Statistics_Toolbox'))
        pause(30);
    end

    A = corr(t);
    
    % setting diagonal elements to zero
    A = A - (diag(diag(A)));
end
