% Function to compute connectivity based on Spearman correlation
% INPUT: time series data as a Txp matrix, where T=number of time points,
% p=number of brain regions
% OUTPUT: Symmetric, weighted pxp adjacency matrix

function A = FC_spearman(t)
    
    % pausing until statistics toolbox is available
    while (~license('checkout', 'Statistics_Toolbox'))
        pause(30);
    end

    A = corr(t, 'Type', 'Spearman');    
    
    % setting diagonal elements to zero
    A = A - (diag(diag(A)));
end
