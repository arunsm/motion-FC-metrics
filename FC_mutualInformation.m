% Function to compute connectivity based on mutual information in the frequency
% domain. Implementation used is the mutualinf() function described in
% Zhou, Dongli, Wesley K. Thompson, and Greg Siegle. "MATLAB toolbox for functional connectivity." Neuroimage 47.4 (2009): 1590-1607.

% INPUT: 
% t: time series data as a Txp matrix, where T=number of time points,
% p=number of brain regions; relaxation time
% tr: relaxation time
% f1, f2: frequency range in which to average connectivity values
%
% OUTPUT: Symmetric, weighted pxp adjacency matrix

function A = FC_mutualInformation(t, tr, f1, f2)

    % pause processing until signal toolbox is available
    while (~license('checkout', 'Signal_Toolbox'))
        pause(30);
    end

    nNodes = size(t, 2);
    A = zeros(nNodes);
    for i = 1:nNodes-1
        y1 = t(:, i);
        for j = i+1:nNodes
            y2 = t(:, j);
            A(i, j) = mutualinf(y1, y2, 1/tr, f1, f2);
        end
    end
    
    A = A + A';
end
