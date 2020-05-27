% Function to compute connectivity based on mutual information in the time
% domain. Implementation used is the average mutual information, ami() function
% https://www.mathworks.com/matlabcentral/fileexchange/10040-average-mutual-information

% INPUT: 
% t: time series data as a Txp matrix, where T=number of time points,
% p=number of brain regions
%
% OUTPUT: Symmetric, weighted pxp adjacency matrix

function A = FC_mutualInformationTime(t)

    % pausing until statistics toolbox is available
    while (~license('checkout', 'Statistics_Toolbox'))
        pause(30);
    end

    nNodes = size(t, 2);
    A = zeros(nNodes);
    for i = 1:nNodes-1
        y1 = t(:, i);
        for j = i+1:nNodes
            y2 = t(:, j);
            MI = computeNMI(y1, y2);
            A(i, j) = MI;
        end
    end
    A = A + A';
end
