function edgeWeights = computeEdgeWeights(A)

    edgeWeights = squareform(A);
%     edgeWeights = [];
%     ctr = 1;
%     % finding edge weights
%     for j = 1:size(A, 1)
%         for k = j:size(A, 2) % assuming symmetric matrix
%             % skip diagonal elements
%             if j == k
%                 continue;
%             end
%             
%             edgeWeights(ctr, 1) = j; % node 1 index
%             edgeWeights(ctr, 2) = k; % node 2 index
%             edgeWeights(ctr, 3) = A(j, k); % edge weight
%             
%             ctr = ctr + 1;
%         end
%     end
end