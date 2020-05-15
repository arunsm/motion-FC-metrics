% Function to compute connectivity based on wavelet coherence.
% Implementation based on Grinsted toolbox
% Grinsted, A., J. C. Moore, S. Jevrejeva (2004), Application of the cross wavelet transform and wavelet coherence to geophysical time series, Nonlin. Process. Geophys., 11, 561566

% INPUT:
% t: time series data as a Txp matrix, where T=number of time points,
% p=number of brain regions; relaxation time
% tr: relaxation time
% bandpass: frequency range in which to average coherence values
%
% OUTPUT: Symmetric, weighted pxp adjacency matrix

function A = FC_waveletCoherence(t, bandpass, tr)
tic
nNodes = size(t,2);
cohij = zeros(nNodes);
for i = 1:nNodes-1
    for j = i+1:nNodes
        [Rsq,period,~] = wtc(t(:,i),t(:,j),'mcc',0);
        freq = 1./(period*tr);
        cohij(i,j) = mean(mean(Rsq(freq<bandpass(2)&freq>bandpass(1),:)));
    end
end

A = cohij +cohij';
toc
end
