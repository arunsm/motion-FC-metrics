function AdjacencyMatrix = computeFunctionalConnectivity(timeSeries, FC_method, f1, f2, tr)

switch FC_method
    case "Pearson"
        AdjacencyMatrix = FC_pearson(timeSeries); % Pearson correlations among time series
    case "PartialCorrelation"
        AdjacencyMatrix = FC_partialCorr(timeSeries); % partial correlations among time series
    case "Spearman"
        AdjacencyMatrix = FC_spearman(timeSeries); % Spearman correlations among time series
    case "MutualInformationTime"
        AdjacencyMatrix = FC_mutualInformationTime(timeSeries); % mutual information in time domain
    case "Coherence"
        AdjacencyMatrix = FC_coherence(timeSeries, tr, f1, f2); % average coherence in frequency band [f1, f2]
    case "MutualInformation"
        AdjacencyMatrix = FC_mutualInformation(timeSeries, tr, f1, f2); % mutual information in frequency band [f1, f2]
    case "WaveletCoherence"
        AdjacencyMatrix = FC_waveletCoherence(timeSeries, [f1 f2], tr); % wavelet coherence in frequency band [f1, f2]
    otherwise
        fprintf("Unknown method")
end

%AdjacencyMatrix(AdjacencyMatrix<0) = 0; % setting negative correlations to zero
%AdjacencyMatrix = atanh(AdjacencyMatrix); % Fisher z-transform

end
