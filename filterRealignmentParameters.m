function filteredRP = filterRealignmentParameters(RP)
    % values of cutoff frequency and bandwidth derived from Fair et al,
    % NeuroImage (2020) and Power et al, NeuroImage (2019)
    cutoffFrequency = 0.3; % Hz
    bandwidth = 0.4; % Hz 
    [b, a] = iirnotch(cutoffFrequency, bandwidth);
    filteredRP = filtfilt(b, a, RP);
end