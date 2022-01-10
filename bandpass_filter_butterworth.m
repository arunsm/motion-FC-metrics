function filteredData = bandpass_filter_butterworth(inputData, tr, f1, f2)
    
    % pause processing until signal toolbox is available
    while (~license('checkout', 'Signal_Toolbox'))
	pause(30);
    end
    
    filterOrder = 5;
    fs = 1/tr; % sampling rate in Hz
    nyq = 0.5*fs; % Nyquist frequency
    lowpass = f1/nyq;
    highpass = f2/nyq;
    
    % using b, a coefficients
    [b, a] = butter(filterOrder, [lowpass highpass], 'bandpass');
    
    % using z, p, k method
%     [z,p,k] = butter(filterOrder, [lowpass highpass], 'bandpass');
%     [sos, g] = zp2sos(z,p,k);
    
    filteredData = filtfilt(b, a, inputData);
end