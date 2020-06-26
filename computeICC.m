function IntraClassCorrelation = computeICC(data)

[n,k] = size(data); % n = number of subjects; k = number of scans

% mean per target
mpt = mean(data,2);
% get total mean
tm = mean(data(:));
% within target sum sqrs
tmp = (data - repmat(mpt,1,k)).^2;
WSS = sum(tmp(:));
% within target mean sqrs
WMS = WSS / (n*(k - 1));
% between target sum sqrs
BSS = sum((mpt - tm).^2) * k;
% between targets mean squares
BMS = BSS / (n - 1);

IntraClassCorrelation = (BMS - WMS) / (BMS + (k - 1) * WMS);
end