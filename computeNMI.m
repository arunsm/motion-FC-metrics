% function to compute normalized mutual information between two time series
% x and y

function NMI = computeNMI(x,y)

% converting x and y to column vectors
x=x(:);
y=y(:);
n=length(x);

if n~=length(y)
    error('x and y should be same length');
end

% compute number of bins using the Freedman-Diaconis rule
nbinsx = ceil((max(x) - min(x))/(2*iqr(x)*n^(-1/3)));
nbinsy = ceil((max(y) - min(y))/(2*iqr(y)*n^(-1/3)));
nbins = ceil((nbinsx + nbinsy)/2);

% normalize data to [0 1)
x=x-min(x);
x=x*(1-eps)/max(x);
y=y-min(y);
y=y*(1-eps)/max(y);

Pxy = hist3([x, y], [nbins, nbins]); % histogram counts of joint distribution
Pxy = Pxy/n; % divide by length of time series to obtain probabilities
Pxy = Pxy + eps; % avoid division and log of zero
Px = sum(Pxy,2);
Py = sum(Pxy,1);

Hx = -sum(Px.*log2(Px)); % entropy of x
Hy = -sum(Py.*log2(Py)); % entropy of y
Hxy = -sum(Pxy(:).*log2(Pxy(:))); % entropy of joint distribution
MI = Hx + Hy - Hxy; % mutual information
NMI = MI/sqrt(Hx*Hy); % normalized mutual information

end