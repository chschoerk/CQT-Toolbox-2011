function loc = pickPeaksSparse(X,threshold)
%loc = pickPeaksSparse(X,extraBinsLeft,extraBinsRight,threshold)
%
%Finds peaks in vector X and returns their location.
%
%INPUT:
%   X           ... input vector (e.g. spectrum)
%   threshold   ... minimum amplitude of peaks

%OUTPUT:
%   loc         ... index of detected peaks
%
%Christian Schörkhuber 2011-08


af = abs(X(:,1)); %magnitudes of current frame  

dl1 = (af(2:end-1) - af(1:end-2)) > 0;
dr1 = (af(2:end-1) - af(3:end)) > 0;

val = af(2:end-1);
valcrit = val > threshold;
extraMask = true(size(val));

% if (extraBinsLeft || extraBinsRight) > 0
%     extraMask(1:extraBinsLeft-2) = false;
%     extraMask(end-extraBinsRight+2:end) = false;
% end

loc = find(dl1 & dr1 & valcrit & extraMask) + 1;

