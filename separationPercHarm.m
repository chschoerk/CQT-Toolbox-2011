function [H P] = separationPercHarm(Xcq,lh,lp,octaveNum,a)
%[H P] = separationPercHarm(Xcq,lh,lp,octaveNum,a)
%
%Divides the CQT representation of a time domain signal into a percussive
%and a harmonic part. The separation approach is based on median filtering
%in horizontal and vertical direction. This implementation is an adaption of 
%the STFT based approach proposed by D. FitzGerald in 'Harmonic/Percussive
%Separation using Median Filtering', Dublin Institute of Technology, 2010.
%
%INPUT:
%    Xcq       ... structure returned by function cqt(...)
%    lh        ... length of horizontal median filter
%    lp        ... length of vertical median filter
%    octaveNum ... number of octaves (starting from the lowest octave) that
%                  should be considered (for CQT based time-/pitch-scale
%                  modifications only transients at lower frequencies are
%                  prone to transient smearing)
%    a         ... parameter for mask (usually a= 1 or 2)
%
%OUTPUT:
%   H           ... CQT representation of harmonic part
%   P           ... CQT representation of percussive part
%
%Christian Schörkhuber, 2011-11

bins = Xcq.bins;

X = Xcq.spCQT(1:octaveNum*bins,:);
A = abs(X).^2;

%% generate harmonic enhanced representation
H = zeros(size(A));

for k = 1:size(A,1)
    oct = ceil((size(Xcq.spCQT,1)-k+1) / bins);
    tvec = 1:2^(oct-2):size(A,2);
    slice = A(k,tvec);
    H(k,tvec) = medfilt1(slice,lh);    
end

%% generate percussion enhanced representation
minStep = 2^(Xcq.octaveNr-octaveNum-1);

BVec = 1:(octaveNum*bins);
P = zeros(size(A));

for n = 1:minStep:size(A,2) 
    tv = (n-1)./(minStep.*2.^(0:octaveNum-1));
    tv = tv - floor(tv);
    tv = tv(tv==0);
    lowerOctsNr = length(tv);
    binVec = BVec(end-lowerOctsNr*bins+1:end);   
    frame = A(binVec,n);
    P(binVec,n) = medfilt1(frame,lp);
end

%% generate mask
Mh = H.^a./(H.^a+P.^a+eps);

Mp = P.^a./(H.^a+P.^a+eps);
P = zeros(size(Xcq.spCQT));
P(1:size(X,1),:) = X.*Mp;

H = Xcq.spCQT;
H(1:size(X,1),:) = X.*Mh;

