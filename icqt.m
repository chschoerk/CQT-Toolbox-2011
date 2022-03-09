function y = icqt(Xcqt,Rs,Ra,cut)
%y = icqt(Xcqt) computes the inverse CQT of the CQT coefficients in Xcqt.spCQT
%
%The input structue Xcqt is the structure gained by cqt() and cqtPerfectRast(), respectively. 
%If the CQT coefficients in Xcqt.spCQT are not changed, the output y is the
%reconstructed (near-perfect) time-domain signal of the input signal x
%(cqt(x,...)) withing the frequency range [fmin fmax].
%
%INPUT:
%   Xcqt        ... structure obtained by the function cqt(..)
%   Rs          ... synthesis hop size
%   Ra          ... analysis hop size
%   cut         ... if set to 1, the output signal will be cropped to the
%                   length of the input signal 
%
%OUTPUT:
%   y           ... reconstructed time domain signal
%
%Christian Schörkhuber, Anssi Klapuri 2010-06
%2011-08: allow for synthesis hop size different from analysis hop size

cellCQT = sparse2cell(Xcqt.spCQT,Xcqt.bins,Xcqt.octaveNr,Xcqt.intParams.atomNr);

FFTLen = Xcqt.intParams.fftLEN;
octaveNr = Xcqt.octaveNr;

if Rs==Ra, HOPSZ = Xcqt.intParams.fftHOP;
else HOPSZ = Rs;
end


%% Kernel for inverse transform
Kinv = Xcqt.fKernel;

%% inverse transform
y = [];
for i = octaveNr:-1:1
    cellCQT_oct = cellCQT{i};
    emptyHops = Xcqt.intParams.firstcenter/Xcqt.intParams.atomHOP;
    dropped = emptyHops*2^(octaveNr-i)-emptyHops;
    fill = dropped * Ra; %to synchronize octaves
    
    Y = Kinv * cellCQT_oct; %compute spectrum of reconstructed signal for all coefficients in this octave  
    y_oct_temp = ifft(Y);
    y_oct = 2*real(y_oct_temp); %Y contains no negative frequencies -> keep only real part*2 to 
                                %reconstruct real valued time signal 
    NBLOCKS = size(Y,2);      
    siglen = FFTLen + (NBLOCKS-1)*HOPSZ;
    y = [y;zeros(siglen-length(y),1)];
    for n = 1:NBLOCKS
        y((n-1)*HOPSZ+1+fill:((n-1)*HOPSZ)+FFTLen+fill) = y_oct(:,n) + y((n-1)*HOPSZ+1+fill:((n-1)*HOPSZ)+FFTLen+fill); %overlap-add
    end
    
    if(i~=1) %upsampling by factor two
         y = upsample(y,2); %insert one zero between each sample
         y = filtfilt(Xcqt.intParams.filtCoeffB,Xcqt.intParams.filtCoeffA,y);
         y = y * 2;
    end

end

if cut
    y = y(Xcqt.intParams.preZeros+1:end); %crop introduced zeros at the beginning
    y = y(1:Xcqt.intParams.xlen_init); %crop overhead zeros at the end
end
