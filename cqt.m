function Xcqt = cqt(x,fmin,fmax,bins,fs,varargin) 

%Xcqt = cqt(x,fmin,fmax,bins,fs,varargin)
%
%Computes the constant-Q transform of the input signal x.
%
%INPUT:
%   fmin      ... lowest frequency of interest
%   fmax      ... highest frequency of interest
%   bins      ... frequency bins per octave
%   fs        ... sampling rate
%
%   optional input parameters (parameter name/value pairs):
%
%   'atomHopFactor'   ... overlap of temporal atoms in percent. Default: 0.25.
%    
%   'q'             ... the maximum value for optimal reconstruction is q=1.
%                       For values smaller than 1 the bandwidths of the spectral
%                       atoms (filter) are increased retaining their center
%                       frequencies (frequency 'smearing', frequency domain redundancy 
%                       increases, time resolutin improves). Default: 1.
%   'thresh'        ... all values in the cqt kernel smaller than tresh are
%                       rounded to zero. A high value for thresh yields a
%                       very sparse kernel (fast) but introduces a bigger error. 
%                       The default value is chosen so that the error due to rounding is negligible.
%   'kernel'        ... if the cqt kernel structure has been precomputed
%                       (using function 'genCQTkernel'), the computation of the kernel
%                       will be by-passed below).
%   'oversampTwo'   ... if set to 1, time-domain oversampling of factor 2 is applied. 
%                       That is, twice as many sampling points are computed (in order to
%                       allow for interoctave phase adjustments)
%   'allowSevAtoms' ... if set to 1, the generated kernel may contain
%                       more than one atom per FFT frame thus increasing
%                       the performance of the CQT. If time-scaling using
%                       the CQT phase vocoder method is about to be
%                       applied, allowSevAtoms has to be set to 0 as the
%                       hop size will be altered during resynthesis.
%                       Default: 1.
%   'coeffB',
%   'coeffA'        ... Filter coefficients for the anti-aliasing filter, where
%                       'coeffB' is the numerator and 'coeffA' is the
%                       denominator (listed in descending powers of z). 
%                                                  
%OUTPUT:
%   Xcqt      ... struct that comprises various fields: 
%              spCQT: CQT coefficients in the form of a sparse matrix 
%                    (rasterized, not interpolated)
%              fKernel: spectral Kernel 
%              fmin: frequency of the lowest bin
%              fmax: frequency of the hiqhest bin
%              octaveNr: number of octaves processed
%              bins: number of bins per octave
%              intParams: structure containing additional parameters for the inverse transform   
%
%Christian Schörkhuber, Anssi Klapuri 2010-06
%2011-03: error removed in the calculation of fmin for output in Xcqt structure
%2011-08: - 'allowSevAtoms' parameter added, 
%         - oversampling of factor 2 implemented ('oversampTwo'),
%         - window function option removed: the used window function now is modhann(N)
%         - phase of atoms altered to ensure predictible vertical phase relations

%% input checking
if size(x,2) > 1 && size(x,1) > 1, error('cqt requires one-dimensional input!'); end;
if size(x,2) > 1, x = x(:); end; %column vector

%% input parameters
q = 1; %default value
atomHopFactor = 0.25; %default value
thresh = 0.0005; %default value
allowSevAtoms = 1;
oversampTwo = 0;
perfRast = 0;

for ain = 1:1:length(varargin)
    if strcmp(varargin{ain},'q'), q = varargin{ain+1}; end;
    if strcmp(varargin{ain},'atomHopFactor'), atomHopFactor = varargin{ain+1}; end;
    if strcmp(varargin{ain},'thresh'), thresh = varargin{ain+1}; end;
    if strcmp(varargin{ain},'allowSevAtoms'), allowSevAtoms = varargin{ain+1}; end;
    if strcmp(varargin{ain},'oversampTwo'), oversampTwo = varargin{ain+1}; end;
    if strcmp(varargin{ain},'kernel'), cqtKernel = varargin{ain+1}; end;
    if strcmp(varargin{ain},'coeffB'), B = varargin{ain+1}; end;
    if strcmp(varargin{ain},'coeffA'), A = varargin{ain+1}; end;    
end

%% define
octaveNr = ceil(log2(fmax/fmin));
fmin = (fmax/2^octaveNr) * 2^(1/bins); %set fmin to actual value
xlen_init = length(x);

%% design lowpass filter
if ~exist('B','var') || ~exist('A','var')
    LPorder = 6; %order of the anti-aliasing filter
    cutoff = 0.5;
    [B A] = butter(LPorder,cutoff,'low'); %design f_nyquist/2-lowpass filter
end

%% design kernel for one octave 
if ~exist('cqtKernel','var')
    cqtKernel = genCQTkernel(fmax, bins,fs,'q',q,'atomHopFactor',atomHopFactor, ...
        'thresh',thresh,'allowSevAtoms',allowSevAtoms,'perfRast',perfRast,'oversampTwo',oversampTwo);
end

%% calculate CQT
maxBlock = cqtKernel.fftLEN * 2^(octaveNr-1); %largest FFT Block (virtual)
suffixZeros = maxBlock; 
prefixZeros = maxBlock;
x = [zeros(prefixZeros,1); x; zeros(suffixZeros,1)]; %zeropadding
OVRLP = cqtKernel.fftLEN - cqtKernel.fftHOP;

K = cqtKernel.fKernel'; %conjugate spectral kernel for cqt transformation

%shifted kernel (for oversampling by factor 2)----
shift = cqtKernel.atomHOP/2;
phShiftVec = exp(1i*2*pi.*(0:(cqtKernel.fftLEN-1))'*shift./cqtKernel.fftLEN);
phShiftMat = repmat(phShiftVec.',size(K,1),1);
K2 = K .* phShiftMat;
%-------------------------------------------------

emptyHops = cqtKernel.firstcenter/cqtKernel.atomHOP;
fftBlockNr = ceil((length(x)-cqtKernel.fftLEN) / (cqtKernel.fftLEN-OVRLP))+1;
tFrameNr = fftBlockNr * cqtKernel.atomNr;
spCQT = zeros(bins*octaveNr,tFrameNr);

for i = 1:octaveNr
    binVec = bins*(octaveNr-i)+1:bins*(octaveNr-i+1);
    drop = emptyHops*2^(octaveNr-i)-emptyHops; %first coefficients of all octaves have to be in synchrony
    xx = buffer(x,cqtKernel.fftLEN, OVRLP,'nodelay'); %generating FFT blocks
    XX = fft(xx); %applying fft to each column (each FFT frame)
    
    Xoct = K*XX; %calculating cqt coefficients for all FFT frames for this octave
    
    %reshape if necessary-----
    if cqtKernel.atomNr > 1
        XoctReSh = zeros(bins,size(Xoct,2)*cqtKernel.atomNr-drop);
        for m = 0:bins-1
           Xtemp = Xoct(m*cqtKernel.atomNr+1:(m+1)*cqtKernel.atomNr,:);
           Xtemp2 = reshape(Xtemp,1,size(Xoct,2)*cqtKernel.atomNr);
           XoctReSh(m+1,:) = Xtemp2(1+drop:end);
        end
        Xoct = XoctReSh;
    else
        Xoct = Xoct(:,1+drop:end);
    end
    %-------------------------
    
    tVec = 1:2^(i-1):size(Xoct,2)*2^(i-1); 
    spCQT(binVec,tVec) = Xoct;
    
    if oversampTwo == 1
        if i > 1 %oversampling by factor 2 -> compute coefficients with shifted kernel
            Xoct2 = K2*XX; 

            %reshape if necessary-----
            if cqtKernel.atomNr > 1
                XoctReSh = zeros(bins,size(Xoct2,2)*cqtKernel.atomNr-drop);
                for m = 0:bins-1
                    Xtemp = Xoct2(m*cqtKernel.atomNr+1:(m+1)*cqtKernel.atomNr,:);
                    Xtemp2 = reshape(Xtemp,1,size(Xoct2,2)*cqtKernel.atomNr);
                    XoctReSh(m+1,:) = Xtemp2(1+drop:end);
                end
                Xoct2 = XoctReSh;
            else
                Xoct2 = Xoct2(:,1+drop:end);
            end
            %-------------------------

            tVec = tVec + 2^(i-2); 
            spCQT(binVec,tVec) = Xoct2;    
        end
    end
    
    if i~=octaveNr
        x = filtfilt(B,A,x); %anti aliasing filter
        x = x(1:2:end); %drop samplerate by 2 (not x(1:2:end), this would cause a phase shift between octaves!)
    end
end
spCQT = sparse(spCQT);

%% map to sparse matrix representation
% spCQT = cell2sparse(cellCQT,octaveNr,bins,cqtKernel.firstcenter,cqtKernel.atomHOP,cqtKernel.atomNr);
% this functionality is now applied during the computation of the
% coefficients (

%% return
intParam = struct('sufZeros',suffixZeros,'preZeros',prefixZeros,'xlen_init',xlen_init,'fftLEN',cqtKernel.fftLEN,'fftHOP',cqtKernel.fftHOP,...
    'q',q,'filtCoeffA',A,'filtCoeffB',B,'firstcenter',cqtKernel.firstcenter,'atomHOP',cqtKernel.atomHOP,...
    'atomNr',cqtKernel.atomNr,'Nk_max',cqtKernel.Nk_max,'Q',cqtKernel.Q,'rast',0,'oversampTwo',oversampTwo);

Xcqt = struct('spCQT',spCQT,'fKernel',cqtKernel.fKernel,'fmax',fmax,'fmin',fmin,'octaveNr',octaveNr,'bins',cqtKernel.bins,'intParams',intParam);



