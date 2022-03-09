

%% settings
shiftBins = +24; % desired pitch transpostion in CQT bins

SEPERATE = 0; % set to one to apply harmonic/percussive seperation prior to pitch transposition
METHOD = 1; %method 1 ... without estimation of instantaneous frequencies 
            %method 2 ... with estimation of instantaneous frequencies
PLOT = 0; %if set to 1, the rasterized CQT representation of the input signal will be plotted

%% input signal
[x fs] = wavread('kempff1.wav');

%% parameters
bins = 48; %CQT resolution [bins/octave]
fmax = fs/2; %highest frequency analyzed
fmin = fmax/2^9; %lowest frequency of interest (CQT bins will start immediatly above fmin)

%% resample (in order to anaylize signal up to Nyquist frequency)
x = resample(x,4,3);
fs = fs*4/3;

%% CQT
atomHopFactor = 0.25 /stretch;
if stretch == 1, allowSevAtoms = 1; else allowSevAtoms = 0; end;
Xcq = cqt(x,fmin,fmax,bins,fs,'atomHopFactor',atomHopFactor,'allowSevAtoms',allowSevAtoms,'oversampTwo',1);
X = Xcq.spCQT;
ahop = Xcq.intParams.atomHOP;

%% plotting
if PLOT
   fcomp = 0.7;
   method = 'surf';
   plotCQT(Xcq,fs,fcomp,method); 
end

%% harmonic/percussive separation
if SEPERATE
    lh = 15;
    lp = 25;
    a = 2;
    [H P] = separationPercHarm(Xcq,lh,lp,5,a);
    Xcq_h = Xcq;
    Xcq_h.spCQT = H;
    Xcq_p = Xcq;
    Xcq_p.spCQT = P;
    Xt = Xcq_h;
else
    Xt = Xcq; 
end

%% phase update
factor = 2^(shiftBins/bins);
Y = phaseUpdate(Xt,factor,fs,METHOD);

%% shift
if shiftBins ~=0
    Y = circshift(Y,shiftBins);
    if shiftBins > 0
        Y(1:shiftBins) = 0;
    elseif shiftBins < 0
        Y(end-shiftBins+1:end) = 0;
    end
end

if SEPERATE, Y = Y + P; end;

%% ICQT
shop = ahop;
Xcq.spCQT = Y;

y = icqt(Xcq,shop,ahop,1);

%% resample to original sampling rate
y = resample(y,3,4); y = y/max(abs(y));
x = resample(x,3,4); x = x/max(abs(x));
fs = fs*3/4;
