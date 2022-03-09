function Y = phaseUpdate(Xcq, factor, fs, method)
%Y = phaseUpdate(Xcq, factor, fs, method)
%
%Alters the phases of CQT coefficients in order to retain vertical and
%horizontal phase coherence for later time-scaling (CQT phase vocoder) or
%pitch-scaling (transposition of CQT coefficients)
%
%INPUT:
%   Xcq          ... structure returned by function 'cqt'
%   factor       ... time-/pitch-scaling factor (for the phase update
%                    there is no differnce between time- and pitch-scaling)
%   fs           ... sampling frequency
%   method       ... with (method=2)/without(method=1) estimation of instantaneous frequencies 
%
%OUTPUT:
%   Y   ... sparse representation of CQT coefficients with updated phase values 
%
%Christian Schörkhuber, 2011-08

bins = Xcq.bins;
fmin = Xcq.fmin;
X = Xcq.spCQT;
ahop = Xcq.intParams.atomHOP;



%% phase update (without f-estimation)--------------------------------------------------
if method == 1 %method = 1 -> without f-estimation

    Y = zeros(size(X));
    Y(:,1) = X(:,1);
    BVec = 1:(Xcq.octaveNr*bins);
    octs = Xcq.octaveNr;
    accumulated_rotation_angles=zeros(1,size(X,1));
    for n = 2:size(X,2)

        %--how many octaves do have coefficients for this timeslice?---
        tv = (n-1)./2.^(1:octs-2); 
        tv = tv - floor(tv);
        tv = tv(tv==0);
        lowerOctsNr = length(tv);
        binVec = BVec(end-(2+lowerOctsNr)*bins+1:end);
        frame = X(binVec,n);
        %highest octave and 2nd to highest octave always have coefficitens
        %(oversampling 2). lowerOctsNr is number of how many octaves also have
        %coefficients in this frame
        %--------------------------------------------------------------

        peaks = pickPeaksSparse(frame,0);
        peaks = binVec(peaks);
        if ~isempty(peaks)
            % Find regions of influence around peaks
            regions = round(0.5*(peaks(1:end-1)+peaks(2:end)));  
            regions = [binVec(1), regions, binVec(end)];

            modFrame = zeros(size(X,1),1);

            for u=1:length(peaks)
                thisBin = peaks(u);
                tempOct = ceil(thisBin/bins);
                nJump = 1;
                if tempOct <= octs-1, nJump = 2^(octs-tempOct-1); end;
                thisAHop = ahop * nJump;

                % Compute the rotation angle required, which has 
                % to be cumulated from frame to frame
                shift_freq = fmin*2^(thisBin/bins)*(factor-1);
                norm_shift_freq = 2*pi*shift_freq/fs;
                rotation_angles=accumulated_rotation_angles(thisBin) + thisAHop*norm_shift_freq;

                % Overlap/add the bins around the peak, changing the
                % phases accordingly
                modFrame(regions(u):regions(u+1))= ...
                    modFrame(regions(u):regions(u+1)) + X(regions(u):regions(u+1),n) * exp(1i*rotation_angles);
                accumulated_rotation_angles((regions(u):regions(u+1))) = rotation_angles;   
            end    
       else
           modFrame = X(:,n); %if no peaks are found
       end
       Y(:,n) = modFrame; 
    end

end %if method

%% phase update (with f-estimation)--------------------------------------------------
if method == 2 %method = 2 -> with f-estimation

    Y = zeros(size(X));
    phi = angle(X);
    Y(:,1) = X(:,1);
    BVec = 1:(Xcq.octaveNr*bins);
    FVec = fmin * 2.^(BVec/bins); %bin frequencies
    octs = Xcq.octaveNr;
    Omega_k = 2*pi*FVec/fs; %bin frequencies (normalized) 
    Omega_k = Omega_k';
    accumulated_rotation_angles=zeros(1,size(X,1));
    for n = 2:size(X,2)

        %--how many octaves do have coefficients for this timeslice?---
        tv = (n-1)./2.^(1:octs-2); 
        tv = tv - floor(tv);
        tv = tv(tv==0);
        lowerOctsNr = length(tv);
        binVec = BVec(end-(2+lowerOctsNr)*bins+1:end);
        frame = X(binVec,n);
        %highest octave and 2nd to highest octave always have coefficitens
        %(oversampling 2). lowerOctsNr is number of how many octaves also have
        %coefficients in this frame
        %--------------------------------------------------------------

        peaks = pickPeaksSparse(frame,0);
        peaks = binVec(peaks);
        if ~isempty(peaks)
            % Find regions of influence around peaks
            regions = round(0.5*(peaks(1:end-1)+peaks(2:end)));  
            regions = [binVec(1), regions, binVec(end)];

            modFrame = zeros(size(X,1),1);

            for u=1:length(peaks)
                thisBin = peaks(u);
                tempOct = ceil(thisBin/bins);
                nJump = 1;
                if tempOct <= octs-1, nJump = 2^(octs-tempOct-1); end;
                thisAHop = ahop * nJump;

                %estimate true frequency of peak---------------
                dphi = phi(peaks(u),n) - phi(peaks(u),n-nJump) - thisAHop*Omega_k(peaks(u));   
                dphi_p = wrapToPi(dphi); %principle phase difference; 
                omega_k_est = Omega_k(peaks(u)) + dphi_p/thisAHop; %instantaneous frequency (normalized)
                fk_est = omega_k_est * fs / (2*pi);
                %----------------------------------------------

                %delta = thisSHop*omega_k_est - thisAHop*omega_k_est;
                %shift_freq = fk_est*(2^(shiftBins/bins) - 1);
                shift_freq = fk_est*(factor-1);
                norm_shift_freq = 2*pi*shift_freq/fs;
                rotation_angles=accumulated_rotation_angles(peaks(u)) + thisAHop*norm_shift_freq;
                %-------------------------------------------------------------

                % Overlap/add the bins around the peak, changing the
                % phases accordingly
                modFrame(regions(u):regions(u+1))= ...
                    modFrame(regions(u):regions(u+1)) + X(regions(u):regions(u+1),n) * exp(1i*rotation_angles);
                accumulated_rotation_angles((regions(u):regions(u+1))) = rotation_angles;   
            end    
       else
           modFrame = X(:,n); %if no peaks are found
       end
       Y(:,n) = modFrame; 
    end

end %if method

