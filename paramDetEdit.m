function [rms,peakFr,bw10db,bw3db,F0,dur] = ...
    paramDetEdit(MSN,posClick,rmsClick,fs,tffreq,tfuppc)

%Take timeseries out of existing file, convert from normalized data to
%counts
%1) calculate spectral received levels RL for click and preceding noise: 
%calculate spectra, account for bin width to reach dB re counts^2/Hz, 
%add transfer function, compute peak frequency and bandwidth
%2) calculate RLpp at peak frequency: find min & max value of timeseries,
%convert to dB, add transfer function value of peak frequency (should come
%out to be about 9dB lower than value of spectra at peak frequency)
%3) calculate signal to noise ratio and eliminate clicks that do not have a
%signal to noise ratio larger than 10dB at peak frequency.

%sb 090814
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) calculate spectral received levels RL: calculate spectra
%compute peak frequency
N = fs/1000;
click=zeros(size(MSN,1),N/2+1);

BinHz = 0:1:N/2;

%analysis of one click after the other
for s=1:size(MSN,1)
    %calculate duration of full click as determined by teager energy in us
    %     dur(s) = (posClick(s,2)-posClick(s,1))/fs*1000*1000;
    dur(s) = (rmsClick(s,2)-rmsClick(s,1))/fs*1000*1000;
    
    
    %pick click as defined by teager energy duration + 30 pts on each side
    posStart=posClick(s,1)-30;
    posEnd=posClick(s,2)+30;
    
    if posStart <1
        posStart = 1;
    end
    if (posEnd > N+2)
        posEnd = N+2;
    end
    
    if ~isnan(posStart)
        x = MSN(s,posStart:posEnd);
       
        wind=hann(posEnd-posStart+1);
        wClick=x.*wind.';
        
        %zero pad to length of N
        nZero = N-length(x);
        pad = zeros(1,nZero);
        wClick = [wClick pad];
        
        spClick=20*log10(abs(fft(wClick,N)));
        spClick = spClick(1:N/2+1);
        click(s,:)=spClick;

    else
        click(s,:) = ones(1,size(click,2))*NaN;
    end

end

tf = interp1(tffreq, tfuppc, BinHz, 'linear', 'extrap');

%add transfer function
for i=1:size(click,1)
    specClickTf(i,:)=click(i,:)+tf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute parameters of click
peakFr=zeros(size(MSN,1),1);
bw10db=zeros(size(MSN,1),3);
bw3db=zeros(size(MSN,1),3);
F0=zeros(size(MSN,1),1);
rms=zeros(size(MSN,1),1);

%analysis of one click after the other
for n=1:size(specClickTf,1) 
    click=specClickTf(n,:);
    if ~isnan(click)
        %calculate peak frequency
        %max value in the first half samples of the spectrogram within
        %80-100 kHz
        cut = 10;
        fmin = BinHz-cut;
        [fCut posCut]= min(abs(fmin));

        [valMx(n) posMx]=max(click(posCut:end));
        posMx = posMx + posCut - 1;

        %calculation from spectrogram -> from 0 to 100kHz in 256 steps (FFT=512)
        peakFr(n,1)=BinHz(posMx); %peak frequency in kHz

        %calculate rms
        rmsStart = rmsClick(n,1);
        rmsEnd = rmsClick(n,2);
        y2 = MSN(n,rmsStart:rmsEnd).^2;
        rms(n,1) = 10 * log10(sum(y2)/(rmsEnd-rmsStart));
        rms(n,1) = rms(n,1) + tf(peakFr(n,1) + 1);
        %%%%%
        %calculate bandwidth

        %-3dB bandwidth   
        %calculation of -3dB bandwidth - amplitude associated with the halfpower points of a pressure pulse (see Au 1993, p.118); 
        low=valMx(n)-3; %p1/2power = 10log(p^2max/2) = 20log(pmax)-3dB = 0.707*pmax; 1/10^(3/20)=0.707
        %walk along spectrogram until low is reached on either side
        slopeup=fliplr(click(1:posMx));
        slopedown=click(posMx:round(length(click)));
        for e3dB=1:length(slopeup)
           if slopeup(e3dB)<low %stop at value < -3dB: point of lowest frequency
               break
           end
        end
        for o3dB=1:length(slopedown)
           if slopedown(o3dB)<low %stop at value < -3dB: point of highest frequency
               break
           end
        end

        %-10dB bandwidth
        low=valMx(n)-10;
        %walk along spectrogram until low is reached on either side
        slopeup=fliplr(click(1:posMx));
        slopedown=click(posMx:end);
        for e10dB=1:length(slopeup)
           if slopeup(e10dB)<low %stop at value < -10dB: point of lowest frequency
               break
           end
        end
        for o10dB=1:length(slopedown)
           if slopedown(o10dB)<low %stop at value < -10dB: point of highest frequency
               break
           end
        end
        %JAH CORRECTION FOR BUG +1 and -1 8/2016
        %calculation from spectrogram -> from 0 to 100kHz in 256 steps (FFT=512)
        high3dB=BinHz(posMx+o3dB-1); %-3dB highest frequency in kHz
        low3dB=BinHz(posMx-e3dB+1); %-3dB lowest frequency in kHz
        high10dB=BinHz(posMx+o10dB-1); %-10dB highest frequency in kHz
        low10dB=BinHz(posMx-e10dB+1); %-10dB lowest frequency in kHz
        bw3=high3dB-low3dB;
        bw10=high10dB-low10dB;

        bw3db(n,1)=low3dB;
        bw3db(n,2)=high3dB;
        bw3db(n,3)=bw3;

        bw10db(n,1)=low10dB;
        bw10db(n,2)=high10dB;
        bw10db(n,3)=bw10;

        %%%%%
        %calculate further click parameters

        %frequency centroid (or center frequency) in kHz
        linearSpec=10.^(click/20);
        Freq_vec = BinHz(cut:end);
        F0(n)=(sum(Freq_vec.*linearSpec(cut:length(linearSpec)).^2)/sum(linearSpec(cut:length(linearSpec)).^2)); % Au 1993, equation 10-3
    else
        peakFr(n,1)=NaN;
        bw10db(n,:)=ones(1,3)*NaN;
        bw3db(n,:)=ones(1,3)*NaN;
        F0(n,1)=NaN;
    end
end


           