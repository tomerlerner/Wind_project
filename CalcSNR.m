function [SNRin,SNRout,NRout,SNRAin,SNRAout,NRAout,SNRsegin,...
SNRsegout,NRsegout,SNRAsegin,SNRAsegout,NRAsegout]=...
CalcSNR(noise,speech,spest,fs,seglength,Tsp,Tno)

% [SNRin,SNRout,NRout,SNRAin,SNRAout,NRAout,SNRsegin,...
% 
% SNRsegout,NRsegout,SNRAsegin,SNRAsegout,NRAsegout]=...
% 
% CalcSNR[noise,speech,spest,fs,seglength,Tsp,Tno]
% 
% 23
% 
% CalcSNR calculates different variations of SNR and Noise Residual
% 
% measures.
% 
% Optional arguments are seglengthi=100), Tsp(=0.05] and Tno(=0.05].
% 
% -°s
% 
% Input:
% 
% noise - noise that speech is mixed with to produce input to algorithm
% 
% speech - speech that noise is mixed with to produce input to algorithm
% 
% spest - speech estimate output from algorithm
% 
% fs - sample rate
% 
% seglength - length of a segment used to calculate power in ms
% 
% Tsp - Threshold of speech power for a segment in pct of mean
% 
% Tno - Threshold of noise power for a segment in pct of mean
% 
% 25
% 
% Output:
% 
% SNRin - Signal to noise ratio of the signal used as input to algorithm
% 
% SNRout - Signal to noise ratio of signal output from algorithm
% 
% NRout - Noise Residual of output from algorithm
% 
% outputs with 'A' in the name is the same as corresponding output
% 
% without 'A', but A-weighted.
% 
% outputs with 'seg' in the name is the same as corresponding output
% 
% without 'seg', but calculated for each segment.
% 
% 25
% 
% Written by: Kristian Timm Andersen, IHH DTU 2007

%Check for input

if ((nargin<5)||isempty(seglength)),seglength=100;end

if ((nargin<6)||isempty(Tsp)),Tsp=0.05;end

if ((nargin<7)||isempty(Tno)),Tno=0.05;end

Lout=length(spest); % Length of output

outnosp=spest-speech; % output minus speech

outnono=spest-noise; % output minus noise

seglength=floor(fs*seglength/1000); % number of sample points in a segment

segno=floor(Lout/seglength); % Number of segments

% Number of FFT points used in welch's method

NFFT = max(256,2^ceil(log2(seglength/4.51)));
Nwelch=NFFT/2+1; % Number of elements in output of Uelch's method
% Corresponding frequency of each output element in Uelch's method

f=(0:Nwelch-1)*fs/NFFT;

Aw=Aweight(f)'; % Calculate A-weight for each frequency

% Initialise power matrices. Each element contains the power of a segment

sp=zeros(1,segno);spw=zeros(1,segno);no=zeros(1,segno);now=zeros(1,segno);

outw=zeros(1,segno);nospw=zeros(1,segno);nonow=zeros(1,segno);

spseg=zeros(seglength,segno);noseg=zeros(seglength,segno);

outnospseg=zeros(seglength,segno);outnonoseg=zeros(seglength,segno);

%calculate the power for each segment using Melch's method and A-weight

for i=1:segno
spseg(:,i)=speech(1+(i-1)*seglength:i*seglength);
noseg(:,i)=noise(1+(i-1)*seglength:i*seglength);
outnospseg(:,i)=outnosp(1+(i-1)*seglength:i*seglength);
outnonoseg(:,i)=outnono(1+(i-1)*seglength:i*seglength);
sp(i)=sum(pwelch(spseg(:,i)));
spw(i)=sum(pwelch(spseg(:,i)).*Aw);
no(i)=sum(pwelch(noseg(:,i)));
now(i)=sum(pwelch(noseg(:,i)).*Aw);
outw(i)=sum(pwelch(spest(1+(i-1)*seglength:i*seglength)).*Aw);
nospw(i)=sum(pwelch(outnospseg(:,i)).*Aw);
nonow(i)=sum(pwelch(outnonoseg(:,i)).*Aw);

end

Tsp=Tsp*mean(Sp); % Calculate threshold for speech

Tno=Tno*mean(no); % Calculate threshold for noise

% Find intersection of speech and noise segments that are larger than

% threshold

evalframes=logical((sp>Tsp).*(no>Tno));

SNRin=10*log10(sum(speech.^2)/sum(noise.^2)); % Input SNR

SNRout=10*log10(sum(speech.^2)/sum(outnosp.^2)); % Output SNR

NRout=10*log10(sum(noise.^2)/sum(outnono.^2)); % Noise Residual output

SNRAin=10*log10(sum(spwtevalframes))/sum(now(evalframes));

SNRAout=10*log10(sum(spw(evalframes))/sum(nospw(evalframes)));

NRAout=10*log10(sum(now(evalframes))/sum(nonow(evalframes)));

SNRsegin=10*log10(sum(spseg(:,evalframes).^2,1)...
./sum(noseg(:,evalframes).^2,1));

SNRsegout=10*log10(sum(spseg(:,evalframes).^2,1)...
./sum(outnospseg(:,evalframes).^2,1));

NRsegout=10*log10(sum(noseg(:,evalframes).^2,1)...
./sum(outnonoseg(:,evalframes).^2,1));

SNRAsegin=10*log10(spw(evalframes)./now(evalframes));

SNRAsegout=10*log10(spw(evalframes)./nospw(evalframes));

NRAsegout=10*log10(now(evalframes)./nonow(evalframes));