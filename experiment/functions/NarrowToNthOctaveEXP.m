function [sortedData,Fc,Flow,Fhigh] = NarrowToNthOctaveEXP(arrayOfNBCenterFreqs,arrayOfdBToConvert,n)

df=arrayOfNBCenterFreqs(2)-arrayOfNBCenterFreqs(1);

% %Use correctionFactor=1 if you did not window or you used ECF
correctionFactor=1;
% %Use this if you used a hanning window with ACF
% ACF=(1/mean(window(@hann,2^11)));
% ECF=(1/mean(window(@hann,2^11)));
% correctionFactor=1/ACF;
arrayOfdBToConvert=arrayOfdBToConvert*correctionFactor;

%Determine Initial Center frequency (Fc)
previous=1000*2^(1/n);
lowFc=1000; %Start at 1000 Hz and find lowest Fc
while (previous-lowFc) > df
    previous=lowFc;
    lowFc=lowFc/2^(1/n);
end

%Compute center frequencies
ii=1;
Fc(ii)=lowFc;
while Fc(ii) <= max(arrayOfNBCenterFreqs) 
    ii=ii+1;
    Fc(ii)=2^(1/n)*Fc(ii-1);
end

%Compute high, low frequencies, and bandwidth (BW) from center frequencies
    Flow=Fc/sqrt(2^(1/n));
    Fhigh=Fc*sqrt(2^(1/n));
    BW=Fhigh-Flow;

%Sort data in frequency bins
sortedData=zeros(length(Fc),size(arrayOfdBToConvert,2)); zeroFreqValue=0;
for jj=1:size(arrayOfdBToConvert,1) %Do for each value given
    for kk=1:length(Fc) %Check each Fc bin to see where value should be placed
        if arrayOfNBCenterFreqs(jj)>=Flow(kk) && arrayOfNBCenterFreqs(jj)<Fhigh(kk) %Find place. In the rare case a given center freq = flow(ii)/fhigh(ii+1) then sum the value in flow
            if sortedData(kk,:)==0 %if no values has been added to the band then set initial value
                sortedData(kk,:)=arrayOfdBToConvert(jj,:);
                N(kk)=1;
            else %else sum in the value
                sortedData(kk,:)=sqrt(sum(([sortedData(kk,:); arrayOfdBToConvert(jj,:)].^2))); %dB add value to that bin
                N(kk)=N(kk)+1;
            end
        else %else check to see if the value belongs in the 0 Hz band
            if (kk==1) && arrayOfNBCenterFreqs(jj)<Flow(kk)
                if zeroFreqValue==0
                    zeroFreqValue=arrayOfdBToConvert(jj,:);
                else
                    zeroFreqValue=sqrt(sum(([zeroFreqValue; arrayOfdBToConvert(jj,:)].^2))); %dB add value to that bin
                end
            end
        end
    end
end

% Apply bin correction factor for when center frequencies are not integer
% multiples of the bandwidth (BW)
for ll=1:length(N) %Look through each bin that had values added to it
    if (N(ll)~=0) && (floor(BW(ll)/df)<=N(ll)) %If values were added and there are enough values added to make applying the correction appropriate then apply correction
        sortedData(ll,:)= sortedData(ll,:)/(sqrt((df*N(ll))/BW(ll)));
    end
end

%If DC value was given, add to sorted array
if zeroFreqValue~=0
    sortedData=[zeroFreqValue; sortedData];
    Fc=[0 Fc];
    Flow=[0 Flow];
    Fhigh=[0 Fhigh];
end

%Remove unused bins at beginning if only higher freqs were given
removeThisMany=0;
for ll=1:size(sortedData,1)
    if sortedData(ll,:)==0
        removeThisMany=1+removeThisMany;
    else
        break
    end
end
sortedData=sortedData(removeThisMany+1:end,:);
Fc=Fc(removeThisMany+1:end);
Flow=Flow(removeThisMany+1:end);
Fhigh=Fhigh(removeThisMany+1:end);

%Remove bin at end if unused
if sortedData(end)==0
    sortedData=sortedData(1:end-1,:);
    Fc=Fc(1:end-1);
    Flow=Flow(1:end-1);
    Fhigh=Fhigh(1:end-1);
end