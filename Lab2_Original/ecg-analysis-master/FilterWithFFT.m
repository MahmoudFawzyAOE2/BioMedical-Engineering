function filtered = FilterWithFFT(ecg)
% 1. Splits the input signal into 25000 long segments.
% 2. Transforms the segment into frequency domain with FFT.
% 3. Applies band-pass filtering.
% 4. Transforms the segment back to time domain with IFFT.
% 5. Merges the 25000 long segments back into the original input size.

segment = 25000;       
part = zeros(segment,1);
filtered = zeros(numel(ecg),1);
numOfSegments = (numel(ecg)) / segment;
for i = 1:numOfSegments
    lowerBound = ((i-1) * segment)+1;
    part = ecg(lowerBound:i*segment);
    part = filterfft(part);
    filtered(lowerBound:i*segment) = part;
end
end

function filtered = filterfft(y)
% FFT -> Band-pass -> IFFT
    centerfrequency = 1803;
    filterwidth = 2720;
    filtershape = 57 ;

    y=reshape(y,1,length(y));
    fy=fft(y);
    lft1=[1:(length(fy)/2)];
    lft2=[(length(fy)/2+1):length(fy)];
    % Compute filter shape.
        ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
        ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
        ffilter=[ffilter1,ffilter2];

    if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
    ffy=fy.*ffilter;  % Multiply filter by Fourier Transform of signal
    filtered=real(ifft(ffy));
    filtered=reshape(filtered,length(filtered),1);

end

function g = shape(x,pos,wid,n)
% Determines the shape of the band-pass filter
if n==0
    g=ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
else
    g = exp(-((x-pos)./(0.6.*wid)) .^(2*round(n)));
end
end


%subplot(2,1,1)
%plot(data);
%subplot(2,1,2);
%plot(ry);	
