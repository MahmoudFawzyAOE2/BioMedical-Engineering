function [figMain, figDiffLine, figDiffHist, diffLine, locs, BPM, fwhm, outliers] = PeakDetect2(ecg, fs, minPeakHeight, minPeakDist)
%PeakDetect R-peak detection, with plots
%   INPUTS:
%       inputPath - A full path to a .mat file.
%       fs - Frequency, given in Hertz. Mandatory input argument.
%       minPeakHeight - The minimum peak height to detect a peak.
%           Optional argument, (default = 0.5).
%       minPeakDist - The minimum distance between 2 peaks.
%           Optional argument, (default = 50).
%
%   OUTPUTS:
%       figMain, figDiffLine, figDiffHist - references for the figures.
%           Function AnalyzeECG() uses them.
%       diffLine - The difference line diagram,, containing the outliers
%       locs - locations of the peaks
%       BPM - Beat per Minute
%       FWHM - Full-Width at Half-Maximum
%       outliers - the number of outliers found.

%Example
%[figMain, figDiffLine, figDiffHist, ecg, diffMap, locs, BPM, fwhm, outliers] = PeakDetect("001.mat",1000);
%[figMain, figDiffLine, figDiffHist, ecg, diffMap, locs, BPM, fwhm, outliers] = PeakDetect("001.mat",1000,0.5,50);

close all; clc;

%Argument check
% if exist(inputPath,'file') ~= 2
%     error('File not found! %s',inputPath)
% end
if nargin < 2 || isempty(fs)
    error('Input argument fs is required.')
end
if nargin < 3 || isempty(minPeakHeight)
    minPeakHeight = 0.5;
end
if nargin < 4 || isempty(minPeakDist)
    minPeakDist = 50;
end

% The algorithm is calibrated for 1000 Hz input files. This correction makes it compatible with other samplerates as well.
fsCorrection = 1000/fs;   
minPeakDist = minPeakDist/fsCorrection;

% ecg = load(inputPath);
% ecg = ecg.data;

%% If the samplerate is high enough, start with a bandpass filter in frequency domain
if fs >= 500
    ecg = FilterWithFFT(ecg);
end


%% ==================== derivative filter ========================== %%
% ------ H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2)) --------- %
if fs ~= 200
 int_c = (5-1)/(fs*1/40);
 b = interp1(1:5,[1 2 0 -2 -1].*(1/8)*fs,1:int_c:5);
else
 b = [1 2 0 -2 -1].*(1/8)*fs;   
end
%filter increases the peak values from 0.5 -> 3500, must downscale them.
ecg_d = filtfilt(b,1,ecg);
ecg_d = ecg_d/max(ecg_d);

% The result signal that the derivative filter produces (ecg_d) makes it
% much more suitable for correct R-peak detection. However the signal
% looses some information which is relevant for displaying it, so the
% plotting still uses the signal without this filter applied.
[peaks,locs] = findpeaks(ecg_d,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDist);
diffLine = diff(locs);

%%  Outlier filtering 
movavgth=1.0;   % Outlier = differs by more than 100% from the moving average
movavg=movmean(diffLine,100,'Endpoints','fill');
% find outliers from moving average and plot them on the figure
wheretoohigh1=find(diffLine>(1+movavgth)*movavg);
wheretoolow1=find(diffLine<(1-movavgth)*movavg);
outliers = numel(wheretoohigh1) + numel(wheretoolow1);
diffMapfilt1=diffLine;
%Mark the outliers as NaN
diffMapfilt1(wheretoolow1)=NaN;diffMapfilt1(wheretoohigh1)=NaN;
%Remove the peak locs and values where the diffmap is NaN
xdata=(1:length(diffLine))';
diffMapfilt1=interp1(xdata(~isnan(diffMapfilt1)),diffMapfilt1(~isnan(diffMapfilt1)),xdata);
diffLine = diffMapfilt1;

%%  So far the outliers are only removed from the diffMap.
% Get rid of them  from the peak locations as well, so they will appear
% correctly on the ECG plot. Note: It only removes the peaks which are too
% close to each other.
peaksToRemove = zeros(numel(wheretoolow1),1);
numOfFalsePeak = 1;
for i = 1: (numel(wheretoolow1)-1)
    if diff([wheretoolow1(i) wheretoolow1(i+1)]) == 1
        peaksToRemove(numOfFalsePeak) = wheretoolow1(i+1);
        numOfFalsePeak = numOfFalsePeak + 1;
    end
end

peaksToRemove(peaksToRemove == 0) = [];
peaks(peaksToRemove) = [];
locs(peaksToRemove) = [];

%% Peak values got deteriorated during the derivative filter. Need to adjust them.
% Peaks should be the max values in their neighborhood
sideThreshold = 3;
adjPeaks = zeros(numel(peaks),1);
adjLocs = zeros(numel(locs),1);
part = zeros(2*sideThreshold+1,1)';
 %to avoid null indexing
locs(locs <= sideThreshold) = sideThreshold +1;
locs(locs > numel(ecg)-sideThreshold) = numel(ecg)-sideThreshold;

for i = 1:numel(locs)
    part = locs(i)-sideThreshold:locs(i)+sideThreshold;
    [M, I]= max(ecg(part));
    adjPeaks(i) = M;
    adjLocs(i) = part(I);
end

%% Plots
figMain = figure('name','ECG Signal With R peaks');
ax = axes;
t = (1:numel(ecg))';
plot(ax, t*fsCorrection, ecg,'Color',[0,0,0]);
grid(ax);
grid(ax, 'minor');
%axis(ax,'square');
ax.GridColor = [ 1 0 0 ];
ax.GridAlpha = 0.5;
ax.MinorGridColor = [1 0 0];
ax.MinorGridAlpha = 0.5;
hold on
plot(adjLocs*fsCorrection,adjPeaks,'bo', 'MarkerFaceColor',[0 0 1]);

xtickformat('%5.5f ms')
xlabel('Time (ms)')
ylabel('Amplitude (mV)')
hold off

h = pan;
h.ButtonDownFilter = @CalcPeaksOfCurrentRegion;
h = zoom;
h.Enable = 'on';
CalcPeaksOfCurrentRegion;
BPM = CalculateBPM(ecg_d,locs);

%Diff map line diagram = Heart rate
figDiffLine = figure('name','R peak difference line diagram');
plot(diffLine*fsCorrection);
grid on;
grid minor;
meanRRinterval = mean(diffLine) * fsCorrection;

formatSpecRr = "Mean RR interval: %0.1f %s ";
title({'RR interval - Line diagram',  sprintf(formatSpecRr,meanRRinterval,'ms.')})
ylabel('RR intervals (ms)')
xlabel('RR intervals over time')

%% Diff map histogram
[figDiffHist, fwhm] = plotHist(diffLine*fsCorrection, meanRRinterval);

%%
function [flag] = CalcPeaksOfCurrentRegion(obj,event_obj)
%Calculates and plots the number of the detected R peaks, and the BPM based
%on the displayed part of the signal. 
        set(0, 'currentfigure', figMain);
        %1st element = leftmost point, 2nd element = rightmost point of the displayed x axis
        actualPlot = xlim; 
        
        lowerBoundary = round(actualPlot(1)/fsCorrection);     
        if lowerBoundary < 1
            lowerBoundary = 1;
        end
        upperBoundary = round(actualPlot(2)/fsCorrection);  
        if upperBoundary >= numel(ecg_d)
            upperBoundary = numel(ecg_d)-1;
        end
        actualPlot = ecg_d(lowerBoundary:upperBoundary);
        
        [~,peakLocsOfRegion] = findpeaks(actualPlot,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDist);
        BPMofRegion = CalculateBPM(actualPlot, peakLocsOfRegion);
        NumberOf_R_Peaks = numel(peakLocsOfRegion);

        %title('ECG Signal with detected R-peaks')
        title({sprintf('ECG Signal with %d detected R peaks.',NumberOf_R_Peaks),...
            sprintf('Heart rate of current region : %0.1f BPM',BPMofRegion)});
        
        flag = false;        
end

function BPM = CalculateBPM(ecg, peakLocations)
%    - normal heart rate in mice awake: ~ 580 bpm
%    - normal heart rate in mice asleep: ~ 520 bpm
%    - normal RR-interval based on awake heart rate: ~ 102 msec.
%    - pathological tachycardia (too fast): ~ 750 bpm
%    - pathological bradycardia (too slow): ~ 440 bpm    
    duration_in_seconds = numel(ecg) / fs;
    duration_in_minutes = duration_in_seconds/60;
    BPM = numel(peakLocations) / duration_in_minutes;
end

end






