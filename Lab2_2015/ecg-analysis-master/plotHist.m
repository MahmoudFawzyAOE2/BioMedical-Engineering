function [figDiffHist, fwhm]= plotHist(diffLine, meanRRinterval)
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
figDiffHist= figure('name','R peak difference histogram');
h = histogram(diffLine);
hold on;
grid on;
grid minor;

data = h.Values;
% Find the half max value.
halfMax = (min(data) + max(data)) / 2;
% Find where the data first drops below half the max.
index1 = find(data >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(data >= halfMax, 1, 'last');
fwhm = index2-index1 + 1; % FWHM in number of bins -> in ms.

% Plot a line at the half-max of the histogram
xCoords = h.BinLimits;
yCoords = [halfMax halfMax];
line(xCoords, yCoords, 'LineWidth', 2, 'Color', 'r');


formatSpecRr = "Mean RR interval: %0.1f %s";
formatSpecFwhm = "Full-Width at Half-Maximum: %0.0f %s";

title({'RR Interval - Histogram',  sprintf(formatSpecRr,meanRRinterval,'ms.'),...
        sprintf(formatSpecFwhm,fwhm,'ms.')})
%title({'RR interval - Histogram',  sprintf('Mean RR interval = %0.1f',meanRRinterval),...
%        sprintf('Full-Width at Half-Maximum = %0.0f',fwhm)})

ylabel('Number of occasions')
xlabel('RR intervals (ms)')
hold off;



