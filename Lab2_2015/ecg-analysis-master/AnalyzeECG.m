function AnalyzeECG(folderName, fs, movavgth, minPeakHeight, minPeakDist)
% Applies R-peak detection on multiple files, prints the plots to PDFs,
% and information about the analysis to a text file.

%   Collects the .mat files from the specified folder, applies
%   peakdetection, then prints ectopic beats to separate pdf
%   files, and the difference histogram and line diagram.
%   Prints information about the analysis to a text file:
%   -BPM
%   -Number of R peaks
%   -Mean RR interval
%   -Full-Width at Half-Maximum
%   -Number of outliers
%   -The input arguments of the analysis
%
%   INPUTS:
%       folderName - A full path to the folder containing the .mat files.
%       fs - Frequency, given in Hertz. Mandatory input argument.
%       movavgth - Ectopic beats will be classified as peaks differing more
%   than movavgth% from the moving average of the RR interval. Set it to
%   0.3, for a 30% threshold.
%
%   Example: AnalyzeECG('C:\FolderWithFiles', 1000, 0.3, 0.5, 50);
%            AnalyzeECG('C:\FolderWithFiles', 1000);



%% Argument check
if nargin < 2 || isempty(fs)
    error('Input argument fs is required.')
end
if nargin < 3 || isempty(movavgth)
    movavgth = 0.3;
end
if nargin < 4 || isempty(minPeakHeight)
    minPeakHeight = 0.5;
end
if nargin < 5 || isempty(minPeakDist)
    minPeakDist = 50;
end

if 7 ~= exist(folderName,'dir')
    error('Cannot locate the folder: %s',folderName)
end

%% Processing files
files = dir(strcat(folderName,'\*.mat'));
numOfFiles = numel(files);

for i = 1: numOfFiles
    fullFileName = files(i).name;
    split = strsplit(fullFileName, '.');
    fileName = split{1};
    status = mkdir(strcat(folderName,'\',fileName));
    if status == 0
        error('Could not create directory for file: %s',fileName)
    end
    
    %calling peakDetect
    [figMain, figDiffLine, figDiffHist, diffLine, locs, BPM, fwhm, outliers] = PeakDetect(strcat(folderName,'\',fullFileName), fs, minPeakHeight, minPeakDist);
    
    %% Ectopic beat detection based on aberrant RR-intervals.
    % Scanning through the moving average RR-interval of the entire ECG trace.
    
    movavgDiff=movmean(diffLine,100,'Endpoints','fill');
    % ectopicBeats = peaks differing more than movavgth from the moving average of the RR interval
    ectopicBeats = find(diffLine>(1+movavgth)*movavgDiff | diffLine<(1-movavgth)*movavgDiff);
    if numel(ectopicBeats) > 0
        distOutliers = diff(ectopicBeats);
        clusters = zeros(size(ectopicBeats));
        n = 1;
        clusters(n) = ectopicBeats(1);
        n = n+1; 
        % Cluster: the distance between 2 outlier is below a threshold (10)
        % In that case, only keep 1 element of that cluster.
        for i = 1:numel(distOutliers)
            if (distOutliers(i) > 10)
               clusters(n) = ectopicBeats(i+1);
               n = n+1; 
            end    
        end
        clusters(clusters == 0) = [];
                
        %% Print the Region of interests (clusters) into separate PDF files.
        set(0, 'currentfigure', figMain)
        figure(figMain)
        for ii=1:numel(clusters)
            title(['ROI #',num2str(ii),' of file: ',fileName]);
            ROI = locs(clusters(ii));
            xlim([ROI-300 ROI+300]);
            print(figMain,'-dpdf', [folderName '\' fileName '\ROI' '_' num2str(ii)],'-bestfit');
        end
    end

    
    %% Print the Difference line and histogram.
    figure(figDiffLine)
    print(figDiffLine,'-dpdf', [folderName '\' fileName '\DifferenceLine'],'-bestfit');
    close(figDiffLine)
    figure(figDiffHist)
    print(figDiffHist,'-dpdf', [folderName '\' fileName '\DifferenceHistogram'],'-bestfit');
    close(figDiffHist)
    
    %% Write BPM, Number of R peaks, Mean RR interval, Number of outliers into a text file.
    NumberOf_R_Peaks = numel(locs);
    meanRRinterval = mean(diffLine);
    fileID = fopen([folderName '\' fileName '\' 'DataSummary.txt'],'w');
    
    % Printing the result of the analysis
    fprintf(fileID,' Summary of ECG analysis \n\n');
    formatSpec = 'BPM = %3.2f \n';
    fprintf(fileID,formatSpec,BPM);
    formatSpec = 'Number of R peaks = %d \n';
    fprintf(fileID,formatSpec,NumberOf_R_Peaks);
    formatSpec = 'Mean RR interval = %f \n';
    fprintf(fileID,formatSpec,meanRRinterval);
    formatSpec = 'Full-Width at Half-Maximum = %d \n';
    fprintf(fileID,formatSpec,fwhm);
    formatSpec = 'Number of outliers = %d \n\n';
    fprintf(fileID,formatSpec,outliers);
    
    % Printing the input arguments
    fprintf(fileID,' Input arguments of the analysis \n\n');
    formatSpec = 'Frequency = %d Hz \n';
    fprintf(fileID,formatSpec,fs);
    formatSpec = 'Moving average thershold = %3.2f \n';
    fprintf(fileID,formatSpec,movavgth);
    formatSpec = 'minPeakHeight = %3.3f \n';
    fprintf(fileID,formatSpec,minPeakHeight);
    formatSpec = 'minPeakDist = %d \n';
    fprintf(fileID,formatSpec,minPeakDist);
    
    fclose(fileID);
    close all
end
end