clear all
clc
 
loaddirectory = 'C:\Users\MateToth\Documents\ECG\Input files\200Hz\Unprocessed';
savedirectory = 'C:\Users\MateToth\Documents\ECG\Input files\200Hz\Unprocessed';
cd(loaddirectory);
files = dir(loaddirectory);
 
numpnts = 20 * 60 * 2000;
numchan = 30;
 
numberfiles = length(files);
filenumber = 0;
samplerate = 2000;
 downsamplefactor = 1; % downsamples, takes every nth point
 nyquist = samplerate/2;
 filterorder = 2;
 highpassfreq = 5; % cutoff frequency in Hz
 highpass = highpassfreq/nyquist; % frequencies higher than 5 Hz are let through.
 [b,a] = butter(filterorder, highpass, 'high');
 
for i = 3:numberfiles;
        %filetype = splitstring(files(i,1).name, '.');
        filetype = strsplit(files(i,1).name, '.');

        if strcmp(filetype{2},'SEG') == true;
            filenumber = filenumber + 1;
            [fid,errormessage] = fopen(files(i,1).name,'r');
            data = fread(fid,[numchan numpnts],'int16=>int16')';
            fclose(fid);
            data = double(data(:,29))*(20/2^16); % choose which trace

            filtered_data = zeros(size(data,1),1);
            for j = 1:size(data,2); % loop for filtering
               filtered_data(:,j) = filtfilt(b,a,data(:,j));
            end

            filtered_downed = zeros(size(data,1)/downsamplefactor,size(data,2),1);
            datadowned = zeros(size(data,1)/downsamplefactor,size(data,2),1);
            for j = 1:size(data,2); % loop for downsampling
                filtered_downed(:,j) = downsample(filtered_data(:,j),downsamplefactor);
                datadowned(:,j) = downsample(data(:,j),downsamplefactor);
            end

            
            %tiledlayout(3,1)
            %ax1 = nexttile;
            %plot(data)
            %ax2 = nexttile;
            %plot(datadowned)
            %ax3 = nexttile;
            %plot(filtered_downed)
            %linkaxes([ax1 ax2 ax3],'xy')
            
            
            data = datadowned;
            clear datadowned

            if filenumber < 10;
               writematrix(data,[savedirectory '00' num2str(filenumber) '.txt']);
     %          xlswrite([savedirectory 'filt00' num2str(filenumber) '.xls'], filtered_downed);
             save([savedirectory '00' num2str(filenumber) '.mat'], 'data');
     %           save([savedirectory 'filt00' num2str(filenumber) '.mat'], 'filtered_downed');
            elseif filenumber > 10 && filenumber < 100;        
               writematrix(data,[savedirectory '0' num2str(filenumber) '.txt']);
     %          xlswrite([savedirectory 'filt0' num2str(filenumber) '.xls'], filtered_downed);
              save([savedirectory '0' num2str(filenumber) '.mat'], 'data');
     %          save([savedirectory 'filt0' num2str(filenumber) '.mat'], 'filtered_downed');
            elseif filenumber > 100;        
               writematrix(data,[savedirectory '' num2str(filenumber) '.txt']);
     %          xlswrite([savedirectory 'filt' num2str(filenumber) '.xls'], filtered_downed);
             save([savedirectory '0' num2str(filenumber) '.mat'], 'data');
     %           save([savedirectory 'filt0' num2str(filenumber) '.mat'], 'filtered_downed');
            end
            disp([files(i,1).name 'loaded']);
        end
end

%dt=1/(samplerate/downsamplefactor);
%time=dt:dt:(20*60);
%plot(time,data, 'DisplayName','Non-filtered')
%hold on
%plot(time,filtered_downed, 'DisplayName',['HighPass = ' num2str(highpassfreq) ' Hz'], 'Color','r')