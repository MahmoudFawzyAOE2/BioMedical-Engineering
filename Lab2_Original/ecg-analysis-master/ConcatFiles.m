function ConcatFiles(folderName)
% Concats the .mat files found in the folder specified by 'folderName'
% It should contain a full path to a folder with .mat files in it.

if 7 ~= exist(folderName,'dir')
    error('Cannot locate the folder: %s',folderName)
end

fds = fileDatastore(folderName, 'ReadFcn', @importdata);
fullFileNames = fds.Files;
numFiles = length(fullFileNames);
value = [];
for i = 1 : numFiles
    ecg = load(fullFileNames{i});
    ecg = ecg.data;
    value = cat(1, value, ecg);
end
save(strcat(folderName,'\FullSegment.mat'), 'value');

