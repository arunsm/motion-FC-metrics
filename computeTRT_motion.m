% Script to compute test-retest reliability of motion

clear all
nScans = 4;
saveResultsFolder = strcat('Results/', date, '/');

if ~exist(saveResultsFolder)
    mkdir(saveResultsFolder);
end

%% read in motion data
taskTypes = {'REST1_LR', 'REST1_RL','REST2_LR', 'REST2_RL'};
TRT_matrix = zeros(nSubjects_completeData, nScans);

% read in motion time series
folderPath = strcat('../Data/Motion_S1200/rfMRI_', currentTaskType, '/');
d = dir(strcat(folderPath, '*RelativeRMS_mean.txt'));
fnme = {d.name};
nSubjects_motion = numel(fnme); % number of subjects for whom motion data is available
subject_ID_motion = zeros(nSubjects_motion, 2); % matrix containing subject ID and relative RMS motion

for i = 1:nSubjects_motion
    currentMotionFileName = fnme{i};
    k = strfind(currentMotionFileName, '_');
    currentSubjectID = str2double(currentMotionFileName(1:k(1)-1)); % pulling out subject ID (numeric sequence before first '_')
    currentFilePath = strcat(folderPath, fnme{i});
    currentMeanRelativeMotion = importdata(currentFilePath);
    subject_ID_motion(i, :) = [currentSubjectID currentMeanRelativeMotion];
end