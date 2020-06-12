% Script to compute functional connectivity from time series data for
% one subject

function engine_FunctionalConnectivity_gordon(subjectID, taskIndicator)

% If obtaining taskIndicator from array job

% switch taskIndicator
%     case 1
% 	taskType = 'REST1_LR';
%     case 2
% 	taskType = 'REST1_RL';
%     case 3
% 	taskType = 'REST2_LR';
%     case 4
% 	taskType = 'REST2_RL';
% end
%
% fprintf('Processing data for %s\n', taskType)

taskTypes = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
%taskTypes = {'REST2_LR'};

%% add necessary toolboxes to path
addpath(genpath('/cbica/home/mahadeva/matlab/wavelet-coherence-master'));
addpath(genpath('/cbica/home/mahadeva/matlab/Functional Connectivity Toolbox_updated'));
addpath(genpath('/cbica/home/mahadeva/matlab/npy-matlab-master'));

%% set parameters

tr = 0.72; % relaxation time in seconds
f1 = 0.009; % lower range for frequency domain connectivity metrics (in Hz)
f2 = 0.08; % upper range for frequency domain connectivity metrics (in Hz)

FC_methods = {'Pearson', 'PartialCorrelation', 'Spearman', 'Coherence', 'WaveletCoherence', 'MutualInformation', 'MutualInformationTime'};
%FC_methods = {'MutualInformationTime'};
motionCorrectionMethods = {'ts'};
nPipelines = numel(motionCorrectionMethods);
atlasType = 'gordon';
resultsFolder = '/cbica/home/mahadeva/motion-FC-metrics/data/FunctionalConnectivityMatrices_nogsr_filter';
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder)
end

%% run FC analysis for all pipelines and all methods
for t = 1:numel(taskTypes)
    taskType = taskTypes{t};
    fprintf(taskType); fprintf('\n')
    
    s = importdata(strcat('../data/subjectsList_', taskType, '.csv'));
    subjectsList = s.data;
    if ~ismember(subjectID, subjectsList)
        fprintf('skipping subject %d, not in subjects list\n', subjectID);
        continue;
    end
    
    for fc = 1:numel(FC_methods)
        currentFCmethod = FC_methods{fc};
        fprintf(currentFCmethod); fprintf('\n')
        for p= 1:nPipelines
            currentPipeline = motionCorrectionMethods{p};
            currentFilePath = strcat('/cbica/home/mahadeva/motion-FC-metrics/data/ICAFIX_matrices_nobp/', currentPipeline, filesep, atlasType, '_', num2str(subjectID), '_', taskType, '_nogsr.npy');
            if exist(currentFilePath, 'file')
                savePath = strcat(resultsFolder, filesep, atlasType, '_', num2str(subjectID), '_', taskType, '_', currentPipeline, '_', currentFCmethod, '.mat');
                if exist(savePath, 'file')
                    fprintf('Adjacency matrix for subject %d, parcellation %s, preprocessing pipeline %s, FC method %s already exists\n', subjectID, atlasType, currentPipeline, currentFCmethod)
                    continue;
                else
                    timeSeriesData = readNPY(currentFilePath)';
                    
                    % filtering data
                    timeSeriesData = bandpass_filter_butterworth(timeSeriesData, tr, f1, f2); 
                    
                    fprintf('Computing functional connectivity for subject %d, parcellation %s, preprocessing pipeline %s, using %s\n', subjectID, atlasType, currentPipeline, currentFCmethod)
                    AdjMat = computeFunctionalConnectivity(timeSeriesData, currentFCmethod, f1, f2, tr);
                    save(savePath, 'AdjMat');
                end
            else
                fprintf('Time series data for subject %d not available\n', subjectID);
            end
        end
    end
end

end
