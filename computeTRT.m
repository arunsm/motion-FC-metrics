% Script to compute test-retest reliability

clear all
nScans = 4;
saveResultsFolder = strcat('Results/', date, '/');

if ~exist(saveResultsFolder)
    mkdir(saveResultsFolder);
end

%% combine all subject IDs
taskTypes = {'REST1_LR', 'REST1_RL','REST2_LR', 'REST2_RL'};
subject_ID_list = []; % array containing all subject IDs

for t = 1:numel(taskTypes)
    currentTaskType = taskTypes{t};
    folderPath = strcat('../Data/Motion_S1200/rfMRI_', currentTaskType, '/');
    d = dir(strcat(folderPath, '*RelativeRMS_mean.txt'));
    fnme = {d.name};
    nSubjects_motion = numel(fnme); % number of subjects for whom motion data is available
    for i = 1:nSubjects_motion
        currentMotionFileName = fnme{i};
        k = strfind(currentMotionFileName, '_');
        currentSubjectID = str2double(currentMotionFileName(1:k(1)-1)); % pulling out subject ID (numeric sequence before first '_')
        subject_ID_list = [subject_ID_list; currentSubjectID];
    end
end

subject_ID_list = unique(subject_ID_list);
nSubjects = numel(subject_ID_list);

%% compute ICC for edge weights and motion across all 4 scans

folderPath_FCmatrices = '../Data/FunctionalConnectivityMatrices/';
folderPath_motion = strcat('../Data/Motion_S1200/');
FC_methods = {'Pearson', 'Spearman', 'Coherence', 'WaveletCoherence', 'MutualInformation', 'MutualInformationTime'};
atlasTypes = {'gordon', 'yeo_100'};
motionCorrectionMethods = {'CompCor_matrices', 'FIX_matrices'};
nPipelines = numel(motionCorrectionMethods);

for fc = 1:numel(FC_methods)
    current_FC_method = FC_methods{fc};
    fprintf('FC method %s\n', current_FC_method)
    
    for a = 1:numel(atlasTypes)
        currentAtlasType = atlasTypes{a};
        fprintf('\tatlas type %s\n', currentAtlasType)
        
        for p = 1:nPipelines
            currentPipeline = motionCorrectionMethods{p};
            fprintf('\t\tpipeline %s\n', currentPipeline)
            
            ew_allScans_allSubjects = {};
            motion_allSubjects = {};
            subject_ctr = 1;
            for i = 1:nSubjects
                currentSubjectID = num2str(subject_ID_list(i));
                %fprintf('\t\t\tsubject %s\n', currentSubjectID)
                
                % read in adjacency matrices for each subject
                d = dir(strcat(folderPath_FCmatrices, currentAtlasType, '_', currentSubjectID, '_*_', currentPipeline, '_', current_FC_method, '.mat'));
                fnme_FC = {d.name};
                nFiles_FC = numel(fnme_FC); % number of files found
                %fprintf('\tfound %d files for %s, %s, %s\n', nFiles, currentAtlasType, currentPipeline, current_FC_method)
                
                % read in motion values
                d = dir(strcat(folderPath_motion, 'rFMRI_*', filesep, currentSubjectID, '_Movement_RelativeRMS_mean.txt'));
                fnme_motion = {d.name};
                nFiles_motion = numel(fnme_motion); % number of files found
                
                if nFiles_FC < nScans
                    %fprintf('\t\t\tLess than 4 scans found; skipping this subject\n')
                    continue;
                    
                elseif nFiles_motion < nScans
                    %fprintf('\t\t\tMotion data for all scans not found; skipping this subject\n')
                    continue;
                        
                else
                    motion_allSubjects{subject_ctr, 1} = currentSubjectID;
                    ew_allScans_allSubjects{subject_ctr, 1} = currentSubjectID;
                    
                    for t = 1:numel(taskTypes)
                        currentTaskType = taskTypes{t};
                        currentFilePath = strcat(folderPath_motion, 'rFMRI_', currentTaskType, filesep, currentSubjectID, '_Movement_RelativeRMS_mean.txt');
                        currentMeanRelativeMotion = importdata(currentFilePath);
                        motion_allSubjects{subject_ctr, t+1} = currentMeanRelativeMotion;
                    end
                    
                    for j = 1:nFiles_FC
                        currentFile = strcat(folderPath_FCmatrices, fnme_FC{j});
                        load(currentFile); % load adjacency matrix
                        
                        AdjMat = atanh(AdjMat); % Fisher z-transform
                        
                        % converting upper triangular matrix to symmetric matrix
                        if istriu(AdjMat)
                            AdjMat = (AdjMat+AdjMat' - eye(size(AdjMat,1)).*diag(AdjMat));
                        end
                        
                        nNodes = size(AdjMat, 1);
                        edgeWeights = computeEdgeWeights(AdjMat); % computing edge weights
                        nEdges = size(edgeWeights, 2);
                        ew_allScans_allSubjects{subject_ctr, j+1} = edgeWeights;
                    end
                    
                    subject_ctr = subject_ctr + 1;
                end
            end
            
            fprintf('\t\t\t\tComputing ICC for edge weights across all scans... \n')
            nSubjects_completeData = size(ew_allScans_allSubjects, 1);
            ICC_allEdges = zeros(1, nEdges);
            for e = 1:nEdges
                % build a nSubjects x nScans matrix for each edge
                TRT_matrix = zeros(nSubjects_completeData, nScans);
                for i = 1:nSubjects_completeData
                    for j = 1:nScans
                        TRT_matrix(i, j) = ew_allScans_allSubjects{i, j+1}(e);
                    end
                end
                
                % compute intra-class correlation
                ICC = computeICC(TRT_matrix);
                ICC_allEdges(e) = ICC;
            end
            
            ICC_stats_FCedgeWeights{1, 1} = 'intra-class correlations for all edges';
            ICC_stats_FCedgeWeights{1, 2} = 'number of subjects (with all 4 scans)';
            ICC_stats_FCedgeWeights{2, 1} = ICC_allEdges;
            ICC_stats_FCedgeWeights{2, 2} = nSubjects_completeData;
            
            % save ICC values for each pipeline,
            save(strcat(saveResultsFolder, 'IntraClassCorrelation_edgeWeights', current_FC_method, '_', currentAtlasType, '_', currentPipeline, '.mat'), 'ICC_stats_FCedgeWeights');
            %save(strcat(saveResultsFolder, 'EdgeWeights_allScans', current_FC_method, '_', currentAtlasType, '_', currentPipeline, '.mat'), 'ew_allScans_allSubjects');
            
            fprintf('\t\t\t\tComputing ICC for motion across all scans... \n')
            % build a nSubjects x nScans matrix for motion
            TRT_matrix = zeros(nSubjects_completeData, nScans);
            for i = 1:nSubjects_completeData
                for j = 1:nScans
                    TRT_matrix(i, j) = motion_allSubjects{i, j+1};
                end
            end
            
            % compute intra-class correlation
            ICC_motion = computeICC(TRT_matrix);
            
            ICC_stats_meanRelativeRMSmotion{1, 1} = 'intra-class correlations for mean relative RMS motion';
            ICC_stats_meanRelativeRMSmotion{1, 2} = 'number of subjects (with all 4 scans)';
            ICC_stats_meanRelativeRMSmotion{2, 1} = ICC_motion;
            ICC_stats_meanRelativeRMSmotion{2, 2} = nSubjects_completeData;
            
            % save ICC values for each pipeline,
            save(strcat(saveResultsFolder, 'IntraClassCorrelation_meanRelativeRMSmotion', current_FC_method, '_', currentAtlasType, '_', currentPipeline, '.mat'), 'ICC_stats_meanRelativeRMSmotion');
            save(strcat(saveResultsFolder, 'meanRelativeRMSmotion_allScans_allSubjects', current_FC_method, '_', currentAtlasType, '_', currentPipeline, '.mat'), 'motion_allSubjects');
        end
    end
end

function IntraClassCorrelation = computeICC(data)

[n,k] = size(data); % n = number of subjects; k = number of scans

% mean per target
mpt = mean(data,2);
% get total mean
tm = mean(data(:));
% within target sum sqrs
tmp = (data - repmat(mpt,1,k)).^2;
WSS = sum(tmp(:));
% within target mean sqrs
WMS = WSS / (n*(k - 1));
% between target sum sqrs
BSS = sum((mpt - tm).^2) * k;
% between targets mean squares
BMS = BSS / (n - 1);

IntraClassCorrelation = (BMS - WMS) / (BMS + (k - 1) * WMS);
end