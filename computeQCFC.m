% Script to compute QC-FC (quality control-functional connectivity) metrics
% fraction significant edges, median absolute correlation

saveResultsFolder = strcat('Results/', date, '/');
if ~exist(saveResultsFolder)
    mkdir(saveResultsFolder);
end

FC_methods = {'Pearson', 'Spearman', 'Coherence', 'WaveletCoherence', 'MutualInformation', 'MutualInformationTime'};
atlasTypes = {'gordon', 'yeo_100'};
taskTypes = {'REST1_LR', 'REST1_RL','REST2_LR', 'REST2_RL'};
motionCorrectionMethods = {'CompCor_matrices', 'FIX_matrices'};
nPipelines = numel(motionCorrectionMethods);

for fc = 1:numel(FC_methods)
    current_FC_method = FC_methods{fc};
    fprintf(current_FC_method); fprintf('\n')
    
    for a = 1:numel(atlasTypes)
        currentAtlasType = atlasTypes{a};
        fprintf(currentAtlasType); fprintf('\n');

        for t = 1:numel(taskTypes)
            currentTaskType = taskTypes{t};
            fprintf(currentTaskType); fprintf('\n'); 
	        fractionSignificantEdges = cell(3, nPipelines);
            medianAbsoluteCorrelation = cell(3, nPipelines);

            for p = 1:nPipelines
                currentPipeline = motionCorrectionMethods{p};
                
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
                
                % read in adjacency matrices for each subject
                folderPath = '../Data/FunctionalConnectivityMatrices/';
                d = dir(strcat(folderPath, currentAtlasType, '_', '*', currentTaskType, '_', currentPipeline, '_', current_FC_method, '.mat'));
                fnme = {d.name};
                nSubjects_ts = numel(fnme); % number of subjects for whom functional connectivity is available
                
                motion_ew_allSubjects = {}; % main cell array with fields; {subjectID, meanRelativeMotion, edge weights} - one for each FC method
                
                for i = 1:nSubjects_ts
                    currentFile = fnme{i};
                    k = strfind(currentFile(length(currentAtlasType)+1:end), '_'); k = k+length(currentAtlasType); % finding '_' characters, skipping atlas name
                    currentSubjectID = str2double(currentFile(k(1)+1:k(2)-1)); % pulling out subject ID (numeric sequence between first and second '_')
                    if sum(ismember(subject_ID_motion(:, 1), currentSubjectID)) == 0
                        fprintf('motion info for subject %d not available\n', currentSubjectID);
                        continue;
                    else
                        currentMeanRelativeMotion = subject_ID_motion(subject_ID_motion(:, 1) == currentSubjectID, 2); % reading in motion data - relative mean RMS motion
                    end
                    
                    currentFilePath = strcat(folderPath, filesep, currentFile); % reading in adjacency matrix, variable name 'AdjMat'
                    load(currentFilePath);
                    
                    %AdjacencyMatrix(AdjacencyMatrix<0) = 0; % setting negative correlations to zero
                    AdjMat = atanh(AdjMat); % Fisher z-transform
                    % converting upper triangular matrix to symmetric matrix
                    if istriu(AdjMat)
                        AdjMat = (AdjMat+AdjMat' - eye(size(AdjMat,1)).*diag(AdjMat));
                    end
                    
                    nNodes = size(AdjMat, 1);
                    edgeWeights = computeEdgeWeights(AdjMat); % computing edge weights
                    nEdges = size(edgeWeights, 2);
                    motion_ew_allSubjects = [motion_ew_allSubjects; {currentSubjectID, currentMeanRelativeMotion, edgeWeights}];
                end
                
                % compute correlations between mean RMS motion and edge strengths
                edgeMotionCorr = [];
                edgeMotionCorr_significance = [];
                
                nSubjects_completeData = size(motion_ew_allSubjects, 1);
                for e = 1:nEdges % loop over all edges
                    ew = [];
                    for i = 1:nSubjects_completeData
                        ew = [ew; motion_ew_allSubjects{i, 3}(e)];
                    end
                    
                    motion = cell2mat(motion_ew_allSubjects(:, 2));
                    
                    [rho, pValue] = corr(ew, motion, 'Rows', 'complete');
                    edgeMotionCorr = [edgeMotionCorr rho];
                    edgeMotionCorr_significance = [edgeMotionCorr_significance pValue];
                end
                
                % correct for multiple comparisons
                FDR = mafdr(edgeMotionCorr_significance, 'BHFDR', true);
                fractionSignificantEdges{1, p} = currentPipeline;
                fractionSignificantEdges{2, p} = sum(FDR<0.05)/nEdges;
                fractionSignificantEdges{3, p} = nSubjects_completeData;
                medianAbsoluteCorrelation{1, p} = currentPipeline;
                medianAbsoluteCorrelation{2, p} = median(abs(edgeMotionCorr(~isnan(edgeMotionCorr))));
                medianAbsoluteCorrelation{3, p} = nSubjects_completeData;
                medianAbsoluteCorrelation{4, p} = edgeMotionCorr;
            end
            
            save(strcat(saveResultsFolder, 'fractionSignificantEdges_', current_FC_method, '_', currentAtlasType, '_', currentTaskType, '_allPipelines.mat'), 'fractionSignificantEdges');
            save(strcat(saveResultsFolder, 'medianAbsoluteCorrelation_', current_FC_method, '_', currentAtlasType, '_', currentTaskType, '_allPipelines.mat'), 'medianAbsoluteCorrelation');
        end
    end
end
