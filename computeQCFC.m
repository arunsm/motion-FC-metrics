% Script to compute QC-FC (quality control-functional connectivity) metrics
% fraction significant edges, median absolute correlation

%% set parameters

preprocessingVariants = {'gsr_filter', 'gsr_nofilter', 'nogsr_filter', 'nogsr_nofilter'};
FC_methods = {'Pearson', 'Spearman', 'PartialCorrelation', 'TikhonovPartialCorrelation', 'Coherence', 'WaveletCoherence', 'MutualInformation', 'MutualInformationTime'};
atlasTypes = {'gordon', 'yeo_100'};
restingStateScans = {'REST1_LR', 'REST1_RL','REST2_LR', 'REST2_RL'};

subjectDemographics = readtable('../data/Covariates/S1200_Release_Subjects_Demographics.csv'); % reading demographic information
nSubjects_total = numel(subjectDemographics.Subject);

%% looping through all preprocessing variants, FC methods, atlases and scans
for p = 1:numel(preprocessingVariants)
    currentPreprocessingVariant = preprocessingVariants{p};
    fprintf(currentPreprocessingVariant); fprintf('\n')
    
    saveResultsFolder = strcat('Results/', date, filesep, currentPreprocessingVariant, filesep);
    if ~exist(saveResultsFolder)
        mkdir(saveResultsFolder);
    end
    
    path2FC_matrices = strcat('../data/FunctionalConnectivityMatrices_', currentPreprocessingVariant, filesep);
    currentPipeline = 'FIX_matrices'; % computing QC-FC only for ICA_FIX pipeline
    
    for fc = 1:numel(FC_methods)
        current_FC_method = FC_methods{fc};
        fprintf(current_FC_method); fprintf('\n')
        
        for a = 1:numel(atlasTypes)
            currentAtlasType = atlasTypes{a};
            fprintf(currentAtlasType); fprintf('\n');
            
            TRT_parameters = {}; % cell array with test-retest reliability parameters
            edgeWeights_allScans = cell(nSubjects_total, 5); % cell array of edge weights for all scans
            motion_allScans = cell(nSubjects_total, 5); % cell array of motion for all scans
            edgeWeights_allScans(:, 1) = num2cell(subjectDemographics.Subject);
            motion_allScans(:, 1) = num2cell(subjectDemographics.Subject);
            
            for t = 1:numel(restingStateScans)
                currentRestingStateScan = restingStateScans{t};
                fprintf(currentRestingStateScan); fprintf('\n');
                
                QCFC_parameters = {}; % cell array with QC-FC parameters
                
                %% compile motion and demographic information
                path2MotionFiles = strcat('../data/Motion_S1200/rfMRI_', currentRestingStateScan, '/');
                path2MovementRegressorFiles = strcat('../data/MovementRegressors_S1200/rfMRI_', currentRestingStateScan, '/');
                d = dir(strcat(path2MotionFiles, '*_Movement_RelativeRMS_mean.txt'));
                fnme = {d.name};
                nSubjects_motion = numel(fnme); % number of subjects for whom motion data is available
                subject_ID_motion = zeros(nSubjects_motion, 5); % matrix containing [subject ID, relative RMS motion, mean FD (after filtering), age, gender]
                
                for i = 1:nSubjects_motion
                    currentMotionFileName = fnme{i};
                    k = strfind(currentMotionFileName, '_');
                    currentSubjectID = str2double(currentMotionFileName(1:k(1)-1)); % pulling out subject ID (numeric sequence before first '_')
                    
                    % import mean relative RMS motion
                    currentFilePath = strcat(path2MotionFiles, fnme{i});
                    currentMeanRelativeRMSmotion = importdata(currentFilePath);
                    
                    % calculate mean FD from filtered realignment
                    % parameters
                    currentFilePath = strcat(path2MovementRegressorFiles, strrep(fnme{i}, 'RelativeRMS_mean', 'Regressors'));
                    currentRP = importdata(currentFilePath); % import realignment parameters
                    currentRP = currentRP(:, 1:6); % read in 6 motion parameters
                    currentRP(:, 4:6) = deg2rad(currentRP(:, 4:6)); % convert rotational parameters from degrees to radians
                    filteredRP = filterRealignmentParameters(currentRP); % filter realignment parameters according to Fair et al, NeuroImage (2020)
                    FD = calculatePowerFD(filteredRP);
                    meanFD = mean(FD);
                    
                    currentAgeString = subjectDemographics.Age{subjectDemographics.Subject == currentSubjectID};
                    switch currentAgeString
                        case '22-25'
                            currentAge = (22+25)/2;
                        case '26-30'
                            currentAge = (26+30)/2;
                        case '31-35'
                            currentAge = (31+35)/2;
                        case '36+'
                            currentAge = 36;
                    end
                    
                    currentGenderString = subjectDemographics.Gender{subjectDemographics.Subject == currentSubjectID};
                    switch currentGenderString
                        case 'F'
                            currentGender = 1;
                        case 'M'
                            currentGender = 0;
                    end
                    
                    subject_ID_motion(i, :) = [currentSubjectID currentMeanRelativeRMSmotion meanFD currentAge currentGender];
                end
                
                %% read in adjacency matrices for each subject
                currentReadPath = strcat(path2FC_matrices, currentAtlasType, '*', currentRestingStateScan, '_', currentPipeline, '_', current_FC_method, '.mat');
                d = dir(currentReadPath);
                fnme = {d.name};
                nSubjects_FC = numel(fnme); % number of subjects for whom functional connectivity is available
                
                motion_edgeWeight_allSubjects = {}; % main cell array with fields; {subjectID, meanRelativeMotion, edge weights, age, gender} - one for each FC method
                
                for i = 1:nSubjects_FC
                    currentFile = fnme{i};
                    k = strfind(currentFile(length(currentAtlasType)+1:end), '_'); k = k+length(currentAtlasType); % finding '_' characters, skipping atlas name
                    currentSubjectID = str2double(currentFile(k(1)+1:k(2)-1)); % pulling out subject ID (numeric sequence between first and second '_')
                    if sum(ismember(subject_ID_motion(:, 1), currentSubjectID)) == 0
                        fprintf('motion info for subject %d not available\n', currentSubjectID);
                        continue;
                    else
                        currentMeanRelativeRMSmotion = subject_ID_motion(subject_ID_motion(:, 1) == currentSubjectID, 2); % reading in mean relative RMS motion
                        currentMeanFD = subject_ID_motion(subject_ID_motion(:, 1) == currentSubjectID, 3); % reading in FD after filtering realignment parameters
                        currentAge = subject_ID_motion(subject_ID_motion(:, 1) == currentSubjectID, 4); % reading in age
                        currentGender = subject_ID_motion(subject_ID_motion(:, 1) == currentSubjectID, 5); % reading in gender
                    end
                    
                    currentFilePath = strcat(path2FC_matrices, filesep, currentFile); % reading in adjacency matrix, variable name 'AdjMat'
                    load(currentFilePath, 'AdjMat');
                    
                    %AdjMat = atanh(AdjMat); % applying Fisher z-transform (implement only for correlation based metrics)
                        
                    % converting upper triangular matrix to symmetric matrix
                    if istriu(AdjMat)
                        AdjMat = (AdjMat+AdjMat' - eye(size(AdjMat,1)).*diag(AdjMat));
                    end
                    
                    nNodes = size(AdjMat, 1);
                    edgeWeights = computeEdgeWeights(AdjMat); % computing edge weights
                    nEdges = size(edgeWeights, 2);
                    motion_edgeWeight_allSubjects = [motion_edgeWeight_allSubjects; {currentSubjectID, currentMeanRelativeRMSmotion, currentMeanFD, edgeWeights, currentAge, currentGender}];
                    
                    motion_allScans{subjectDemographics.Subject == currentSubjectID, t+1} = currentMeanRelativeRMSmotion;
                    edgeWeights_allScans{subjectDemographics.Subject == currentSubjectID, t+1} = edgeWeights;
                end
                
                nSubjects_completeData = size(motion_edgeWeight_allSubjects, 1); % number of subjects with motion and FC data
                
                %% compute average edge weights across subjects
                averageEdgeWeights = zeros(size(motion_edgeWeight_allSubjects{1, 4}));
                for i = 1:nSubjects_completeData
                    averageEdgeWeights = averageEdgeWeights + motion_edgeWeight_allSubjects{i, 4};
                end
                averageEdgeWeights = averageEdgeWeights/nSubjects_completeData;
                
                %% compute correlations between motion and edge strengths
                QCFC_correlations = [];
                QCFC_correlations_meanFD = [];
                QCFC_correlations_noCorrections = [];
                QCFC_correlations_significance = [];
                QCFC_correlations_meanFD_significance = [];
                QCFC_correlations_significance_noCorrections = [];
                QCFC_correlations_significance_absoluteValues = [];
                QCFC_correlations_significance_zeroedOut = [];
                
                for e = 1:nEdges % loop over all edges
                    edgeWeights = [];
                    for i = 1:nSubjects_completeData
                        edgeWeights = [edgeWeights; motion_edgeWeight_allSubjects{i, 4}(e)];
                    end
                    
                    meanRelativeRMSmotion = cell2mat(motion_edgeWeight_allSubjects(:, 2));
                    meanFD = cell2mat(motion_edgeWeight_allSubjects(:, 3));
                    age = cell2mat(motion_edgeWeight_allSubjects(:, 5));
                    gender = cell2mat(motion_edgeWeight_allSubjects(:, 6));
                    
                    [rho, pValue] = partialcorr(edgeWeights, meanFD, [age gender], 'Rows', 'complete');
                    QCFC_correlations_meanFD = [QCFC_correlations_meanFD rho];
                    QCFC_correlations_meanFD_significance = [QCFC_correlations_meanFD_significance pValue];
                    
                    [rho, pValue] = corr(edgeWeights, meanRelativeRMSmotion, 'Rows', 'complete');
                    QCFC_correlations_noCorrections = [QCFC_correlations_noCorrections rho];
                    QCFC_correlations_significance_noCorrections = [QCFC_correlations_significance_noCorrections pValue];
                                        
                    [rho, pValue] = partialcorr(edgeWeights, meanRelativeRMSmotion, [age gender], 'Rows', 'complete');
                    QCFC_correlations = [QCFC_correlations rho];
                    QCFC_correlations_significance = [QCFC_correlations_significance pValue];
                    
                    edgeWeights_absoluteValues = abs(edgeWeights);
                    [~, pValue] = partialcorr(edgeWeights_absoluteValues, meanRelativeRMSmotion, [age gender], 'Rows', 'complete');
                    QCFC_correlations_significance_absoluteValues = [QCFC_correlations_significance_absoluteValues pValue];
                    
                    edgeWeights_zeroedOut = edgeWeights;
                    edgeWeights_zeroedOut(edgeWeights_zeroedOut<0) = 0;
                    [~, pValue] = partialcorr(edgeWeights_zeroedOut, meanRelativeRMSmotion, [age gender], 'Rows', 'complete');
                    QCFC_correlations_significance_zeroedOut = [QCFC_correlations_significance_zeroedOut pValue];
                end
                
                % correction for multiple comparisons
                FDR = mafdr(QCFC_correlations_significance, 'BHFDR', true);
                fractionSignificantEdges_FDR = sum(FDR<0.05)/numel(FDR);
                
                fractionSignificantEdges_meanFD = sum(QCFC_correlations_meanFD_significance<0.05)/numel(QCFC_correlations_meanFD_significance);
                fractionSignificantEdges_noFDR = sum(QCFC_correlations_significance<0.05)/numel(QCFC_correlations_significance);
                fractionSignificantEdges_noCorrections_noFDR = sum(QCFC_correlations_significance_noCorrections<0.05)/numel(QCFC_correlations_significance_noCorrections);
                fractionSignificantEdges_absoluteValues_noFDR = sum(QCFC_correlations_significance_absoluteValues<0.05)/numel(QCFC_correlations_significance_absoluteValues);
                fractionSignificantEdges_zeroedOut_noFDR = sum(QCFC_correlations_significance_zeroedOut<0.05)/numel(QCFC_correlations_significance_zeroedOut);
                subjectCovariates = [cell2mat(motion_edgeWeight_allSubjects(:, 1)), age, gender, meanRelativeRMSmotion, meanFD];
                
                QCFC_parameters = {fractionSignificantEdges_FDR, fractionSignificantEdges_noFDR, fractionSignificantEdges_meanFD, ...
                    fractionSignificantEdges_noCorrections_noFDR, fractionSignificantEdges_absoluteValues_noFDR, ...
                    fractionSignificantEdges_zeroedOut_noFDR, QCFC_correlations, QCFC_correlations_noCorrections, ...
                    averageEdgeWeights, nSubjects_completeData, subjectCovariates};
                
                QCFC_parameters = cell2table(QCFC_parameters, 'VariableNames', {'fractionSignificantEdges_FDR', 'fractionSignificantEdges_noFDR', 'fractionSignificantEdges_meanFD', ...
                    'fractionSignificantEdges_noCorrections_noFDR', 'fractionSignificantEdges_absoluteValues_noFDR', 'fractionSignificantEdges_zeroedOut_noFDR', ...
                    'QCFC_correlations', 'QCFC_correlations_noCorrections', 'averageEdgeWeights', 'nSubjects', 'subjectCovariates'});
                save(strcat(saveResultsFolder, 'QCFC_parameters_', current_FC_method, '_', currentAtlasType, '_', currentRestingStateScan, '.mat'), 'QCFC_parameters');
            end
            
            %% computing test-retest reliability
            
            nScans = 4;
            % reducing data to subjects with all 4 scans
            flag_all4Scans = ones(nSubjects_total, 1);
            for i = 1:nSubjects_total
                if isempty(edgeWeights_allScans{i, 2}) || isempty(edgeWeights_allScans{i, 3}) || isempty(edgeWeights_allScans{i, 4}) || isempty(edgeWeights_allScans{i, 5})
                    flag_all4Scans(i) = 0;
                end
            end
            idx_all4Scans = find(flag_all4Scans);
            edgeWeights_allScans = edgeWeights_allScans(idx_all4Scans, :);
            motion_allScans = motion_allScans(idx_all4Scans, :);
            nSubjects_all4scans = size(edgeWeights_allScans, 1);
             
            % compute fingerprinting accuracy for all combinations of scans      
            fingerprintingAccuracy = zeros(1, 12); % 12 combinations of target and database from 4 resting state scans
            ctr = 1;
            for t = 1:nScans
                targetEdgeWeights = edgeWeights_allScans(:, t+1);
                for d = 1:nScans
                    if t ~= d
                        databaseEdgeWeights = edgeWeights_allScans(:, d+1);
                        idx_match = zeros(1, nSubjects_all4scans);
                        for i = 1:nSubjects_all4scans
                            currentTargetEdgeWeights = targetEdgeWeights{i};
                            idx_match(i) = calculateFingerprintAccuracy(databaseEdgeWeights, currentTargetEdgeWeights, i);
                        end
                        fingerprintingAccuracy(ctr) = sum(idx_match)/nSubjects_all4scans;
                        ctr = ctr + 1;
                    end
                end
            end
            
            % compute intra-class correlation for edge weights
            ICC_allEdges = zeros(1, nEdges);
            for e = 1:nEdges
                % build a nSubjects x nScans matrix for each edge
                TRT_matrix_edges = zeros(nSubjects_all4scans, nScans);
                for i = 1:nSubjects_all4scans
                    for j = 1:nScans
                        TRT_matrix_edges(i, j) = edgeWeights_allScans{i, j+1}(e);
                    end
                end
                
                ICC = computeICC(TRT_matrix_edges);
                ICC_allEdges(e) = ICC;                
            end
            
            % compute intra-class correlation for subject motion
            TRT_matrix_motion = zeros(nSubjects_all4scans, nScans);
            for i = 1:nSubjects_all4scans
                for j = 1:nScans
                    TRT_matrix_motion(i, j) = motion_allScans{i, j+1};
                end
            end
            
            ICC_motion = computeICC(TRT_matrix_motion);
            
            TRT_parameters = {fingerprintingAccuracy, ICC_allEdges, ICC_motion, nSubjects_all4scans};
            TRT_parameters = cell2table(TRT_parameters, 'VariableNames', {'fingerprintingAccuracy', 'ICC_allEdges', 'ICC_motion', 'nSubjects_all4scans'});
            save(strcat(saveResultsFolder, 'TRT_parameters_', current_FC_method, '_', currentAtlasType, '.mat'), 'TRT_parameters');
            
        end
    end
end