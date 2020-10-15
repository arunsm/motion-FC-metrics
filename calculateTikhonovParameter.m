% Script to calculate L2-regularization parameter through minimizing root
% mean square distance between regularized subject precision matrices and
% group average of non-regularized subject precision matrices;
%
% This procedure is described in Pervaiz, Usama, et al.
% "Optimising network modelling methods for fMRI." NeuroImage 211 (2020): 116604.

%% set parameters

addpath(genpath('/cbica/home/mahadeva/matlab/npy-matlab-master'));

tr = 0.72; % relaxation time in seconds
f1 = 0.009; % lower range for frequency domain connectivity metrics (in Hz)
f2 = 0.08; % upper range for frequency domain connectivity metrics (in Hz)

taskTypes = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
atlasTypes = {'gordon', 'yeo_100'};
alphaValues = [0.5, 1, 2, 5, 10, 50, 100, 200];

resultsDir = '/cbica/home/mahadeva/motion-FC-metrics/code/Results/regularizationParameterOptimization/gsr_filter/';

resultsFile = strcat(resultsDir, 'optimal_regularization_parameters.csv');
fid = fopen(resultsFile, 'a');

%% iterate over alpha values and calculate root mean square distance; repeat for all atlases and resting-state scans
for t = 1:numel(taskTypes)
    currentTaskType = taskTypes{t};
    fprintf(currentTaskType); fprintf('\n')
    
    s = importdata(strcat('../data/subjectsList_', currentTaskType, '.csv'));
    subjectsList = s.data;
    nSubjects = numel(subjectsList);
    
    for a = 1:numel(atlasTypes)
        currentAtlasType = atlasTypes{a};
        fprintf(currentAtlasType); fprintf('\n');
        
        switch currentAtlasType
            case 'gordon'
                nParcels = 333;
            case 'yeo_100'
                nParcels = 100;
        end
        
        root_mean_square_distance = zeros(size(alphaValues));
        for i = 1:numel(alphaValues)
            alphai = alphaValues(i);
            fprintf('alpha = %.1f\n', alphai);
            
            precision_matrices_regularized = zeros(nParcels, nParcels, nSubjects);
            precision_matrices_unregularized = zeros(nParcels, nParcels, nSubjects);
            
            for s = 1:nSubjects
                currentSubjectID = subjectsList(s);
                currentFilePath = strcat('/cbica/home/mahadeva/motion-FC-metrics/data/ICAFIX_matrices_nobp/ts/', currentAtlasType, '_', num2str(currentSubjectID), '_', currentTaskType, '_gsr.npy');
                if exist(currentFilePath, 'file')
                    timeSeriesData = readNPY(currentFilePath)';
                    
                    % filter data
                    timeSeriesData = bandpass_filter_butterworth(timeSeriesData, tr, f1, f2);
                    
                    % calculate regularized and unregularized precision
                    % matrices
                    precision_matrices_regularized(:, :, s) = inv(cov(timeSeriesData) + alphai*eye(nParcels));
                    precision_matrices_unregularized(:, :, s) = inv(cov(timeSeriesData));
                    
                else
                    fprintf('Time series data for subject %d not available\n', currentSubjectID);
                end
                
                % calculate root mean square distance between regularized and
                % unregularized matrices
                avge_precision_matrices_unregularized = mean(precision_matrices_unregularized, 3);
                square_distance = sum((precision_matrices_regularized - avge_precision_matrices_unregularized).^2, 3);
                upper_triangular_square_distance = triu(square_distance);
                root_mean_square_distance(i) = sqrt(sum(upper_triangular_square_distance(:)));
            end
        end
        
        % plot alpha values against root mean square distance
        f = figure('Visible', 'Off');
        plot(alphaValues, root_mean_square_distance, 'LineWidth', 2);
        title(strcat(currentAtlasType, ', ', currentTaskType));
        xlabel('regularization parameter (\alpha)');
        ylabel('root mean square distance');
        savePath = strcat(resultsDir, 'regularization_optimization_', currentAtlasType, '_', currentTaskType, '.svg');
        saveas(f, savePath);
        
        % save alpha value associated with minimum root mean square distance
        [~, idx] = min(root_mean_square_distance);
        fprintf(fid, '%s,%s,%d\n', currentAtlasType, currentTaskType, alphaValues(idx));
    end
end

fclose(fid);