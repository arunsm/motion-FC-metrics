% script to generate figures 1-7 in the main text and figures 1-7 in the supplementary information 

%% set parameters and add dependencies to path

clear all;
close all;
taskTypes = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
FC_methods = {'Pearson', 'Spearman', 'MutualInformationTime', 'Coherence', 'WaveletCoherence', 'MutualInformation'}; nFC_methods = numel(FC_methods);
parcellations = {'gordon', 'yeo_100'};
pipelines = {'CompCor_matrices', 'FIX_matrices'}; nPipelines = numel(pipelines);
%pipelines = {'FIX_matrices'}; nPipelines = numel(pipelines);

labels = FC_methods;

addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/violin'));
addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/BCT'));
addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/beeswarm-master'));

%% table 1 - summary statistics on motion in cohort

clc
resultsFolder = 'Results/08-Apr-2020/';

load(strcat(resultsFolder, 'QCFC_parameters_Spearman_gordon_REST1_LR_allPipelines.mat'));
covariates = cell2table(QCFC_parameters{14, 3}, 'VariableNames', {'subjectID', 'gender', 'age', 'relativeRMSmotion'});
fprintf('REST1_LR (n=%i)\n', QCFC_parameters{13, 3})
fprintf('mean relative RMS motion = %0.3f\n', mean(covariates.relativeRMSmotion))
fprintf('std relative RMS motion = %0.3f\n', std(covariates.relativeRMSmotion))

load(strcat(resultsFolder, 'QCFC_parameters_Spearman_gordon_REST1_RL_allPipelines.mat'));
covariates = cell2table(QCFC_parameters{14, 3}, 'VariableNames', {'subjectID', 'gender', 'age', 'relativeRMSmotion'});
fprintf('REST1_RL (n=%i)\n', QCFC_parameters{13, 3})
fprintf('mean relative RMS motion = %0.3f\n', mean(covariates.relativeRMSmotion))
fprintf('std relative RMS motion = %0.3f\n', std(covariates.relativeRMSmotion))

load(strcat(resultsFolder, 'QCFC_parameters_Spearman_gordon_REST2_LR_allPipelines.mat'));
covariates = cell2table(QCFC_parameters{14, 3}, 'VariableNames', {'subjectID', 'gender', 'age', 'relativeRMSmotion'});
fprintf('REST2_LR (n=%i)\n', QCFC_parameters{13, 3})
fprintf('mean relative RMS motion = %0.3f\n', mean(covariates.relativeRMSmotion))
fprintf('std relative RMS motion = %0.3f\n', std(covariates.relativeRMSmotion))

load(strcat(resultsFolder, 'QCFC_parameters_Spearman_gordon_REST2_RL_allPipelines.mat'));
covariates = cell2table(QCFC_parameters{14, 3}, 'VariableNames', {'subjectID', 'gender', 'age', 'relativeRMSmotion'});
fprintf('REST2_RL (n=%i)\n', QCFC_parameters{13, 3})
fprintf('mean relative RMS motion = %0.3f\n', mean(covariates.relativeRMSmotion))
fprintf('std relative RMS motion = %0.3f\n', std(covariates.relativeRMSmotion))

%% figure 1 - pairwise correlation plots between edge weights obtained using all 6 metrics

resultsFolder = 'Results/02-Aug-2019/';

for p = 1:numel(parcellations)
    currentParcellation = parcellations{p};
    
    for i = 1:nPipelines
        currentPipeline = pipelines{i};
        
        for t = 1:numel(taskTypes)
            currentTaskType = taskTypes{t};
            edgeWeights = [];
            
            for fc = 1:nFC_methods
                currentFC_method = FC_methods{fc};
                currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
                load(currentResultsFile);
                
                edgeWeights = [edgeWeights, QCFC_parameters{8, i+1}'];
            end
            
            edgeWeightCorr = zeros(nFC_methods, nFC_methods); % matrix of correlation values between edge weights from different FC methods
            
            f = figure('Visible', 'off'); set(gcf, 'color', 'w');
            f.PaperUnits = 'inches';
            f.PaperPosition = [0 0 6 5];
            
            ctr = 1;
            
            for j = 1:nFC_methods
                for k = 1:nFC_methods
                    if j > k
                        ctr = ctr + 1;
                        continue;
                    else
                        subplot(nFC_methods, nFC_methods, ctr);
                        ax = gca;
                        ax.FontSize = 12;
                        if j == k
                            histogram(edgeWeights(:, j), 50, 'Normalization','probability');
                            ax.YTick = [];
                            if j < 3
                                ax.XLim = [-1, 1];
                                ax.XTick = [-1, 1];
                            else
                                ax.XLim = [0, 1];
                                ax.XTick = [0, 1];
                            end
                        else
                            edgeWeightCorr(k, j) = corr(edgeWeights(:, j), edgeWeights(:, k));
                            plot(edgeWeights(:, j), edgeWeights(:, k), '.', 'MarkerSize', 2);
                            if j < 3
                                ax.XLim = [-1, 1];
                                ax.XTick = [-1, 1];
                            else
                                ax.XLim = [0, 1];
                                ax.XTick = [0, 1];
                            end
                            
                            if k < 3
                                ax.YLim = [-1, 1];
                                ax.YTick = [-1, 1];
                            else
                                ax.YLim = [0, 1];
                                ax.YTick = [0, 1];
                            end
                        end
                        ctr = ctr + 1;
                    end
                end
            end
            currentSaveFileName = strcat(resultsFolder, 'plots/corrPlots_avgeEdgeWeights_', currentParcellation, '_', currentPipeline, '_', currentTaskType, '.svg');
            saveas(f, currentSaveFileName);
            close(f);
            
            g = figure('Visible', 'off'); set(gcf, 'color', 'w');
            g.PaperUnits = 'inches';
            g.PaperPosition = [0 0 6 5];
            im = imagesc(edgeWeightCorr, [0.4, 1]); im.AlphaData = tril(ones(nFC_methods, nFC_methods), -1); cmap = redbluecmap; cmap = cmap(6:end, :); colormap(cmap); colorbar; axis off;
            currentSaveFileName = strcat(resultsFolder, 'plots/corrHeatmaps_avgeEdgeWeights_', currentParcellation, '_', currentPipeline, '_', currentTaskType, '.svg');
            saveas(g, currentSaveFileName);
            close(g);
        end
    end
end

%% pairwise correlation plots between edge weights obtained using all 6 metrics, plotted for individual subjects

baseDir = '../Data/FunctionalConnectivityMatrices/';
resultsFolder = 'Results/02-Aug-2019/';
subjectID = '285345';

for p = 1:numel(parcellations)
    currentParcellation = parcellations{p};
    
    for i = 1:nPipelines
        currentPipeline = pipelines{i};
        
        for t = 1:numel(taskTypes)
            currentTaskType = taskTypes{t};
            edgeWeights = [];
            
            for fc = 1:nFC_methods
                currentFC_method = FC_methods{fc};
                currentResultsFile = strcat(baseDir, currentParcellation, '_', subjectID, '_', currentTaskType, '_', currentPipeline, '_', currentFC_method,  '.mat');
                load(currentResultsFile);
                
                edgeWeights = [edgeWeights, computeEdgeWeights(AdjMat)'];
            end
            
            f = figure('Visible', 'off'); set(gcf, 'color', 'w');
            f.PaperUnits = 'inches';
            f.PaperPosition = [0 0 6 5];
            
            ctr = 1;
            for j = 1:nFC_methods
                for k = 1:nFC_methods
                    if j > k
                        ctr = ctr + 1;
                        continue;
                    else
                        subplot(nFC_methods, nFC_methods, ctr);
                        ax = gca;
                        ax.FontSize = 12;
                        if j == k
                            histogram(edgeWeights(:, j), 50, 'Normalization','probability');
                            ax.YTick = [];
                            if j < 3
                                ax.XLim = [-1, 1];
                                ax.XTick = [-1, 1];
                            else
                                ax.XLim = [0, 1];
                                ax.XTick = [0, 1];
                            end
                        else
                            plot(edgeWeights(:, j), edgeWeights(:, k), '.', 'MarkerSize', 2);
                            if j < 3
                                ax.XLim = [-1, 1];
                                ax.XTick = [-1, 1];
                            else
                                ax.XLim = [0, 1];
                                ax.XTick = [0, 1];
                            end
                            
                            if k < 3
                                ax.YLim = [-1, 1];
                                ax.YTick = [-1, 1];
                            else
                                ax.YLim = [0, 1];
                                ax.YTick = [0, 1];
                            end
                        end
                        ctr = ctr + 1;
                    end
                end
            end
            currentSaveFileName = strcat(resultsFolder, 'plots/corrPlots_avgeEdgeWeights_', subjectID, '_', currentParcellation, '_', currentPipeline, '_', currentTaskType, '.svg');
            saveas(f, currentSaveFileName);
            close(f);
        end
    end
end

%% figure 2 - connectivity matrices for all FC metrics with sub-network labels

resultsFolder = 'Results/02-Aug-2019/';

currentParcellation = 'gordon';
M = importdata('../Data/GordonParcels/Parcels.xlsx');

% create community indices
nodeCommunityLabels = M.textdata(2:end, 5);
nodeCommunityIndices = zeros(size(nodeCommunityLabels, 1), 1);
uniqueLabels = unique(nodeCommunityLabels);

for i = 1:size(nodeCommunityLabels, 1)
    nodeCommunityIndices(i) = find(strcmp(uniqueLabels, nodeCommunityLabels{i}));
end

for t = 1:numel(taskTypes)
    currentTaskType = taskTypes{t};
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
        load(currentResultsFile);
        for i = 1:nPipelines
            currentPipeline = pipelines{i};
            edgeWeights_allEdges = QCFC_parameters{8, i+1};
            edgeWeights_matrix = squareform(edgeWeights_allEdges);
            
            [X,Y,INDSORT] = grid_communities(nodeCommunityIndices); % function from BCT to pool together communities for visualization
            
            f = figure('visible', 'off'); set(gcf, 'color', 'w');
            f.PaperUnits = 'inches';
            f.PaperPosition = [0 0 5 4];
            
            if strcmp(currentFC_method, 'Pearson') || strcmp(currentFC_method, 'Spearman')
                limits = [-0.5, 0.8];
            else
                limits = [0, 0.8];
            end
            
            imagesc(edgeWeights_matrix(INDSORT,INDSORT), limits);
            colormap('redbluecmap');
            colorbar;
            axis off;
            
            hold on;
            plot(X,Y,'k','linewidth',1);
            
            ax = gca;
            ax.FontSize = 20;
            legend off;
            currentSaveFileName = strcat(resultsFolder, 'plots/FC_matrices_', currentParcellation, '_', currentFC_method, '_', currentTaskType, '_', currentPipeline, '.svg');
            saveas(f, currentSaveFileName);
            close(f);
            
        end
    end
end

%% figure 3 - plot fraction of significant edges without FDR

resultsFolder = 'Results/02-Aug-2019/';

for p = 1:numel(parcellations)
    currentParcellation = parcellations{p};
    fractionSignificantEdges = cell(nPipelines, nFC_methods);
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        for t = 1:numel(taskTypes)
            currentTaskType = taskTypes{t};
            currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
            load(currentResultsFile);
            for i = 1:nPipelines
                fractionSignificantEdges{i, fc} = [fractionSignificantEdges{i, fc}; QCFC_parameters{3, i+1}];
            end
        end
    end
    
    for i = 1:nPipelines
        currentPipeline = pipelines{i};
        f = figure('Visible', 'off'); set(gcf, 'color', 'w');
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 6 5];
        
        x = []; y = []; x_idx = 1:1:nFC_methods;
        for fc = 1:nFC_methods
            y = [y; fractionSignificantEdges{i, fc}];
            x = [x; fc*ones(numel(fractionSignificantEdges{i, fc}), 1)];
        end
        
        beeswarm(x, 100*y, 'nosort', 'none', 2, 'sd');
        
        ax = gca;
        ax.FontSize = 20;
        if i==1
            ax.YLim = [0 100];
        elseif i==2
            ax.YLim = [0 50];
        end
        ax.XTickLabel = [];
        currentSaveFileName = strcat(resultsFolder, 'plots/fractionSignificantEdges_withoutFDR_', currentParcellation, '_', currentPipeline, '.svg');
        saveas(f, currentSaveFileName);
        close(f);
    end
end

%% figure 3 - histograms of QC-FC correlations

resultsFolder = 'Results/02-Aug-2019/';
cmap = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], ...
    [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560], ...
    [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330]};

for p = 1:numel(parcellations)
    currentParcellation = parcellations{p};
    for t = 1:numel(taskTypes)
        currentTaskType = taskTypes{t};
        for i = 1:nPipelines
            currentPipeline = pipelines{i};
            for fc = 1:nFC_methods
                currentFC_method = FC_methods{fc};
                currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
                load(currentResultsFile);
                
                QCFC_correlations = QCFC_parameters{6, i+1};
                
                f = figure('visible', 'off'); set(gcf, 'color', 'w');
                f.PaperUnits = 'inches';
                f.PaperPosition = [0 0 5 4];
                
                histogram(QCFC_correlations, 50, 'Normalization', 'probability', 'FaceColor', cmap{fc}, 'FaceAlpha', 0.9);
                
                ax = gca;
                ax.FontSize = 20;
                ax.YTick = [];
                Xlim([-0.3, 0.3]);
                currentSaveFileName = strcat(resultsFolder, 'plots/histogram_QCFC_correlations_', currentParcellation, '_', currentFC_method, '_', currentPipeline, '_', currentTaskType, '.svg');
                saveas(f, currentSaveFileName);
                close(f);
            end
        end
    end
end

%% supplementary figure x - scatter plots of edge weights versus QC-FC correlations

resultsFolder = 'Results/02-Aug-2019/';
cmap = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], ...
    [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560], ...
    [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330]};

for p = 1:numel(parcellations)
    currentParcellation = parcellations{p};
    for t = 1:numel(taskTypes)
        currentTaskType = taskTypes{t};
        for i = 1:nPipelines
            currentPipeline = pipelines{i};
            for fc = 1:nFC_methods
                currentFC_method = FC_methods{fc};
                currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
                load(currentResultsFile);
                
                QCFC_correlations = abs(QCFC_parameters{6, i+1});
                edgeWeights = QCFC_parameters{8, i+1};
                
                f = figure('visible', 'off'); set(gcf, 'color', 'w');
                f.PaperUnits = 'inches';
                f.PaperPosition = [0 0 5 4];
                
                scatter(edgeWeights, QCFC_correlations, 2, cmap{fc});
                lsline;
                
                ax = gca;
                ax.FontSize = 20;
                if strcmp(currentFC_method, 'Pearson') || strcmp(currentFC_method, 'Spearman')
                    ax.XLim = [-0.8, 0.8];
                else
                    ax.XLim = [0, 0.8];
                end
                ax.YLim = [0, 0.3];
                
                currentSaveFileName = strcat(resultsFolder, 'plots/edgeWeight_QCFC_correlations_', currentParcellation, '_', currentFC_method, '_', currentPipeline, '_', currentTaskType, '.svg');
                saveas(f, currentSaveFileName);
                close(f);
            end
        end
    end
end

%% supplementary figure x - plot fraction of significant edges with FDR

resultsFolder = 'Results/02-Aug-2019/';
for p = 1:numel(parcellations)
    currentParcellation = parcellations{p};
    fractionSignificantEdges = cell(nPipelines, nFC_methods);
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        for t = 1:numel(taskTypes)
            currentTaskType = taskTypes{t};
            currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
            load(currentResultsFile);
            for i = 1:nPipelines
                fractionSignificantEdges{i, fc} = [fractionSignificantEdges{i, fc}; QCFC_parameters{2, i+1}];
            end
        end
    end
    
    for i = 1:nPipelines
        currentPipeline = pipelines{i};
        f = figure('Visible', 'off'); set(gcf, 'color', 'w');
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 6 5];
        
        x = []; y = []; x_idx = 1:1:nFC_methods;
        for fc = 1:nFC_methods
            y = [y; fractionSignificantEdges{i, fc}];
            x = [x; fc*ones(numel(fractionSignificantEdges{i, fc}), 1)];
        end
        
        beeswarm(x, 100*y, 'nosort', 'none', 2, 'sd');
        
        ax = gca;
        ax.FontSize = 20;
        if i==1
            ax.YLim = [0 100];
        elseif i==2
            ax.YLim = [0 30];
        end
        ax.XTickLabel = [];
        currentSaveFileName = strcat(resultsFolder, 'plots/fractionSignificantEdges_withFDR_', currentParcellation, '_', currentPipeline, '.svg');
        saveas(f, currentSaveFileName);
        close(f);
    end
end

%% supplementary figure 2 - plot fraction of significant edges (absolute value) without FDR

resultsFolder = 'Results/02-Oct-2019/';
for p = 1:numel(parcellations)
    currentParcellation = parcellations{p};
    fractionSignificantEdges = cell(nPipelines, nFC_methods);
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        for t = 1:numel(taskTypes)
            currentTaskType = taskTypes{t};
            currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
            load(currentResultsFile);
            for i = 1:nPipelines
                fractionSignificantEdges{i, fc} = [fractionSignificantEdges{i, fc}; QCFC_parameters{5, i+1}];
            end
        end
    end
    
    for i = 1:nPipelines
        currentPipeline = pipelines{i};
        f = figure('Visible', 'off'); set(gcf, 'color', 'w');
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 6 5];
        
        x = []; y = []; x_idx = 1:1:nFC_methods;
        for fc = 1:nFC_methods
            y = [y; fractionSignificantEdges{i, fc}];
            x = [x; fc*ones(numel(fractionSignificantEdges{i, fc}), 1)];
        end
        
        beeswarm(x, 100*y, 'nosort', 'none', 2, 'sd');
        
        ax = gca;
        ax.FontSize = 20;
        if i==1
            ax.YLim = [0 100];
        elseif i==2
            ax.YLim = [0 50];
        end
        ax.XTickLabel = [];
        currentSaveFileName = strcat(resultsFolder, 'plots/fractionSignificantEdges_absoluteValue_withoutFDR_', currentParcellation, '_', currentPipeline, '.svg');
        saveas(f, currentSaveFileName);
        close(f);
    end
end

%% supplementary figure 3 - plot fraction of significant edges (zeroed out) without FDR

resultsFolder = 'Results/02-Oct-2019/';
for p = 1:numel(parcellations)
    currentParcellation = parcellations{p};
    fractionSignificantEdges = cell(nPipelines, nFC_methods);
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        for t = 1:numel(taskTypes)
            currentTaskType = taskTypes{t};
            currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
            load(currentResultsFile);
            for i = 1:nPipelines
                fractionSignificantEdges{i, fc} = [fractionSignificantEdges{i, fc}; QCFC_parameters{10, i+1}];
            end
        end
    end
    
    for i = 1:nPipelines
        currentPipeline = pipelines{i};
        f = figure('Visible', 'off'); set(gcf, 'color', 'w');
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 6 5];
        
        x = []; y = []; x_idx = 1:1:nFC_methods;
        for fc = 1:nFC_methods
            y = [y; fractionSignificantEdges{i, fc}];
            x = [x; fc*ones(numel(fractionSignificantEdges{i, fc}), 1)];
        end
        
        beeswarm(x, 100*y, 'nosort', 'none', 2, 'sd');
        
        ax = gca;
        ax.FontSize = 20;
        ax.YLim = [0 50];
        ax.XTickLabel = [];
        currentSaveFileName = strcat(resultsFolder, 'plots/fractionSignificantEdges_zeroedOut_withoutFDR_', currentParcellation, '_', currentPipeline, '.svg');
        saveas(f, currentSaveFileName);
        close(f);
    end
end

%% figure 4 - plot heatmaps of QC-FC correlations for all edges with network labels

resultsFolder = 'Results/02-Aug-2019/';

currentParcellation = 'gordon';
M = importdata('../Data/GordonParcels/Parcels.xlsx');

% create community indices
nodeCommunityLabels = M.textdata(2:end, 5);
nodeCommunityIndices = zeros(size(nodeCommunityLabels, 1), 1);
uniqueLabels = unique(nodeCommunityLabels);

for i = 1:size(nodeCommunityLabels, 1)
    nodeCommunityIndices(i) = find(strcmp(uniqueLabels, nodeCommunityLabels{i}));
end

for t = 1:numel(taskTypes)
    currentTaskType = taskTypes{t};
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
        load(currentResultsFile);
        for i = 1:nPipelines
            currentPipeline = pipelines{i};
            QCFC_correlations_allEdges = abs(QCFC_parameters{6, i+1}); % taking absolute value of QC-FC correlations
            QCFC_correlations_matrix = squareform(QCFC_correlations_allEdges);
            
            [X,Y,INDSORT] = grid_communities(nodeCommunityIndices); % function from BCT to pool together communities for visualization
            
            f = figure('visible', 'off'); set(gcf, 'color', 'w');
            f.PaperUnits = 'inches';
            f.PaperPosition = [0 0 5 4];
            imagesc(QCFC_correlations_matrix(INDSORT,INDSORT), [0, 0.3]);
            cmap = redbluecmap;
            cmap = cmap(6:end, :);
            colormap(cmap);
            colorbar;
            axis off;
            
            hold on;
            plot(X,Y,'k','linewidth',1);
            
            ax = gca;
            ax.FontSize = 20;
            legend off;
            currentSaveFileName = strcat(resultsFolder, 'plots/Heatmap_QCFC_correlations_', currentParcellation, '_', currentFC_method, '_', currentTaskType, '_', currentPipeline, '.svg');
            saveas(f, currentSaveFileName);
            close(f);
            
        end
    end
end

%% supplementary figure 5 - plot heatmaps of QC-FC correlations with mean values for each system with network labels

resultsFolder = 'Results/02-Aug-2019/';

currentParcellation = 'gordon';
M = importdata('../Data/GordonParcels/Parcels.xlsx');

% create community indices
nodeCommunityLabels = M.textdata(2:end, 5);
nodeCommunityIndices = zeros(size(nodeCommunityLabels, 1), 1);
uniqueLabels = unique(nodeCommunityLabels);
nUniqueLabels = numel(uniqueLabels);

for i = 1:size(nodeCommunityLabels, 1)
    nodeCommunityIndices(i) = find(strcmp(uniqueLabels, nodeCommunityLabels{i}));
end

for t = 1:numel(taskTypes)
    currentTaskType = taskTypes{t};
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
        load(currentResultsFile);
        for i = 1:nPipelines
            QCFC_correlations_matrix_mean = zeros(nUniqueLabels, nUniqueLabels); % creating a matrix of mean QCFC values for each system-system pair
            currentPipeline = pipelines{i};
            QCFC_correlations_allEdges = abs(QCFC_parameters{6, i+1}); % taking absolute value of QC-FC correlations
            QCFC_correlations_matrix = squareform(QCFC_correlations_allEdges);
            for j = 1:nUniqueLabels
                for k = 1:nUniqueLabels
                    QCFC_correlations_matrix_mean(j, k) = mean2(QCFC_correlations_matrix(nodeCommunityIndices==j, nodeCommunityIndices==k));
                end
            end
                    
            f = figure('visible', 'off'); set(gcf, 'color', 'w');
            f.PaperUnits = 'inches';
            f.PaperPosition = [0 0 5 4];
            %imagesc(QCFC_correlations_matrix_mean, [0, 0.3]);
            imagesc(QCFC_correlations_matrix_mean, [0, 0.08]);
            cmap = redbluecmap;
            cmap = cmap(6:end, :);
            colormap(cmap);
            colorbar;
            
            ax = gca;
            ax.FontSize = 8;
            ax.XTick = 1:nUniqueLabels;            
            ax.YTick = 1:nUniqueLabels;
            ax.XTickLabel = uniqueLabels;            
            ax.YTickLabel = uniqueLabels;
            ax.XTickLabelRotation = 90;
            legend off;
            currentSaveFileName = strcat(resultsFolder, 'plots/Heatmap_QCFC_correlations_meanSystem_', currentParcellation, '_', currentFC_method, '_', currentTaskType, '_', currentPipeline, '.svg');
            saveas(f, currentSaveFileName);
            close(f);
            
        end
    end
end

%% figure 5 and supplementary figure 6- compute all pairwise inter-community QC-FC correlations and plot the top ranking values for each FC metric as boxplots

resultsFolder = 'Results/02-Aug-2019/';

currentParcellation = 'gordon';
nNodes = 333;
M = importdata('../Data/GordonParcels/Parcels.xlsx');
nPairings2Plot = 6;

% create community indices
nodeCommunityLabels = M.textdata(2:end, 5);
nodeCommunityIndices = zeros(size(nodeCommunityLabels, 1), 1);
uniqueLabels = unique(nodeCommunityLabels);
nCommunities = numel(uniqueLabels);

uniqueLabels_abbreviated = cell(size(uniqueLabels));
for i = 1:nCommunities
    currentCommunity = uniqueLabels{i};
    switch currentCommunity
        case 'None'
            uniqueLabels_abbreviated{i} = 'N';
        case 'Default'
            uniqueLabels_abbreviated{i} = 'D';
        case 'Visual'
            uniqueLabels_abbreviated{i} = 'V';
        case 'Auditory'
            uniqueLabels_abbreviated{i} = 'A';
        case 'SMhand'
            uniqueLabels_abbreviated{i} = 'SH';
        case 'SMmouth'
            uniqueLabels_abbreviated{i} = 'SM';
        case 'DorsalAttn'
            uniqueLabels_abbreviated{i} = 'DA';
        case 'VentralAttn'
            uniqueLabels_abbreviated{i} = 'VA';
        case 'CinguloOperc'
            uniqueLabels_abbreviated{i} = 'CO';
        case 'Salience'
            uniqueLabels_abbreviated{i} = 'S';
        case 'FrontoParietal'
            uniqueLabels_abbreviated{i} = 'FP';
        case 'CinguloParietal'
            uniqueLabels_abbreviated{i} = 'CP';
        case 'RetrosplenialTemporal'
            uniqueLabels_abbreviated{i} = 'RT';
    end
end

for i = 1:size(nodeCommunityLabels, 1)
    nodeCommunityIndices(i) = find(strcmp(uniqueLabels, nodeCommunityLabels{i}));
end

% % assigning each community a color for plotting
% uniqueLabels_cmap = zeros(nCommunities, 3);
% for i = 1:nCommunities
%     currentCommunity = uniqueLabels{i};
%     switch currentCommunity
%         case 'None'
%             uniqueLabels_cmap(i, :) = [230/255 230/255 230/255]; % white
%         case 'Default'
%             uniqueLabels_cmap(i, :) = [255/255 0 0]; % red
%         case 'Visual'
%             uniqueLabels_cmap(i, :) = [0 0 255/255]; % blue
%         case 'Auditory'
%             uniqueLabels_cmap(i, :) = [255/255 100/255 100/255]; % salmon pink
%         case 'SMhand'
%             uniqueLabels_cmap(i, :) = [0 255/255 255/255]; % cyan
%         case 'SMmouth'
%             uniqueLabels_cmap(i, :) = [255/255 100/255 0]; % orange
%         case 'DorsalAttn'
%             uniqueLabels_cmap(i, :) = [0 100/255 50/255]; % dark green
%         case 'VentralAttn'
%             uniqueLabels_cmap(i, :) = [0 166/255 156/255]; % turquoise
%         case 'CinguloOperc'
%             uniqueLabels_cmap(i, :) = [101/255 45/255 144/255]; % purple
%         case 'Salience'
%             uniqueLabels_cmap(i, :) = [0 0 0]; % black
%         case 'FrontoParietal'
%             uniqueLabels_cmap(i, :) = [255/255 255/255 0]; % yellow
%         case 'CinguloParietal'
%             uniqueLabels_cmap(i, :) = [235/255 0 139/255]; % magenta
%         case 'RetrosplenialTemporal'
%             uniqueLabels_cmap(i, :) = [0 165/255 81/255]; % green
%     end
% end

for t = 1:numel(taskTypes)
    currentTaskType = taskTypes{t};
    for i = 1:nPipelines
        currentPipeline = pipelines{i};
        for fc = 1:nFC_methods
            currentFC_method = FC_methods{fc};
            currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
            load(currentResultsFile);
            QCFC_correlations_allEdges = QCFC_parameters{6, i+1}; % taking absolute value of QC-FC correlations
            QCFC_correlations_matrix = squareform(QCFC_correlations_allEdges);
            
            % storing all inter-community QC-FC edges in a 13x13 cell array
            interCommunityEdges = cell(nCommunities, nCommunities);
            for j = 1:nNodes
                community_node1 = nodeCommunityIndices(j);
                for k = j+1:nNodes
                    community_node2 = nodeCommunityIndices(k);
                    interCommunityEdges{community_node1, community_node2} = [interCommunityEdges{community_node1, community_node2}; QCFC_correlations_matrix(j, k)];
                end
            end
            
            %% plotting rankings of community pairs by median inter-network QCFC correlations
            totalCombinations = nCommunities*(nCommunities-1)/2 + nCommunities; % total possible combinations of communities (including intra-community)
            interCommunityEdges_flattened = cell(1, totalCombinations);
            interCommunityEdges_labels = cell(1, totalCombinations);
            interCommunityEdges_median = zeros(1, totalCombinations);
            ctr = 1;
            for j = 1:nCommunities
                for k = j:nCommunities
                    interCommunityEdges_flattened{ctr} = [interCommunityEdges{j, k}; interCommunityEdges{k, j}]; % including edges going from j to k and k to j
                    interCommunityEdges_labels{ctr} = strcat(uniqueLabels_abbreviated{j}, ' - ', uniqueLabels_abbreviated{k});
                    interCommunityEdges_median(ctr) = median([interCommunityEdges{j, k}; interCommunityEdges{k, j}]);
                    ctr = ctr + 1;
                end
            end
            
            % storing top ranking inter-community QC-FC correlations
            [sortedMedians, idx] = sort(interCommunityEdges_median, 'descend');
            interCommunityEdges_toPlot = [];
            interCommunityEdges_toPlot_groupings = [];
            interCommunityEdges_toPlot_labels = {};
            for j = 1:nPairings2Plot
                current_idx = idx(j);
                interCommunityEdges_toPlot = [interCommunityEdges_toPlot; interCommunityEdges_flattened{current_idx}];
                interCommunityEdges_toPlot_groupings = [interCommunityEdges_toPlot_groupings; j*ones(size(interCommunityEdges_flattened{current_idx}))];
                interCommunityEdges_toPlot_labels{j} = interCommunityEdges_labels{current_idx};
            end
            
            f = figure('visible', 'off'); set(gcf, 'color', 'w');
            %             f.PaperUnits = 'inches';
            %             f.PaperPosition = [0 0 5 4];
            
            boxplot(interCommunityEdges_toPlot, interCommunityEdges_toPlot_groupings, 'Labels', interCommunityEdges_toPlot_labels, 'LabelOrientation', 'inline');
            legend off;
            ax = gca;
            ax.FontSize = 20;
            ylim([0, 0.3]);
            currentSaveFileName = strcat(resultsFolder, 'plots/interCommunity_QCFC_correlations_', currentParcellation, '_', currentFC_method, '_', currentTaskType, '_', currentPipeline, '.svg');
            saveas(f, currentSaveFileName);
            close(f);
                   
            %% plotting rankings of communities by median inter- and intra-network QCFC correlations
            individualCommunityEdges = cell(1, nCommunities);
            individualCommunityEdges_labels = cell(1, nCommunities);
            individualCommunityEdges_median = zeros(1, nCommunities);
            for j = 1:nCommunities
                allEdges_currentCommunity = interCommunityEdges{j, :};
                individualCommunityEdges{j} = allEdges_currentCommunity; % including all edges from current network
                individualCommunityEdges_labels{j} = uniqueLabels_abbreviated{j};
                individualCommunityEdges_median(j) = median(allEdges_currentCommunity);
            end
            
            % ranking communities by median QCFC correlation
            [sortedMedians, idx] = sort(individualCommunityEdges_median, 'descend');
            individualCommunityEdges_toPlot = [];
            individualCommunityEdges_toPlot_groupings = [];
            individualCommunityEdges_toPlot_labels = {};
            for j = 1:nCommunities
                current_idx = idx(j);
                individualCommunityEdges_toPlot = [individualCommunityEdges_toPlot; individualCommunityEdges{current_idx}];
                individualCommunityEdges_toPlot_groupings = [individualCommunityEdges_toPlot_groupings; j*ones(size(individualCommunityEdges{current_idx}))];
                individualCommunityEdges_toPlot_labels{j} = individualCommunityEdges_labels{current_idx};
            end
            
            f = figure('visible', 'off'); set(gcf, 'color', 'w');
            %             f.PaperUnits = 'inches';
            %             f.PaperPosition = [0 0 5 4];
            
            boxplot(individualCommunityEdges_toPlot, individualCommunityEdges_toPlot_groupings, 'Labels', individualCommunityEdges_toPlot_labels, 'LabelOrientation', 'inline');
            legend off;
            ax = gca;
            ax.FontSize = 20;
            %ylim([0, 0.3]);
            currentSaveFileName = strcat(resultsFolder, 'plots/individualCommunity_QCFC_correlations_', currentParcellation, '_', currentFC_method, '_', currentTaskType, '_', currentPipeline, '.svg');
            saveas(f, currentSaveFileName);
            close(f);
            
            %% plotting rankings of communities by median intra-network QCFC correlations only
            individualCommunityEdges = cell(1, nCommunities);
            individualCommunityEdges_labels = cell(1, nCommunities);
            individualCommunityEdges_median = zeros(1, nCommunities);
            for j = 1:nCommunities
                allEdges_currentCommunity = interCommunityEdges{j, j};
                individualCommunityEdges{j} = allEdges_currentCommunity; % including all edges from current network
                individualCommunityEdges_labels{j} = uniqueLabels_abbreviated{j};
                individualCommunityEdges_median(j) = median(allEdges_currentCommunity);
            end
            
            % ranking communities by median QCFC correlation
            [sortedMedians, idx] = sort(individualCommunityEdges_median, 'descend');
            individualCommunityEdges_toPlot = [];
            individualCommunityEdges_toPlot_groupings = [];
            individualCommunityEdges_toPlot_labels = {};
            for j = 1:nCommunities
                current_idx = idx(j);
                individualCommunityEdges_toPlot = [individualCommunityEdges_toPlot; individualCommunityEdges{current_idx}];
                individualCommunityEdges_toPlot_groupings = [individualCommunityEdges_toPlot_groupings; j*ones(size(individualCommunityEdges{current_idx}))];
                individualCommunityEdges_toPlot_labels{j} = individualCommunityEdges_labels{current_idx};
            end
            
            f = figure('visible', 'off'); set(gcf, 'color', 'w');
            %             f.PaperUnits = 'inches';
            %             f.PaperPosition = [0 0 5 4];
            
            boxplot(individualCommunityEdges_toPlot, individualCommunityEdges_toPlot_groupings, 'Labels', individualCommunityEdges_toPlot_labels, 'LabelOrientation', 'inline');
            legend off;
            ax = gca;
            ax.FontSize = 20;
            %ylim([0, 0.3]);
            currentSaveFileName = strcat(resultsFolder, 'plots/individualCommunity_intraNetwork_QCFC_correlations_', currentParcellation, '_', currentFC_method, '_', currentTaskType, '_', currentPipeline, '.svg');
            saveas(f, currentSaveFileName);
            close(f);
        end
    end
end

%% figure 6 - plot distance-dependence of QC-FC correlations

resultsFolder = 'Results/02-Aug-2019/';

for p = 1:numel(parcellations)
    currentParcellation = parcellations{p};
    
    % extracting centroids and computing pairwise node distances for
    % different parcellations
    switch currentParcellation
        case 'gordon'
            path = '/Users/ArunMahadevan/Documents/hcp_Max/Data/GordonParcels/Parcels.xlsx';
            distanceData = readtable(path);
            centroids = distanceData.Centroid_MNI_;
            nNodes = size(centroids, 1);
            centroidsNumbers = zeros(nNodes, 3);
            for i = 1:nNodes
                currentCentroid = convertCharsToStrings(centroids(i));
                spaceLocations = strfind(currentCentroid, ' ');
                centroidsNumbers(i, 1) = str2double(extractBetween(currentCentroid, 1, spaceLocations(1)-1));
                centroidsNumbers(i, 2) = str2double(extractBetween(currentCentroid, spaceLocations(1)+1, spaceLocations(2)-1));
                centroidsNumbers(i, 3) = str2double(extractBetween(currentCentroid, spaceLocations(2)+1, strlength(currentCentroid)));
            end
            
            nodeDistanceMatrix = pdist(centroidsNumbers)';
            
        case 'yeo_100'
            path = '/Users/ArunMahadevan/Documents/hcp_Max/Data/SchaeferParcels/MNI/Centroid_coordinates/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv'; % using 17-network parcellations in 2mm MNI space
            centroidsNumbers = csvread(path, 2, 2);
            nNodes = size(centroidsNumbers, 1);
            nodeDistanceMatrix = pdist(centroidsNumbers)';
    end
    
    QCFC_distance_correlations = cell(nPipelines, nFC_methods);
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        for t = 1:numel(taskTypes)
            currentTaskType = taskTypes{t};
            currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
            load(currentResultsFile);
            for i = 1:nPipelines
                QCFC_correlations = abs(QCFC_parameters{6, i+1})'; % taking absolute value of QC-FC correlations
                QCFC_distance_correlations{i, fc} = [QCFC_distance_correlations{i, fc}; corr(QCFC_correlations, nodeDistanceMatrix)];
            end
        end
    end
    
    for i = 1:nPipelines
        currentPipeline = pipelines{i};
        f = figure('visible', 'off'); set(gcf, 'color', 'w');
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 6 5];
        
        x = []; y = []; x_idx = 1:1:nFC_methods;
        for fc = 1:nFC_methods
            y = [y; QCFC_distance_correlations{i, fc}];
            x = [x; fc*ones(numel(QCFC_distance_correlations{i, fc}), 1)];
        end
        
        beeswarm(x, y, 'nosort', 'none', 2, 'sd');
        
        ax = gca;
        ax.FontSize = 20;
        if i==1
            ax.YLim = [-0.1 0.1];
        elseif i==2
            ax.YLim = [-0.05 0.11];
        end
        ax.XTickLabel = [];
        currentSaveFileName = strcat(resultsFolder, 'plots/QCFC_distanceDependence_', currentParcellation, '_', currentPipeline, '.svg');
        saveas(f, currentSaveFileName);
        close(f);
    end
end

%% figure 7 - plot intra-class correlation for all edges and for edges unaffected by motion

resultsFolder = 'Results/15-Aug-2019/';

for p = 1:numel(parcellations)
    currentParcellation = parcellations{p};
    
    %for i = 1:nPipelines
    %    currentPipeline = pipelines{i};
    currentPipeline = 'FIX_matrices';
    
    ICC_allFCmethods = []; ICC_noMotionEdges_allFCmethods = [];
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        
        currentResultsFile = strcat(resultsFolder, 'IntraClassCorrelation_edgeWeights', currentFC_method, '_', currentParcellation, '_', currentPipeline, '.mat');
        load(currentResultsFile);
        
        ICC = ICC_stats_FCedgeWeights{2, 1};
        ICC_noMotionEdges = ICC_stats_FCedgeWeights{2, 2};
        
        ICC_allFCmethods(:, fc) = ICC;
        ICC_noMotionEdges_allFCmethods(:, fc) = ICC_noMotionEdges;
    end
    
    f = figure('visible', 'off'); f.PaperUnits = 'inches'; f.PaperPosition = [0 0 6 5];
    set(gcf, 'color', 'w');
    violin(ICC_allFCmethods, 'facecolor',[0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
        0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; ...
        0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330], 'facealpha', 1, 'medc','k','mc','');
    ax = gca; ax.FontSize = 20; ax.YLim = [-0.25 1]; ax.XTickLabel = '';
    legend off;
    currentSaveFileName = strcat(resultsFolder, 'IntraClassCorrelation_', currentParcellation, '_', currentPipeline, '.svg');
    saveas(f, currentSaveFileName);
    close(f);
    
    f = figure('visible', 'off'); f.PaperUnits = 'inches'; f.PaperPosition = [0 0 6 5];
    set(gcf, 'color', 'w');
    violin(ICC_noMotionEdges_allFCmethods, 'facecolor',[0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
        0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; ...
        0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330], 'facealpha', 1, 'medc','k','mc','');
    ax = gca; ax.FontSize = 20; ax.YLim = [-0.25 1]; ax.XTickLabel = '';
    legend off;
    currentSaveFileName = strcat(resultsFolder, 'IntraClassCorrelation_noMotionEdges_', currentParcellation, '_', currentPipeline, '.svg');
    saveas(f, currentSaveFileName);
    close(f);
    %end
end

%% figure x - connectivity heatmaps for topological matrices with sub-network labels

resultsFolder = 'Results/02-Oct-2019/';

currentParcellation = 'gordon';
M = importdata('../Data/GordonParcels/Parcels.xlsx');

% create community indices
nodeCommunityLabels = M.textdata(2:end, 5);
nodeCommunityIndices = zeros(size(nodeCommunityLabels, 1), 1);
uniqueLabels = unique(nodeCommunityLabels);

for i = 1:size(nodeCommunityLabels, 1)
    nodeCommunityIndices(i) = find(strcmp(uniqueLabels, nodeCommunityLabels{i}));
end

[X,Y,INDSORT] = grid_communities(nodeCommunityIndices); % function from BCT to pool together communities for visualization

for t = 1:numel(taskTypes)
    currentTaskType = taskTypes{t};
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_topologicalSimilarityMatrix_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
        load(currentResultsFile);
        
        currentPipeline = 'FIX_matrices';
        edgeWeights_correlation_allEdges = QCFC_parameters_topologicalSimilarityMatrix{7, 2};
        edgeWeights_correlation_matrix = squareform(edgeWeights_correlation_allEdges);
        
        f = figure('visible', 'off'); set(gcf, 'color', 'w');
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 5 4];
        
        imagesc(edgeWeights_correlation_matrix(INDSORT,INDSORT), [-1, 1]);
        colormap('redbluecmap');
        colorbar;
        axis off;
        
        hold on;
        plot(X,Y,'k','linewidth',1);
        
        ax = gca;
        ax.FontSize = 20;
        legend off;
        currentSaveFileName = strcat(resultsFolder, 'plots/topologicalMatrices_', currentParcellation, '_', currentFC_method, '_', currentTaskType, '_', currentPipeline, '.svg');
        saveas(f, currentSaveFileName);
        close(f);
    end
end

%% figure x - plot fraction of significant edges in topological matrices without FDR

resultsFolder = 'Results/02-Oct-2019/';

for p = 1:numel(parcellations)
    currentParcellation = parcellations{p};
    fractionSignificantEdges = cell(1, nFC_methods);
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        for t = 1:numel(taskTypes)
            currentTaskType = taskTypes{t};
            currentResultsFile = strcat(resultsFolder, 'QCFC_parameters_topologicalSimilarityMatrix_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
            load(currentResultsFile);
            
            fractionSignificantEdges{1, fc} = [fractionSignificantEdges{1, fc}; QCFC_parameters_topologicalSimilarityMatrix{3, 2}];
        end
    end
    
    currentPipeline = 'FIX_matrices';
    f = figure('Visible', 'off'); set(gcf, 'color', 'w');
    f.PaperUnits = 'inches';
    f.PaperPosition = [0 0 6 5];
    
    x = []; y = []; x_idx = 1:1:nFC_methods;
    for fc = 1:nFC_methods
        y = [y; fractionSignificantEdges{1, fc}];
        x = [x; fc*ones(numel(fractionSignificantEdges{1, fc}), 1)];
    end
    
    beeswarm(x, 100*y, 'nosort', 'none', 2, 'sd');
    
    ax = gca;
    ax.FontSize = 20;
    ax.YLim = [0 50];
    ax.XTickLabel = [];
    currentSaveFileName = strcat(resultsFolder, 'plots/fractionSignificantEdges_topologicalMatrices_', currentParcellation, '_', currentPipeline, '.svg');
    saveas(f, currentSaveFileName);
    close(f);
    
end