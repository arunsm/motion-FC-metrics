clear all;
close all;
taskTypes = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
FC_methods = {'Pearson', 'Spearman', 'MutualInformationTime', 'Coherence', 'WaveletCoherence', 'MutualInformation'}; nFC_methods = numel(FC_methods);
parcellationTypes = {'gordon', 'yeo_100'};
pipelines = {'CompCor_matrices', 'FIX_matrices'}; nPipelines = numel(pipelines);

labels = FC_methods;

addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/violin'));
addpath(genpath('/Users/ArunMahadevan/Documents/MATLAB/BCT'));

%% plot fraction of significant edges

resultsFolder = 'Results/03-Jun-2019/';

for p = 1:numel(parcellationTypes)
    currentParcellation = parcellationTypes{p};
    fractionSignificantEdges_numbers = zeros(nPipelines, nFC_methods);
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        for t = 1:numel(taskTypes)
            currentTaskType = taskTypes{t};
            currentResultsFile = strcat(resultsFolder, 'fractionSignificantEdges_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
            load(currentResultsFile);
            for i = 1:nPipelines
                fractionSignificantEdges_numbers(i, fc) = fractionSignificantEdges_numbers(i, fc) + fractionSignificantEdges{2, i};
            end
        end
    end
    fractionSignificantEdges_numbers = fractionSignificantEdges_numbers/numel(taskTypes); % averaging across resting state scans
    
    %         [sortedFractonSignificantEdges, idx] = sort(fractionSignificantEdges_numbers, 'descend');
    %         labels = FC_methods(idx);
    
    for i = 1:nPipelines
        currentPipeline = pipelines{i};
        f = figure('visible', 'off'); set(gcf, 'color', 'w');
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 6 5];
        x = 1:nFC_methods;
        bar(100*fractionSignificantEdges_numbers(i, :), 'c', 'EdgeColor', 'c');
        for j = 1:nFC_methods
            h = text(x(j), 8, num2str(100*fractionSignificantEdges_numbers(i, j), '%0.2f'), 'FontSize', 18, 'HorizontalAlignment', 'center');
            set(h,'Rotation',90);
        end
        ax = gca;
        ax.FontSize = 20;
        ax.YLim = [0 80];
        ax.XTickLabel = [];
        %ax.XTickLabelRotation = 90;
        %ylabel("Edges related to motion (%)");
        %title(strcat(currentTaskType, '-', currentParcellation, '-', pipelines(i)));
        %currentSaveFileName = strcat(resultsFolder, currentTaskType, '_', currentParcellation, '_', currentPipeline, '.pdf');
        %saveas(f, currentSaveFileName, 'pdf');
        currentSaveFileName = strcat(resultsFolder, 'fractionSignificantEdges_', currentParcellation, '_', currentPipeline, '.svg');
        saveas(f, currentSaveFileName);
        close(f);
    end
end

%% plot median absolute correlation

resultsFolder = 'Results/03-Jun-2019/';

for p = 1:numel(parcellationTypes)
    currentParcellation = parcellationTypes{p};
    for t = 1:numel(taskTypes)
        currentTaskType = taskTypes{t};
        medianAbsoluteCorrelation_numbers = zeros(nPipelines, nFC_methods);
        for fc = 1:nFC_methods
            currentFC_method = FC_methods{fc};
            currentResultsFile = strcat(resultsFolder, 'medianAbsoluteCorrelation_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
            load(currentResultsFile);
            for i = 1:nPipelines
                medianAbsoluteCorrelation_numbers(i, fc) = medianAbsoluteCorrelation{2, i};
            end
        end
        
        %         [sortedFractonSignificantEdges, idx] = sort(fractionSignificantEdges_numbers, 'descend');
        %         labels = FC_methods(idx);
        
        for i = 1:nPipelines
            currentPipeline = pipelines{i};
            f = figure('visible', 'off'); set(gcf, 'color', 'w');
            x = 1:nFC_methods;
            bar(medianAbsoluteCorrelation_numbers(i, :), 'c', 'EdgeColor', 'c');
            for j = 1:nFC_methods
                h = text(x(j), 0.01, num2str(medianAbsoluteCorrelation_numbers(i, j), '%0.2f'), 'FontSize', 18, 'HorizontalAlignment', 'center');
                set(h,'Rotation',90);
            end
            ax = gca;
            ax.FontSize = 20;
            ax.YLim = [0 0.1];
            %ax.XTickLabel = labels;
            ax.XTickLabel = [];
            %ax.XTickLabelRotation = 90;
            %ylabel("Median absolute correlation");
            %title(strcat(currentTaskType, '-', currentParcellation, '-', pipelines(i)));
            %currentSaveFileName = strcat(resultsFolder, currentTaskType, '_', currentParcellation, '_', currentPipeline, '.pdf');
            %saveas(f, currentSaveFileName, 'pdf');
            currentSaveFileName = strcat(resultsFolder, 'medianAbsoluteCorrelation_', currentTaskType, '_', currentParcellation, '_', currentPipeline, '.png');
            saveas(f, currentSaveFileName);
            close(f);
        end
    end
end

%% plot intra-class correlation for edge weights

resultsFolder = 'Results/03-Jun-2019/';

for p = 1:numel(parcellationTypes)
    currentParcellation = parcellationTypes{p};
    
    for i = 1:nPipelines
        currentPipeline = pipelines{i};
        
        ICC_allFCmethods = [];
        for fc = 1:nFC_methods
            currentFC_method = FC_methods{fc};
            
            currentResultsFile = strcat(resultsFolder, 'IntraClassCorrelation_edgeWeights', currentFC_method, '_', currentParcellation, '_', currentPipeline, '.mat');
            load(currentResultsFile);
            
            ICC = ICC_stats_FCedgeWeights{2, 1};
            ICC_allFCmethods(:, fc) = ICC;
        end
        
        %f = figure('visible', 'off');
        f = figure;
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 6 5];
        set(gcf, 'color', 'w');
        violin(ICC_allFCmethods, 'medc','k','mc','');
        ax = gca;
        ax.FontSize = 20;
        ax.YLim = [-0.25 1];
        ax.XTickLabel = labels;
        ax.XTickLabel = [];
        ax.XTickLabelRotation = 90;
        %ylabel("Intra-class correlation (edge weights)");
        %title(strcat(currentTaskType, '-', currentParcellation, '-', pipelines(i)));
        %currentSaveFileName = strcat(resultsFolder, currentTaskType, '_', currentParcellation, '_', currentPipeline, '.pdf');
        %saveas(f, currentSaveFileName, 'pdf');
        legend off;
        currentSaveFileName = strcat(resultsFolder, 'IntraClassCorrelation_', currentParcellation, '_', currentPipeline, '.svg');
        saveas(f, currentSaveFileName);
        close(f);
    end
end

%% plot heatmaps of QC-FC correlations for all edges with network labels

resultsFolder = 'Results/03-Jun-2019/';

currentParcellation = 'gordon';
nNodes = 333;
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
        currentResultsFile = strcat(resultsFolder, 'medianAbsoluteCorrelation_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
        load(currentResultsFile);
        for i = 1:nPipelines
            currentPipeline = pipelines{i};
            QCFC_correlations_allEdges = medianAbsoluteCorrelation{4, i};
            QCFC_correlations_matrix = squareform(QCFC_correlations_allEdges);
            [X,Y,INDSORT] = grid_communities(nodeCommunityIndices); % function from BCT to pool together communities for visualization
            
            % find average QC-FC correlation for each node
            avgeCorrelation_allNodes = zeros(1, nNodes);
            for n = 1:nNodes
                avgeCorrelation_allNodes(n) = mean(QCFC_correlations_allEdges(:, n)); % averaging across all edges to current node
            end
            
            currentSaveFileName = strcat(resultsFolder, 'avgeCorrelation_allNodes_', currentParcellation, '_', currentFC_method, '_', currentTaskType, '_', currentPipeline, '.mat');
            save(currentSaveFileName, 'avgeCorrelation_allNodes');
            
            f = figure('visible', 'off');
            f.PaperUnits = 'inches';
            f.PaperPosition = [0 0 5 4];
            %f = figure;
            set(gcf, 'color', 'w');
            imagesc(QCFC_correlations_matrix(INDSORT,INDSORT), [-0.3, 0.3]);
            %imagesc(QCFC_correlations_matrix(INDSORT,INDSORT));
            colormap('redbluecmap');
            colorbar;
            axis off;
            
            hold on;
            plot(X,Y,'k','linewidth',1);
            
            ax = gca;
            ax.FontSize = 20;
            legend off;
            currentSaveFileName = strcat(resultsFolder, 'Heatmap_QCFC_correlations_', currentParcellation, '_', currentFC_method, '_', currentTaskType, '_', currentPipeline, '.svg');
            saveas(f, currentSaveFileName);
            close(f);
        end
    end
    
    %         [sortedFractonSignificantEdges, idx] = sort(fractionSignificantEdges_numbers, 'descend');
    %         labels = FC_methods(idx);
    
end

%% plot distance-dependence of QC-FC correlations

resultsFolder = 'Results/03-Jun-2019/';

for p = 1:numel(parcellationTypes)
    currentParcellation = parcellationTypes{p};
    
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
            fprintf('boo')
            continue;
    end
    
    QCFC_distance_correlations = zeros(nPipelines, nFC_methods);
    for fc = 1:nFC_methods
        currentFC_method = FC_methods{fc};
        for t = 1:numel(taskTypes)
            currentTaskType = taskTypes{t};
            currentResultsFile = strcat(resultsFolder, 'medianAbsoluteCorrelation_', currentFC_method, '_', currentParcellation, '_', currentTaskType, '_allPipelines.mat');
            load(currentResultsFile);
            for i = 1:nPipelines
                QCFC_correlations = medianAbsoluteCorrelation{4, i}';
                QCFC_distance_correlations(i, fc) = QCFC_distance_correlations(i, fc) + corr(QCFC_correlations, nodeDistanceMatrix);
            end
        end
    end
    
    QCFC_distance_correlations = QCFC_distance_correlations/numel(taskTypes); % averaging across resting state scans
    
    for i = 1:nPipelines
        currentPipeline = pipelines{i};
        f = figure('visible', 'off'); set(gcf, 'color', 'w');
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 6 5];
        x = 1:nFC_methods;
        bar(QCFC_distance_correlations(i, :), 'c', 'EdgeColor', 'c');
        for j = 1:nFC_methods
            h = text(x(j), 0.02, num2str(QCFC_distance_correlations(i, j), '%0.2f'), 'FontSize', 18, 'HorizontalAlignment', 'center');
            set(h,'Rotation',90);
        end
        ax = gca;
        ax.FontSize = 20;
        ax.YLim = [-0.15 0.1];
        ax.XTickLabel = [];
        %ax.XTickLabelRotation = 90;
        %ylabel("Edges related to motion (%)");
        %title(strcat(currentTaskType, '-', currentParcellation, '-', pipelines(i)));
        %currentSaveFileName = strcat(resultsFolder, currentTaskType, '_', currentParcellation, '_', currentPipeline, '.pdf');
        %saveas(f, currentSaveFileName, 'pdf');
        currentSaveFileName = strcat(resultsFolder, 'QCFC_distanceDependence_', currentParcellation, '_', currentPipeline, '.svg');
        saveas(f, currentSaveFileName);
        close(f);
    end
end