%% THIS SCRIPT RUNS MODULARITY ANALYSIS ON 4 RUNS of EACH SUBJECT, LOW-MOTION EDGES ONLY, USING YEO100 or GORDON PARCELLATION
% ICA-FIX MATRICES ONLY

%% SETUP
%working directly on the cluster
subj_dir='/data/jag/bassett-lab/hcp_Max/Data/Covariates'
data_dir='/data/jag/bassett-lab/hcp_Max/Data/QCFC_correlationMatrices'
outdir='/data/jux/mackey_group/Ursula/projects/in_progress/arun_fc_metrics_motion/output/data/Gordon_ICA_FIX'
%subject list
subjList=readtable(fullfile(subj_dir, 'S1200_Release_Subjects_Demographics.csv'));
subjList=subjList.Subject;

rest_runs={'_REST1_LR_', '_REST1_RL_','_REST2_LR_','_REST2_RL_'};
fc_metrics={'Coherence', 'MutualInformation', 'MutualInformationTime','Pearson','Spearman', 'WaveletCoherence'};
%% initialize vectors
%%%%%
%for each of four runs
%%%%%%
%find the edges across all 6 metrics that are the least overall affected by
%motion
for i=1:4
    matrix=zeros(100,100);
    run=rest_runs{i};
    run_name=strcat('run',run, 'l');
for j=1:6 %only Pearson and Spearman
    metric=fc_metrics{j}
    %read in the average matrix, take the absolute value of correlation,
    %and rank it
    file=fullfile(data_dir,strcat('QCFC_correlationMatrix_',metric,'_gordon',run,'FIX_matrices.mat'));
    load(file);
    QCFC_correlationMatrix_abs=abs(QCFC_correlationMatrix);
    [~, ranks]=ismember(QCFC_correlationMatrix_abs,unique(sort(QCFC_correlationMatrix_abs,'descend')));
    matrix=matrix+ranks; %add each ranked matrix to the other.
end
%threshold matrix to the 20% of edges with the smallest values
thresholded=unique(sort(matrix(:),'ascend'))       
[binary_nomotion_edges, ranks]=ismember(matrix,thresholded(7:ceil(length(thresholded)*0.2))); %the diagonal is 6, start at 7
%save out this matrix as a mask for each run
save(fullfile(outdir, strcat('/noMotion_edges/noMotion_edges_','gordon', run,'FIX_matrices')), 'binary_nomotion_edges')
end

%%%%%%%%%%%%%%%%
%% Run modularity analysis with only these edges for each metric
%%%%%%%%%%%%%%
%mounted cluster locally
data_dir='~/Desktop/cluster/jag/bassett-lab/hcp_Max/Data/FunctionalConnectivityMatrices'
outdir='/Users/utooley/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/arun_fc_metrics_motion/output/data/Gordon_ICA_FIX/nomotion_edges'
edges_mat='~/Desktop/cluster/jag/bassett-lab/hcp_Max/Data/noMotion_edges'

%%%%%%
%for each of four runs
%%%%%%
for i=1:4
    run=rest_runs{i};
    run_name=strcat('run',run, 'l');
    load(fullfile(edges_mat, strcat('noMotion_edges_','gordon', run,'FIX_matrices'))) %load the mask of low-motion edges for this run)
%%%%%%
%for each FC metric
%%%%%%
for j=1:6 
    metric=fc_metrics{j}
    modul.(metric)=zeros(length(subjList),1);
    avgweight.(metric)=zeros(length(subjList),1);
    num_communities.(metric)=zeros(length(subjList),1); %set up num communities

%%%%%%
%for each subject
%%%%%%
for n=1:length(subjList);
    sub=subjList(n);
    try
    file=fullfile(data_dir,strcat('gordon_',num2str(sub),run,'FIX_matrices_', metric,'.mat'));
    load(file);
    AdjMat=AdjMat.*binary_nomotion_edges; %mask AdjMat with the low-motion edges only
    %AdjMat=threshold_absolute(AdjMat,0); % uncomment for absolute value only.
    avgweight.(metric)(n,1)=mean(AdjMat(AdjMat~=0)); %get average weight for each metric
%% calculate the modularity quality index raw on each metric
if (j==4 | j == 5) %Weighted the negative connections asymmetrically, Q* as
    %recommended by Rubinov & Sporns
    [M Q]=community_louvain(AdjMat, [], [], 'negative_asym');
    modul.(metric)(n,1)=Q;
    num_communities.(metric)(n,1)=length(unique(M)); %how many communities were output
 else
    %use the default modularity
    [M Q]=community_louvain(AdjMat, []);
    modul.(metric)(n,1)=Q;
    num_communities.(metric)(n,1)=length(unique(M));
end
    catch
    disp('This subject not found')
    end
end
avgnumcommunities=mean(num_communities.(metric)(num_communities.(metric)~=0))
end
%% Save outfiles for each run
save(fullfile(outdir, strcat('modul_nomotionedges_run',int2str(i))), 'modul')
save(fullfile(outdir, strcat('numcommunities_nomotionedges_run',int2str(i))), 'num_communities')
save(fullfile(outdir, strcat('avgnumcommunities_nomotionedges_run',int2str(i))), 'avgnumcommunities')
save(fullfile(outdir, strcat('avgweight_nomotionedges_run',int2str(i))), 'avgweight')
  
end

