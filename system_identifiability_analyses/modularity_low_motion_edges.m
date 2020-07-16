%% THIS SCRIPT RUNS MODULARITY ANALYSIS ON 4 RUNS of EACH SUBJECT, LOW-MOTION EDGES ONLY, USING YEO100 or GORDON PARCELLATION
% ICA-FIX MATRICES ONLY

%%%%%%%%%%%%%%%%
%% Run modularity analysis with only these edges for each metric
%%%%%%%%%%%%%%
%mounted cluster locally
subj_dir='/cbica/home/mahadeva/motion-FC-metrics/data/Covariates'
data_dir='/cbica/home/mahadeva/motion-FC-metrics/data/FunctionalConnectivityMatrices_gsr_filter'
outdir='//cbica/home/tooleyu/arun_fc_metrics_motion/output/Schaefer_100_ICA_FIX/nomotion_edges'
edges_mat='/cbica/home/mahadeva/motion-FC-metrics/data/noMotion_edges'
addpath(genpath('/cbica/home/tooleyu/motion-FC-metrics/system_identifiability_analyses/code/functions'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/'))
%subject list
subjList=readtable(fullfile(subj_dir, 'S1200_Release_Subjects_Demographics.csv'));
subjList=subjList.Subject;

rest_runs={'_REST1_LR_', '_REST1_RL_','_REST2_LR_','_REST2_RL_'};
fc_metrics={'Coherence', 'MutualInformation', 'MutualInformationTime','Pearson','PartialCorrelation', 'WaveletCoherence'};

%%%%%%
%for each of four runs
%%%%%%
for i=1:4
    run=rest_runs{i};
    run_name=strcat('run',run, 'l');
%%%%%%
%for each FC metric
%%%%%%
for j=1:6 
    metric=fc_metrics{j}
    load(fullfile(edges_mat, strcat('lowMotionEdges_yeo_100_', metric,'.mat'))) %load the mask of low-motion edges for this run)
    mask=squareform(idx_lowMotionEdges_currentFCmethod_allScans);
    modul.(metric)=zeros(length(subjList),1);
    avgweight.(metric)=zeros(length(subjList),1);
    num_communities.(metric)=zeros(length(subjList),1); %set up num communities
%%%%%%
%for each subject
%%%%%%
for n=1:length(subjList);
    sub=subjList(n);
    try
    file=fullfile(data_dir,strcat('yeo_100_',num2str(sub),run,'FIX_matrices_', metric,'.mat'));
    load(file);
    AdjMat=AdjMat.*mask; %mask AdjMat with the low-motion edges only
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
try
    outfile=dataset(avgweight.Pearson, avgweight.PartialCorrelation, avgweight.Coherence, avgweight.WaveletCoherence, avgweight.MutualInformation, avgweight.MutualInformationTime, modul.Pearson, modul.PartialCorrelation, modul.Coherence, modul.WaveletCoherence, modul.MutualInformation, modul.MutualInformationTime)
    header={'avgweight_Pearson',	'avgweight_PartialCorrelation',	'avgweight_Coherence',	'avgweight_WaveletCoherence',	'avgweight_MutualInformation',	'avgweight_MutualInformationTime',	'modul_Pearson',	'modul_PartialCorrelation',	'modul_Coherence',	'modul_WaveletCoherence',	'modul_MutualInformation',	'modul_MutualInformationTime'}
    outfile.Properties.VarNames=header
    filename=strcat('modularity_nomotionedges_yeo_',run, '071420.csv') %need to figure out how to get the headers to work.
    export(outfile,'File',fullfile(outdir,filename),'Delimiter',',')
catch
    
    disp('saving csvs didnt work')
end

end

