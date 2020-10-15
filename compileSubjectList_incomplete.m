% script to compile a list of subjects for whom processing is not yet complete

taskTypes = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
subjectsList_incomplete_allTasks = [];

for t = 1:numel(taskTypes)
    currentTaskType = taskTypes{t}
    s = importdata(strcat('../data/subjectsList_', currentTaskType, '.csv'));
    subjectsList = s.data;
    
    % gordon
    d = dir(strcat('../data/FunctionalConnectivityMatrices_gsr_nofilter/gordon_*', currentTaskType, '_ts_TikhonovPartialCorrelation.mat'));
    fnme = {d.name};
    subjectsList_complete = [];
    for i = 1:numel(fnme)
        currentFileName = fnme{i};
        k = strfind(currentFileName, '_');
        currentSubjectID = str2double(currentFileName((k(1)+1):k(2)-1)); % pulling out subject ID (numeric sequence between first and second '_')
        subjectsList_complete = [subjectsList_complete; currentSubjectID];
    end
    
    subjectsList_incomplete = setdiff(subjectsList, subjectsList_complete)
    subjectsList_incomplete_allTasks = [subjectsList_incomplete_allTasks; subjectsList_incomplete];
    
    % yeo
    d = dir(strcat('../data/FunctionalConnectivityMatrices_gsr_nofilter/yeo_*', currentTaskType, '_ts_TikhonovPartialCorrelation.mat'));
    fnme = {d.name};
    subjectsList_complete = [];
    for i = 1:numel(fnme)
        currentFileName = fnme{i};
        k = strfind(currentFileName, '_');
        currentSubjectID = str2double(currentFileName((k(2)+1):k(3)-1)); % pulling out subject ID (numeric sequence between second and third '_')
        subjectsList_complete = [subjectsList_complete; currentSubjectID];
    end
    
    subjectsList_incomplete = setdiff(subjectsList, subjectsList_complete)
    subjectsList_incomplete_allTasks = [subjectsList_incomplete_allTasks; subjectsList_incomplete];
end

subjectsList_incomplete_allTasks = unique(subjectsList_incomplete_allTasks);

fid = fopen('subjectsList_incomplete.csv', 'w');
fprintf(fid, '%d\n', subjectsList_incomplete_allTasks);
fclose(fid);