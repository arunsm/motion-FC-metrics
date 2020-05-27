% heuristic to compile subject list for analysis - all subjects for whom
% nogsr data is available

taskTypes = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};

for t = 1:numel(taskTypes)
    
    currentTaskType = taskTypes{t};
    
    fileName = strcat('../data/subjectsList_', currentTaskType, '.csv');
    fid = fopen(fileName, 'a');
    
    fprintf(fid, '%s\n', currentTaskType);
    
    d = dir(strcat('../data/ICAFIX_matrices_nobp/ts/gordon*', currentTaskType, '_nogsr.npy'));
    fnme = {d.name};
    
    for i = 1:numel(fnme)
        currentFileName = fnme{i};
        k = strfind(currentFileName, '_');
        currentSubjectID = str2double(currentFileName((k(1)+1):k(2)-1)); % pulling out subject ID (numeric sequence between first and second '_')
        fprintf(fid, '%d\n', currentSubjectID);
    end
    
    fclose(fid);
end
