folderPath = "../xcpEngine-QCFC/FC_matrices/AROMA+GSR/Gordon/";
%folderPath = "../motion/subjects";
d = dir(strcat(folderPath, filesep, "*.txt"));
fnme = {d.name};
nSubjects = numel(fnme);

fid = fopen("../xcpEngine-QCFC/filePaths.csv", 'a');

for i = 1:nSubjects
    currentFilePath = strcat("/home/amahad/xcpEngine-QCFC/FC_matrices/AROMA+GSR/Gordon/", fnme{i});
    fprintf(fid, '%s\n', currentFilePath);
end

fclose(fid);