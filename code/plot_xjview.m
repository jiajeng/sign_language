basepath = 'C:\Users\user\Desktop\jeng\sign_language\sign_langa_TW\Data\ISCdata';
folderN = {dir(basepath).name};
folderN = folderN(contains(folderN,'stats'));
opth = 'C:\Users\user\Desktop\jeng\sign_language\sign_langa_TW\Data\ISCdata\result';
if ~exist(opth,'dir'), mkdir(opth);  end


for i = 1:length(folderN)
    try
        filepath = fullfile(basepath,folderN{i},'thres');
        niifile = {dir(filepath).name};
        niifile = niifile(contains(niifile,'.nii'));
        for j = 1:length(niifile)
            xjview(fullfile(filepath,niifile{j}));
            filename = split(niifile{j},'.');
            filename = filename{1};
            filename = split(niifile{j},'_');

            filename = [folderN{i},'_',filename{end-1},'_',filename{end}];
            saveas(gcf,fullfile(opth,[filename,'.fig']))
            close all;
            close(findall(groot,'type','figure','tag','Msgbox_Warning Dialog'))
        end
    catch ME
        disp(folderN{i})
    end 
end
