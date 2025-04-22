subpath = 'C:\Users\user\Desktop\jeng\sign_language\sign_langa_TW\Data\rawdata\Learner';
round = {'AT','BT','AO','BO','LOCALIZER','REST'};
subject = {dir(subpath).name};
subject = subject(~cellfun(@(x) contains(x,'.'),subject));
wkdir = pwd;

log = ones(length(subject),length(round));

for nsub = 1:length(subject)
    foldName = {dir(fullfile(subpath,subject{nsub})).name};
    for nroud = 1:length(round)
        try
            srcfold = foldName{contains(foldName,round{nroud})};
            destfold = round{nroud};
            subfolder = fullfile(subpath,subject{nsub});
            if string(srcfold)~=string(destfold)
                movefile(fullfile(subfolder,srcfold),fullfile(subfolder,destfold));
            end
        catch
            log(nsub,nroud) = 0;
        end
    end
end
log = array2table(log);
log.Properties.VariableNames = round;
log.Properties.RowNames = subject;
save("Learner_Log.mat","log")