subpath = 'C:\Users\user\Desktop\jeng\sign_language\sign_langa_TW\Data\rawdata\Interpreter';
subject = {dir(subpath).name};
subject = subject(~cellfun(@(x) contains(x,'.'),subject));
nsubject = cellfun(@(x) del(x),subject,'UniformOutput',false);
for i = 1:length(subject)
    srcname = subject{i};
    detname = nsubject{i};
    if string(srcname)~=string(detname)
        movefile(fullfile(subpath,srcname),fullfile(subpath,detname));
    end
end
function x = del(x)
    x(x=='P') = [];
    x = cat(2,'SUB',x);
end