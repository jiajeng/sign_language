function step1_convert2nii(sub,sess,round,varargin)
% load data from nas(/rawdata) and store it to nas(/process/file_prep)
% folder nest -- 
% subject
%   session
%       round
%           task
%           REST
%           T1
%           ...
in = finputcheck(varargin, ...
{'folder_nest'  'cell'      []      [];
 'subpath'      'string'    []      [];
 'outsubpath'   'string'    []      [];
 'ftpServer'    'struct'    []      [];
 });

% ------------define use nas get raw data or local folder

% -------------------------------------------------------
ftpServer = in.ftpServer;

% --------------------------define log variable
if ~exist('log.mat','file')
    log = cell2table(cell(length(unique(sub)),length(round))); 
    log.Properties.VariableNames = strcat(round,'folder');
    RowNames = 1:length(unique(sub));
    log.Properties.RowNames = convertStringsToChars(string(num2str(RowNames')));
else
    load('log.mat');
end

infoVar = {'subject','Fpath','voxelSize','ImageSize','sliceorder','RT'};
infotab = cell(length(sub),length(infoVar));
infotab = cell2table(infotab);
infotab.Properties.VariableNames = infoVar;
info = struct();
for nround = 1:length(round)
    info.(round{nround}) = infotab;
end


% --------------------------------------------

Ifold = in.folder_nest;
Ifold(cellfun(@(x) x=="sess",in.folder_nest)) = {'char(sess(nsess))'};
Ifold(cellfun(@(x) x=="round",in.folder_nest)) = {'round{nround}'};
idx = find(~(cellfun(@(x) x=="sess",in.folder_nest) | cellfun(@(x) x=="round",in.folder_nest)));
for i = idx
    Ifold(i) = {['char("',Ifold{i},'")']};
end
Ifold = ['Ifolder = fullfile(',strjoin(Ifold,','),');'];


% ------------------------------------main

for nsub = 1:length(sub)
    for nsess = 1:length(sess)
        for nround = 1:length(round)
            try
            eval(Ifold);
            Ofolder = fullfile('');
            if round{nround} == "T1", fT1 = true; else, fT1 = false; end
            [log,info.(round{nround})] = DICOM2nii(char(sub(nsub)),Ifolder,Ofolder,log,info.(round{nround}), ...
                'ftpServer',ftpServer, ...
                'subpath',in.subpath, ...
                'outsubpath',in.outsubpath, ...
                'T1',fT1, ...
                'subNum',nsub);
            catch ME
            end
        end
    end
end

save('niifldr_log.mat',"log");
for nround= 1:length(round)
    writetable(info.(round{nround}),fullfile(in.outsubpath,"MRinfo.xlsx"),"Sheet",round{nround});
end
% -------------------------------------------------------
end

%% define function
function [log,info] = DICOM2nii(sub,input_folder,output_folder,log,info,varargin)
    in = finputcheck(varargin, ...
    {'sub'              'string'    []              sub; ...
     'subNum'           'real'      []              [];
     'input_folder'     'string'    []              input_folder; ...
     'output_folder'    'string'    []              output_folder; ...
     'ftpServer'        'struct'    []              struct();
     'subpath'          'string'    []              [];
     'bufferfolder'     'string'    []              'buffer';
     'outsubpath'       'string'    []              [];
     'T1'               'boolean'   []              [];
     'samefileM'        'string'    []              'filetailMark';
     });

    ftpServer = in.ftpServer;
    if isempty(fieldnames(ftpServer))
        in.bufferfolder = in.subpath;
        dpath = in.bufferfolder;
        opath = fullfile(in.outsubpath);
    else
        dpath = fullfile(in.bufferfolder,ftpServer.infolder);
        opath = fullfile(in.bufferfolder,ftpServer.outfolder);
    end
 
    rename = in.sub;
    % rename = split(in.sub,'-');
    % rename = ['SUB',rename{end}];

    if ~isempty(fieldnames(ftpServer))
        ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
        cd(ftpobj,ftpServer.infolder);
        SUBname = string({dir(ftpobj).name}');
        SUBname = SUBname(contains(SUBname,in.sub));
        close(ftpobj);
        % check if has same subject name 
        if length(unique(SUBname)) > 1
            flag = ones(1,length(SUBname));
            % check subject name does have wanted folder
            for nsub = 1:length(SUBname)
                ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
                try
                    in.input_folder(in.input_folder=='\') = '/';
                    cd(ftpobj,[ftpServer.infolder,'/',char(SUBname(nsub)),'/',in.input_folder]);
                catch
                    flag(nsub) = 0;
                end
                close(ftpobj);
            end
            % if all subject name has same folder, then get the latest one(_n)
           if sum(flag) > 1
                tmp = split(SUBname,'_');
                switch in.samefileM
                    case "Date"
                        tmp = str2double(tmp(:,1));
                    case "filetailMark"
                        tmp = str2double(tmp(:,end));
                end
                SUBname = SUBname(tmp == max(tmp));
            else    
                SUBname = SUBname(logical(flag));
            end
        end
    else
        wd = pwd;
        ftpobj = in.subpath;
        SUBname = string({dir(ftpobj).name}');
        SUBname = SUBname(contains(SUBname,in.sub));
        % check if has same subject name 
        if length(unique(SUBname)) > 1
            flag = ones(1,length(SUBname));
            % check subject name does have wanted folder
            for nsub = 1:length(SUBname)
                try
                    cd(fullfile(ftpobj,char(SUBname(nsub)),in.input_folder));
                catch
                    flag(nsub) = 0;
                end
                cd(wd);
            end
            % if all subject name has same folder, then get the latest one
            if sum(flag) > 1
                tmp = split(SUBname,'_');
                switch in.samefileM
                    case "Date"
                        tmp = str2double(tmp(:,1));
                    case "filetailMark"
                        tmp = str2double(tmp(:,end));
                end
                SUBname = SUBname(tmp == max(tmp));
            else    
                SUBname = SUBname(logical(flag));
            end
        end
        cd(wd);
    end
    
    in.sub = char(SUBname);
    datapath = fullfile(dpath,in.sub,in.input_folder);
    if ~exist(datapath,"dir")
        datapath = ftpgetfile(ftpServer,in,in.input_folder);
    end
    dataname = {dir(datapath).name}';
    dataname = string(dataname(3:end));
    dataname = dataname(contains(dataname,'.IMA'));
    dataname = convertStringsToChars(dataname);

    %%%%%%%%%%%%%%%%%%%%%%%%%%% main code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    folderN = split(in.input_folder,filesep);
    folderN = folderN{end};
    log{in.subNum,[folderN,'folder']} = {datapath};

    in.sub = rename;

    % output file 
    outputpath = fullfile(opath,in.sub,in.output_folder,in.input_folder);
    if ~exist(outputpath,'dir')
        mkdir(outputpath);
    end
    
    matlabbatch{1}.spm.util.import.dicom.data = fullfile(datapath,dataname);
    %%
    matlabbatch{1}.spm.util.import.dicom.root = 'flat';
    matlabbatch{1}.spm.util.import.dicom.outdir = {outputpath};
    matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
    matlabbatch{2}.spm.util.cat.vols(1) = cfg_dep('DICOM Import: Converted Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{2}.spm.util.cat.name = [rename,'_4D.nii'];
    matlabbatch{2}.spm.util.cat.dtype = 4;
    spm('defaults', 'FMRI');
    spm_jobman('initcfg')
    spm_jobman('run',matlabbatch);
    
    if in.T1
        wd = pwd;
        cd(outputpath)
        gzip('*.nii');
        cd(wd);
    end
    tmp = niftiinfo(fullfile(outputpath,[rename,'_4D.nii']));
    info(in.subNum,"subject") = {rename};
    info(in.subNum,"Fpath") = {erase(datapath,in.bufferfolder)};
    info(in.subNum,"ImageSize") = {num2str(tmp.ImageSize)};
    info(in.subNum,"voxelSize") = {num2str(tmp.PixelDimensions)};
    tmp = dicominfo(fullfile(datapath,dataname{1}));
    try
        sliceorder = tmp.Private_0019_1029;
        [~,sliceorder] = sort(sliceorder);
        if size(sliceorder,1) > 1,sliceorder = sliceorder';end
    catch
        sliceorder = nan;
    end
    info(in.subNum,"sliceorder") = {num2str(sliceorder)};
    info(in.subNum,"RT") = {num2str(tmp.RepetitionTime/1000)};
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    in.sub = rename;
    % save data to nas
    if ~isempty(fieldnames(ftpServer))
        ftpputfile(ftpServer,in)
        if exist(in.bufferfolder,"dir")
            rmdir(in.bufferfolder,'s')
        end
    end
end


function path = ftpgetfile(ftpServer,in,folder)
    if ~isempty(fieldnames(ftpServer))
        ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
        targetfolder = fullfile(ftpServer.infolder,in.sub,folder);
        targetfolder(targetfolder=='\') = '/';
        mget(ftpobj,targetfolder,in.bufferfolder);
        close(ftpobj)
        path = char(fullfile(in.bufferfolder,ftpServer.infolder,in.sub,folder));
    else
        path = char(fullfile(in.bufferfolder,in.sub,folder));
    end
end


function ftpputfile(ftpServer,in)
    ftpObj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
    folder = split(ftpServer.outfolder,'/');
    folder = fullfile(folder{1:3});
    folder(folder=='\') = '/';
    localfolder = fullfile(in.bufferfolder,ftpServer.outfolder);
    try
        cd(ftpObj,folder);
    catch
        mkdir(ftpObj,folder);
        cd(ftpObj,folder);
    end
    mput(ftpObj,localfolder)
    close(ftpObj);
end
