% define step1 using local or nas folder
local = true;
check = true;
conn_flag = true;
ISC_flag = true;
o = def_nasorlocal(local,check,conn_flag,ISC_flag);

% step1_convert2nii(o.sub,o.sess,o.round, ...
%     'folder_nest',o.folder_nest, ...
%     'subpath',o.subpath, ...
%     'outsubpath',o.outsubpath, ...
%     'ftpServer',o.ftpServer);
if ~local
    step2_niiFile2local(o.sess,o.round,o.ftpServer,o.folder_nest,o.outsubpath)
end
if conn_flag
    step3_create_nproj( o.conn_steps,o.conn_slice_order,...
        'sub',o.sub, ...
        'sess',o.conn_ses, ...
        'round',o.conn_round, ...
        'condition',o.conn_conditon, ...
        'conn_proj_name',o.conn_prjName, ...
        'conn_proj_path',o.conn_prjPath, ...
        'TR',o.conn_TR, ...
        'smooth_kernel',o.conn_fwhm, ...
        'StrucVres',o.conn_strVres, ...
        'funcVres',o.conn_funcVres, ...
        'filter_band',o.conn_filtBand, ...
        'analysisName',o.conn_AnalysisName, ...
        'contrast',o.conn_contrast,...
        'mrifilepath',o.outsubpath,...
        'rawdata',o.conn_rawdata)
end

if ISC_flag
    step4_ISC_analysis(o.ISCinpath,o.ISCoutpath,o.ISCcondition, ...
        'hypoth',o.ISCcon, ...
        'sub',o.sub, ...
        'group',o.ISCgroup, ...
        'ISC_ana',false, ...
        'permutTest',true);
end

function o = def_nasorlocal(local_flag,check_par,connF,ISCF)
    % --------------------------define raw data folder direction
    % nas server parameter
    ftpServer.ip = 'ftp://xxxxxxxxxxx/';
    ftpServer.account = 'xxxxxxxxxxx';
    ftpServer.password = 'xxxxxxxxxxx';
    ftpServer.infolder = 'LabData/jeng/test/RawData'; % Nas Raw Data folder 
                             %   e.x.LabData/廣達人腦健康資料庫
    ftpServer.outfolder = 'LabData/jeng/test/niiData'; % Nas put nii file folder
                              %   e.x.LabData/jeng/test/Data
    % ------------------------------------------------------------

    % -----------------------------define local path
    % local Raw Data folder 
    subpath = '.\..\Data\rawdata\Learner';
    % local put nii file folder
    outsubpath = '.\..\Data\niidata\Learner';
    % ------------------------------------------------
    
    % ------------------------------define folder_nest
    % define folder nest sess and round are defined in variable sess and round, 
    % others is the folder name of user-defined
    % e.x. {'sess','mri','round'} --> sess_folder/'mri'/round_folder
    folder_nest = {'round'};
    % -----------------------------------------------

    % ----------------------------define sess and round
    sess = {' '};
    round = {'AO','AT','BO','BT','T1'};
    % -------------------------------------------------

    % ------------------------get subject Name
    % get age 20~30 subject
    % wordname_list = readtable("wordname_list.xlsx");
    % sub = wordname_list.ID;

    % get subject Name from Raw data folder 
    sub = {dir(outsubpath).name};
    sub = sub(contains(sub,'sub')|contains(sub,'SUB'));
    
    % get subject Name from Nas Raw data folder
    % ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
    % cd(ftpobj,ftpServer.infolder);
    % SUBname = string({dir(ftpobj).name}');
    % SUBname = SUBname(contains(SUBname,'sub'));
    % close(ftpobj);
    % sub = SUBname;

    o.taskconditon = {''};
    o.Restconditon = {'REST'};
    % ----------------------------------------

    % ------------------------------CONN project
    o.conn_ses = '';
    o.conn_round = {'AT'};
    o.conn_prjName = ['conn_',cell2mat(o.conn_round)];
    o.conn_prjPath = fullfile(pwd,'Learner_conn');
    % o.conn_steps = 'Denois_new';
    o.conn_steps =  struct("Setup",1, ...
                           "Preprocessing",1, ...
                           "Denoising",1, ...
                           "fst_Analysis",0, ...
                           "snd_Analysis",0, ...
                           "Add_Roi",0);
    % struct("Setup",0, ...
    %        "Preprocessing",0, ...
    %        "Denoising",0, ...
    %        "fst_Analysis",0, ...
    %        "snd_Analysis",0, ...
    %        "Add_Roi",1)
    if exist(fullfile(outsubpath,'MRinfo.xlsx'),'file')
        if connF
            empsub = cell(1,length(round));
            for roundi = 1:length(round)
                info = readtable(fullfile(outsubpath,'MRinfo.xlsx'),'Sheet',round{roundi}); 
                empsub{roundi} = string(info.subject)=="";
            end
            empsub = any(cell2mat(empsub),2);
            info = readtable(fullfile(outsubpath,'MRinfo.xlsx'),'Sheet',round{1}); 
            info = info(~empsub,:);
            sub = info.subject;

            o.conn_TR = unique(str2double(info.RT));
            o.conn_slice_order = str2num(unique(string(info.sliceorder))); % use str2num instead of str2double, using str2double will get nan
        else
            o.conn_TR = 0;
            o.conn_slice_order = '';
        end
        if size(o.conn_slice_order,1) > 1, o.conn_slice_order = o.conn_slice_order'; end
        if size(o.conn_TR,1) > 1
            error('infoData "TR" has different between subjects');
        end
        if size(o.conn_slice_order,1) > 1
            warning('infoData "sliceorder" has different between subject');
            o.conn_slice_order = '';
        end
    else
        o.conn_TR = 2.4;
        o.conn_slice_order = 'interleaved (bottom-up)';
    end
    
    o.conn_fwhm = 8;
    o.conn_strVres = 1;
    o.conn_funcVres = 2;
    o.conn_filtBand = [0.008, 0.09];
    o.conn_AnalysisName = 'SBC_01';
    o.conn_contrast = {1};
    o.conn_ROIpath = [];
    o.conn_conditon = {};
    o.conn_rawdata = true;

    % ------------------------------------------

    % --------------------------------------ISC
    o.ISCinpath = outsubpath;
    o.ISCoutpath = '.\..\Data\ISCdata\';
    o.ISCgroup = {'Interpreter','Learner'};
    o.ISCcondition = {'AO','AT','BO','BT'};
    o.ISCcon = {[1,-1]};

    % o.ISCcondition = {'AO','AT','BO','BT'};
    % o.ISCcon = {[1,0,1,0]};
    % o.ISCcon = {[1,0,1,0,-1,0,-1,0]};
    % -----------------------------------------

    if local_flag
        o.subpath = subpath;
        o.ftpServer = struct();
    else
        o.ftpServer = ftpServer;
        o.subpath = [];
    end
    o.outsubpath = outsubpath;
    o.folder_nest = folder_nest;
    o.sess = sess;
    o.round = round;
    o.sub = sub;
    if size(o.sess,2) ~= 1, o.sess = o.sess';end
    if size(o.round,2) ~= 1, o.round = o.round';end
    if size(o.sub,2) ~= 1, o.sub = o.sub';end
    o.done = false;
    if check_par
        o = checkGUI([],local_flag,o,connF,ISCF);
    end
end



