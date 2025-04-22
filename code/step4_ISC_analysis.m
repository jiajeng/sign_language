function step4_ISC_analysis(datafolder,resultfolder,condition_label,varargin)
    in = finputcheck(varargin, {
        'sub'       'cell'      []  {};
        'group'     'cell'      []  {};
        'maskdir'   'string'    []  '.\mask\mask.volume.brainmask.nii';
        'fileP'     'string'    []  'dswau*.nii';
        'permTimes' 'real'      []  2000;
        'Pthres'    'real'      []  0.05;
        'rmVst'     'real'      []  1;
        'rmVend'    'real'      []  1;                  
        'hypoth'    'cell'      []  {1};
        'ISC_ana'   'boolean'   []  false;
        'nest_folder' 'string'  []  '';
        'permutTest' 'boolean'  []  true;
         });
    Dens = false;
    % remove the first & last few volumes (unit = TR)
    remove_volume_start = in.rmVst;
    remove_volume_end = in.rmVend;
    
    % import mask
    V_mask = spm_vol(in.maskdir); 
    mask = spm_read_vols(V_mask);
    Brainout = single(find(mask==0));

    % filename pattern
    fileP = in.fileP;
    group = in.group;
   
    if in.ISC_ana
        %% PART 1: calculate inter-subject correlation for each voxel
        disp('ISC analysis')
        % condiG_label = [];
        for cond_idx=1:length(condition_label) % for every condition
            % --------------------------------------define stc variable
            for group_idx = 1:length(group)
                subject = {dir(fullfile(datafolder,group{group_idx})).name};
                subject = subject(contains(subject,'SUB'));
                niifile = dir(fullfile(datafolder,group{group_idx},subject{1},condition_label{cond_idx},fileP));
                tmp = spm_vol(fullfile(niifile.folder,niifile.name));
                tmp = spm_read_vols(tmp);
                vX = size(tmp,1); vY = size(tmp,2); vZ = size(tmp,3); % get image voxel x*y*z (vX--x,vY--y,vZ--z)
                stc = nan(vX*vY*vZ,size(tmp,4),length(subject)); % In certain condition all subject .nii data -- reshape to volumn*slice(time point)*subject
    
                for subj_idx=1:length(subject)
                    niifile = dir(fullfile(datafolder,group{group_idx},subject{subj_idx},condition_label{cond_idx},fileP));
        
                    V_data = spm_vol(fullfile(niifile.folder,niifile.name));  % read nii file attribute(V_data)
                    data = spm_read_vols(V_data); % get 4-D data(data)
                    disp(['loading ',fullfile(niifile.folder,niifile.name)])
            
                    % 1 subject = 4D img to 2D => put every subject in a matrix   ex.img 2*2*2    *10  *20 --> stc 8  *10  *20
                    %                                                                    n*m*slice*time*subject    img*time*subject
                    for t = 1:size(data,4) 
                        thisV = data(:,:,:,t); % find every time point image
                        thisV(Brainout) = NaN; % set vales outside the mask to NaN
                        stc(:,t,subj_idx) = reshape(thisV,vX*vY*vZ,1); % set a image data(3D) to a column(1D) ex. img 2*2*2 --> stc 8*1 (stc--voxel*slice*subject)
                    end
                end
                clear("mask","V_mask","tmp","thisV","data")

                % define output folder
                outfolder = fullfile(resultfolder,group{group_idx});

                % ii = find(isnan(stc)); % find NaN
                stc(isnan(stc))=randn(size(find(isnan(stc)))).*eps; % set NaN to a extreme small value (random number(base on normal distribute)*2^-16)
            
                %remove the first and last few time points
                stc=stc(:,remove_volume_start:end-remove_volume_end,:);
                % -----------------------------------------------------------------
                
                
                if Dens
                    confound_polynomial_order=2;
                    confound_sinusoidal_order=3;
                    % ----------------------------------------------remove global mean
                    for s_idx=1:size(stc,3) % for every subject
                        tmp=squeeze(stc(:,:,s_idx)); % find every subject image (allvoxel*timepoint)
                        ga=mean(tmp,1); % find every img mean(ga)
                
                        tmp=tmp-tmp*ga'*inv(ga*ga')*ga; %?? remove mean
                
                        stc(:,:,s_idx)=tmp; % set new img in same space
                    end
                    % ---------------------------------------------------------------
               
                    % -----------------------------------------------------denoising??
                    timeVec=(1:size(stc,2))'; % time vector(1:number of slice)
                    D_poly=[];
                    D_sinu=[];
                    D_poly=ones(length(timeVec),1); 
                    for i=1:confound_polynomial_order
                        tmp=timeVec.^(i);
                        D_poly(:,i+1)=fmri_scale(tmp(:),1,0); % scale data to 0:1
                    end
                    for i=1:confound_sinusoidal_order
                        D_sinu(:,i*2-1)=sin(timeVec.*i./timeVec(end).*pi);
                        D_sinu(:,i*2)=cos(timeVec.*i./timeVec(end).*pi); % [sin cos sin cos]
                    end
                    D=cat(2,D_poly,D_sinu);
                
                    if(~isempty(D))
                        D_prep=D*inv(D'*D)*D'; % D matrix sub stc(time*subject)
                    else
                        D_prep=[];
                    end
                    % ----------------------------------------------------------------
                else
                    D_prep = [];
                end
                fprintf('\n');
            
                % -----------------------------------------every voxel correlation 
                Zm = zeros(1,size(stc,1));
                Zmd = Zm;
                % save all correlation coef.
                if cond_idx==1
                    % corr_buffer --> condition*voxel*subject*subject
                    % corr_buffer=zeros(length(condition_label)*length(group),size(stc,1),size(stc,3),size(stc,3)); % set a buffer to restore the correlation value
                    Z = zeros(length(condition_label)*length(group),size(stc,1),nchoosek(length(subject),2)); %  cc--(set z in here) (condition,voxel,C(subject,2))
                    corr_tril = zeros(size(stc,1),nchoosek(length(subject),2));
                end
                disp('start parrallel')
                vol_idx = single(1:size(stc,1));
                vol_idx = vol_idx(~ismember(vol_idx,Brainout));
                t = 0;
                for v_idx=vol_idx % find every voxel 
                    t = t+1;
                    if(mod(t,100)==0)
                        fprintf('[%1.1f%%]...\r',v_idx/size(stc,1)*100); % print process rate 
                    end
                    % data --> slice*subject
                    if(~isempty(D_prep))
                        data=squeeze(stc(v_idx,:,:))-D_prep*squeeze(stc(v_idx,:,:)); % data(voxel value) -- slice*subject(for every voxel)
                    else
                        data=squeeze(stc(v_idx,:,:));
                    end
                    if any(isnan(data)), continue; end
                    data(abs(data)<eps)=randn(size(find(abs(data)<eps))).*eps; % if voxel value is lower than 2^-52 set it a new value


                    % corr_buffer(cond_idx,v_idx,:,:)=corrcoef(data); % find evrey voxel in subject correlation
                    %
                    % condition*voxel*subject*subject
                    %                 (corrlation matrix)
                    
                    
                    tmp=tril(corrcoef(data),-1); % lower trangular part of matrix
                    
                    rr=tmp(tmp~=0); % tmp -- subject*subject(upper trangular part is zero) set tmp to (C(subject,2)*1) ---2D to 1D
                    corr_tril(v_idx,:) = rr;
                    if any(rr==1)
                        row = ceil(find(tmp == 1)/length(subject));
                        col = mod(find(tmp == 1),length(subject));
                        error([subject{row},' and ',subject{col},' R value equal to 1.']);
                    end

                    z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3)); % compute z score
                    zm=mean(z);
                    zmd=median(z);

                    if size(z,1) > 1,z = z';end
                    if(find(isnan(z))), z = 0; end % if z have nan give control to keyboard
                    Z(cond_idx,v_idx,:) = z; % set z


                    Zm(v_idx)=zm; % store every mean in z
                    Zmd(v_idx)=zmd; % store every median in z                
                end
                % ----------------------------------------------------------
        
                fprintf('\n');
                if ~exist(fullfile(outfolder,'pairs'),'dir')
                    mkdir(fullfile(outfolder,'pairs'))
                end
            
                if ~exist(fullfile(outfolder,'group'),'dir')
                    mkdir(fullfile(outfolder,'group'))
                end
            
                if ~exist(fullfile(outfolder,'pairs_Z',condition_label{cond_idx}),'dir')
                    mkdir(fullfile(outfolder,'pairs_Z',condition_label{cond_idx}))
                end

                if ~exist(fullfile(outfolder,'paris_corr',condition_label{cond_idx}),'dir')
                    mkdir(fullfile(outfolder,'pairs_corr',condition_label{cond_idx}))
                end

                % save pair  index
                pair_i = find(tril(ones(size(Z,3)),-1));
                pair_pos = zeros(length(pair_i),3);
                pair_pos(:,1) = (1:length(pair_i))';
                for i = 1:length(pair_i)
                    idx = pair_i(i);
                    pair_pos(i,2) = ceil(idx/length(subject)); % column
                    tmp = mod(idx,length(subject)); % row
                    if ~tmp, tmp = length(subject); end
                    pair_pos(i,3) = tmp;
                end
                pair_pos = num2cell(pair_pos);
                pair_pos = cellfun(@num2str,pair_pos,'UniformOutput',false);
                pair_pos(:,1) = cellfun(@(x) ['pair',x],pair_pos(:,1),'UniformOutput',false);
                % writecell(pair_pos,fullfile(outfolder,'pairs','pairName'));
        
                for pair_idx = 1:size(Z,3) % save every pair 
                    % Z
                    V2save = reshape(Z(cond_idx,:,pair_idx),vX,vY,vZ); % V2save -- save all voxel z value in a .nii (all file number will be C(subject,2))
                    V_header = V_data(1);
                    V_header.fname = fullfile(outfolder,'pairs_Z',condition_label{cond_idx},['ISC_',condition_label{cond_idx},'_pairZ_',sprintf('%03d',pair_idx),'.nii']);
                    spm_write_vol(V_header,V2save); 

                    % correlation
                    V2save = reshape(corr_tril(:,pair_idx),vX,vY,vZ);
                    V_header.fname = fullfile(outfolder,'pairs_corr',condition_label{cond_idx},['ISC_',condition_label{cond_idx},'_pairC_',sprintf('%03d',pair_idx),'.nii']);
                    spm_write_vol(V_header,V2save);
                end
                
                % save mean and median Z
                V_header_Zm =  V_data(1);  V_header_Zm.fname = fullfile(outfolder,'group',['ISC_',condition_label{cond_idx},'_mean.nii']);
                V_header_Zm.dt(1)=16;
                V_header_Zmd =  V_data(1); V_header_Zmd.fname = fullfile(outfolder,'group',['ISC_',condition_label{cond_idx},'_median.nii']);
                V_header_Zmd.dt(1)=16;
                
                spm_write_vol(V_header_Zm,reshape(Zm, vX, vY, vZ));% save mean data
                spm_write_vol(V_header_Zmd,reshape(Zmd, vX, vY, vZ)); % save median data
                % condiG_label = cat(1,condiG_label,{[condition_label{cond_idx},'_',group{group_idx}]});

                clear("corr_tril","Zm","Zmd","V2save","V_header","V_header_Zm","V_header_Zmd")

            end
        end
        stcSize = size(stc);
        nsub = length(subject);
        save(fullfile(resultfolder,'var.mat'),"stcSize","nsub");
        cv = string(who);
        cv = convertStringsToChars(cv(cv~="Z"& ...
                                      cv~="stcSize"& ...
                                      cv~="corr_buffer"& ...
                                      cv~="in"& ...
                                      cv~="condiG_label"& ...
                                      cv~="condition_label"& ...
                                      cv~="group"& ...
                                      cv~="subject"& ...
                                      cv~="datafolder"& ...
                                      cv~="resultfolder"& ...
                                      cv~="V_data"& ...
                                      cv~="vX"& ...
                                      cv~="vY"& ...
                                      cv~="vZ"));
        clear(cv{:});
    end
    %% PART 2: calculate the analysis effect and run a permutation test
    
    if in.permutTest
        
        % --------------------------------------------define variable
   
        %specify the number of permutation and stats threshold
        nPerm = in.permTimes;
        p_threshold = in.Pthres;
    
        % define analysis of interest
        hypothesis = in.hypoth;
        
        if ~in.ISC_ana
            load(fullfile(resultfolder,'var.mat'));
            condiG_label = [];
            for group_idx = 1:length(group)
                for cond_idx = 1:length(condition_label)
                     condiG_label = cat(1,condiG_label,{[condition_label{cond_idx},'_',group{group_idx}]});
                end
            end
            
            % get Z
            fprintf('Get Z variable ...\n')
            
            Z = zeros(length(condition_label)*length(group),stcSize(1),nchoosek(nsub,2)); %  cc--(set z in here) (condition,voxel,C(subject,2))
            
            for cond_idx=1:length(condition_label)
                for group_idx = 1:length(group)
                    niifile = {dir(fullfile(resultfolder,group{group_idx},'pairs_Z',condition_label{cond_idx})).name};
                    niifile = niifile(3:end);

                    for pair_idx = 1:length(niifile)
                        disp(['loading ',fullfile(resultfolder,group{group_idx},'pairs_Z',condition_label{cond_idx},niifile{pair_idx}),'...'])
                        V_data = spm_vol(fullfile(resultfolder,group{group_idx},'pairs_Z',condition_label{cond_idx},niifile{pair_idx}));
                        tmp = spm_read_vols(V_data);
                        Z(group_idx+(cond_idx-1)*length(group),:,pair_idx) = reshape(tmp,1,stcSize(1));
                    end
                end
            end
            vX = size(tmp,1); vY = size(tmp,2); vZ = size(tmp,3); % get image voxel x*y*z (vX--x,vY--y,vZ--z)
            % get corr
            fprintf('Get corr variable ...\n')
            corr_buffer = single(zeros(length(condition_label)*length(group),stcSize(1),stcSize(3),stcSize(3))); % set a buffer to restore the correlation value
            for cond_idx=1:length(condition_label)
                for group_idx = 1:length(group)
                    niifile = {dir(fullfile(resultfolder,group{group_idx},'pairs_corr',condition_label{cond_idx})).name};
                    niifile = niifile(3:end);

                    corr_tril = zeros(stcSize(1),length(niifile));
                    for pair_idx = 1:length(niifile)
                        disp(['loading ',fullfile(resultfolder,group{group_idx},'pairs_corr',condition_label{cond_idx},niifile{pair_idx}),'...'])
                        V_data = spm_vol(fullfile(resultfolder,group{group_idx},'pairs_corr',condition_label{cond_idx},niifile{pair_idx}));
                        tmp = spm_read_vols(V_data);
                        corr_tril(:,pair_idx) = reshape(tmp,vX*vY*vZ,1);
                    end
                    for vol_idx = 1:stcSize(1)
                        corr_buffer(group_idx+(cond_idx-1)*length(group),vol_idx,:,:) = get_corr_buffer(nsub,corr_tril(vol_idx,:));
                    end
                end
            end
            

        end
        clear tmp corr_tril;
        % -----------------------------------------main statistical
        disp('statistical')
        effect_mean = single(zeros(length(hypothesis),stcSize(1)));
        effect_median = single(effect_mean);
        fname = cell(1,length(hypothesis));
        for hypothesis_idx=1:length(hypothesis)  % set every hypothesis mean and median and save to .nii
            tmp = ['stats_'];
            for i = 1:length(condiG_label)
                hcon = hypothesis{hypothesis_idx};
                tmp = cat(2,tmp,[condiG_label{i},'(',num2str(hcon(i)),')']);
            end
            fname{hypothesis_idx} = tmp;
            
            if ~exist(fullfile(resultfolder,fname{hypothesis_idx}),'dir')
                mkdir(fullfile(resultfolder,fname{hypothesis_idx}))
            end
            
            
            for v_idx=1:stcSize(1)
                tmp=hypothesis{hypothesis_idx}(:)'*squeeze(Z(:,v_idx,:)); %  cc--(set z in here) (condition,voxel,C(subject,2))
                effect_mean(hypothesis_idx,v_idx)=mean(tmp); % set hypothesis mean and median
                effect_median(hypothesis_idx,v_idx)=median(tmp);
            end
            
            % save meanz and medianz value
            V_header_mean =  V_data(1);  V_header_mean.fname = fullfile(resultfolder,fname{hypothesis_idx},['ISC_hypothesis',sprintf('%02d',hypothesis_idx),'_mean_z.nii']);
            V_header_mean.dt(1)=16;
            V_header_median =  V_data(1); V_header_median.fname = fullfile(resultfolder,fname{hypothesis_idx},['ISC_hypothesis',sprintf('%02d',hypothesis_idx),'_median_z.nii']);
            V_header_median.dt(1)=16;
            
            spm_write_vol(V_header_mean,reshape(effect_mean(hypothesis_idx,:), vX, vY, vZ));
            spm_write_vol(V_header_median,reshape(effect_median(hypothesis_idx,:), vX, vY, vZ));
        end
        clear("Z","tmp","V_header_median","V_header_mean","mask","ans");
       
        for hypothesis_idx=1:length(hypothesis) % ---------hypothesis
            
            outputname = fullfile(resultfolder,fname{hypothesis_idx}, ['ISC_hypothesis',sprintf('%02d',hypothesis_idx)]);
    
            zm_perm_p= single(nan(1,stcSize(1)));
            zmd_perm_p=zm_perm_p;
    
            tril_idx=single(find(tril(ones(stcSize(3)),-1))); % find lower triangular part index
            % cc_perm = zeros(length(condition_label),nchoosek(subject,2)); %cc_perm -- save all condition data(condition*C(subject,2))
            % ccm_perm = zeros(length(condition_label),1);
            % ccmd_perm = ccm_perm;
            null_mean = single(zeros(stcSize(1),nPerm));
            null_median = single(zeros(stcSize(1),nPerm));
            hypo = hypothesis{hypothesis_idx};
            
            par_com = 1;
            if par_com
                parpool(3);
                parfor perm_idx=1:nPerm % ----------permutation
                    [null_mean(:,perm_idx),null_median(:,perm_idx)] = permute_ISC(stcSize,condiG_label,corr_buffer,tril_idx,hypo,Brainout);
                end
                delete(gcp('nocreate'))
            else
                o = ' ';
                for perm_idx=1:nPerm % ----------permutation
                    n = sprintf('permutation [%1.1f%%]',perm_idx/nPerm*100);
                    if string(o) ~= string(n), disp(n); end
                    [null_mean(:,perm_idx),null_median(:,perm_idx)] = permute_ISC(stcSize,condiG_label,corr_buffer,tril_idx,hypo,Brainout);
                    o = n;
                end
            end
            % the ratio of how many permutation mean(null mean) greater than
            % original mean(effect_mean) (one-side right-tail)
            vol_idx = single(1:stcSize(1));
            vol_idx = vol_idx(~ismember(vol_idx,Brainout));
            for v_idx=vol_idx
                zm_perm_p(v_idx)=length(find(null_mean(v_idx,:)>effect_mean(hypothesis_idx,v_idx)))./nPerm;
                zmd_perm_p(v_idx)=length(find(null_median(v_idx,:)>effect_median(hypothesis_idx,v_idx)))./nPerm;
            end
        
            fprintf('\n');
       
            % save meanP and medianP 
            V_header_zm =  V_data(1);  V_header_zm.fname = sprintf('%s_mean_p.nii', outputname);
            V_header_zm.dt(1)=16;
            V_header_zmd =  V_data(1); V_header_zmd.fname = sprintf('%s_median_p.nii', outputname);
            V_header_zmd.dt(1)=16;
            
            spm_write_vol(V_header_zm,reshape(zm_perm_p, vX, vY, vZ));
            spm_write_vol(V_header_zmd,reshape(zmd_perm_p, vX, vY, vZ));
            
            
            % save 1-meanP and 1-medianP    
            V_header_zm =  V_data(1);  V_header_zm.fname = sprintf('%s_mean_lmp.nii', outputname);
            V_header_zm.dt(1)=16;
            V_header_zmd =  V_data(1); V_header_zmd.fname = sprintf('%s_median_lmp.nii', outputname);
            V_header_zmd.dt(1)=16;
            
            
            spm_write_vol(V_header_zm,reshape(1-zm_perm_p, vX, vY, vZ));
            spm_write_vol(V_header_zmd,reshape(1-zmd_perm_p, vX, vY, vZ));
                 
        end
        
        if ~exist(fullfile(resultfolder,fname{hypothesis_idx},'thres'),'dir')
            mkdir(fullfile(resultfolder,fname{hypothesis_idx},'thres'))
        end
        
        ptype = {'lmp','p'};
        for hypothesis_idx=1:length(hypothesis)
            for i = 1:length(ptype)
                hypothesis_z= fullfile(resultfolder,fname{hypothesis_idx},['ISC_hypothesis',sprintf('%02d',hypothesis_idx),'_mean_z.nii']);
                hypothesis_p= fullfile(resultfolder,fname{hypothesis_idx},['ISC_hypothesis',sprintf('%02d',hypothesis_idx),'_mean_',ptype{i},'.nii']);
                
                % read z value and p value
                V = spm_vol(hypothesis_z);
                thisMask = spm_read_vols(spm_vol(in.maskdir));
                thisZ = spm_read_vols(spm_vol(hypothesis_z));
                thisP = spm_read_vols(spm_vol(hypothesis_p));
            
                % use p value to correct the null hypothesis
                q=0.05;    
                [~,h]=fdr(thisP,q); % FDR
                
                % correct the z value
                thisZ(h == 0) = 0; % FDR
                thisZ(thisP>=p_threshold) = 0;
                thisZ(thisMask == 0) = 0;
                
                % save mean z value
                % p_display = split(string(p_threshold),'.');
                newV =  V;  newV.fname = fullfile(resultfolder, fname{hypothesis_idx},'thres',['ISC_hypothesis',sprintf('%02d',hypothesis_idx),'_mean_',ptype{i},'_',num2str(p_threshold),'.nii']);
                spm_write_vol(newV,thisZ);
                figpath = fullfile(resultfolder,'result');if ~exist(figpath,'dir'),mkdir(figpath);end
                filepath = fullfile(resultfolder, fname{hypothesis_idx},'thres');
                filename = ['ISC_hypothesis',sprintf('%02d',hypothesis_idx),'_mean_',ptype{i},'_',num2str(p_threshold)];
                xjview_save(filepath,filename, ...
                            figpath,[fname{hypothesis_idx},'_mean_',ptype{i},'_',num2str(p_threshold)]);
            end
        end
        hypo = array2table(cell2mat(hypothesis));
        hypo.Properties.VariableNames = condiG_label;
        writetable(hypo,fullfile(resultfolder, fname{hypothesis_idx},'thres','hypothesis.txt'));
    end
end

function [null_mean,null_median] = permute_ISC(stcSize,condition_label,corr_buffer,tril_idx,hypo,Brainout)
    null_mean = zeros(stcSize(1),1);
    null_median = zeros(stcSize(1),1);
    n_subj=stcSize(3); % subject number
    H = hypo(:)';
    stcS2 = stcSize(2);
    vol_idx = single(1:stcSize(1));
    vol_idx = vol_idx(~ismember(vol_idx,Brainout));

    for v_idx=vol_idx % -------------voxel
        % if(mod(v_idx,100)==0)
        %     fprintf('permutation [%03d|%03d]...[%1.1f%%]...\r',perm_idx,nPerm,v_idx/size(stc,1)*100);
        % end
  
        cc_perm = zeros(length(condition_label),nchoosek(n_subj,2));
        for cond_idx=1:length(condition_label) % ---------------coditional
            perm_swb_subj=randperm(n_subj); % random suffle the subject number, ex:1,2,3,4,5 --> 5,2,4,3,1 (perm_swb_subj)
            perm_swb_subj=perm_swb_subj(1:round(length(perm_swb_subj)/2)); % get top half subject, ex:5,2,4,3,1 --> 5,2,4 (perm_swb_subj)
            perm_swb_vec=ones(n_subj,1);
            perm_swb_vec(perm_swb_subj)=-1; % ex:-1,-1,-1,1,1 (perm_swb_vec)
            perm_swb_mat=perm_swb_vec*perm_swb_vec'; % ex: -1,-1,-1,1,1 -->  1  1  1 -1 -1
                                                     %                       1  1  1 -1 -1
                                                     %                       1  1  1 -1 -1
                                                     %                      -1 -1 -1  1  1
                                                     %                      -1 -1 -1  1  1
        
            %data_perm=etc_phasescramble(data,'dim',1);
        
            tmp_perm=tril(squeeze(corr_buffer(cond_idx,v_idx,:,:)).*perm_swb_mat,-1); % let some correlation value to negative and get lower triangular area
            rr_perm=tmp_perm(tril_idx); % lower triangular to 1D array
            z_perm=0.5.*log((1+rr_perm)./(1-rr_perm))./(1/sqrt(stcS2/2.34-3));% fisher Z
        
            %cc(v_idx,:)=z(:)';
            cc_perm(cond_idx,:)=z_perm(:)'; %cc_perm -- save all condition data(condition*C(subject,2))
        end
        
        tmp=H*cc_perm(:,:); 
        null_mean(v_idx,1)=mean(tmp); % store all null mean and null median -- num_voxel*num_permutation(271633*5000)
        null_median(v_idx,1)=median(tmp);
    end
end
function out = getVarName(var)
    out = inputname(1);
end

function rr = get_corr_buffer(matSize,tril_array)
    rr = ones(matSize);
    tmp = tril(rr,-1);
    rr(find(tmp)) = tril_array;
    rr = rr';
    rr(find(tmp)) = tril_array;
    rr = rr';
end

function xjview_save(filepath,filename,figpath,figname)
    xjview(fullfile(filepath,[filename,'.nii']));
    saveas(gcf,fullfile(figpath,[figname,'.fig']))
    close gcf;
    close(findall(groot,'type','figure','tag','Msgbox_Warning Dialog'))
end