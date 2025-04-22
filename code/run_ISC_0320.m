  
% close all; clear all;
%% specify parameters
% define path
% addpath('D:\spm12')   % location of spm 
% addpath('D:\fhlin_toolbox-master');  % location of the fhlin-toolbox
datafolder = 'C:\Users\user\Desktop\jeng\fmri_ISC\ISC\data\';  % Data should be sorted into subfolders according to the condition
resultfolder  = 'C:\Users\user\Desktop\jeng\fmri_ISC\ISC\output'; % Output directoy
maskdir = 'C:\Users\user\Desktop\jeng\fmri_ISC\ISC\mask\rmask.volume.brainmask.nii';  % location of the mask

%define subject and condition
%subj_path= {'swauISO_01', 'swauISO_02'};

% all subject
% for i = 1:length(n)
%     subj_path = {}
% end


subj_path= {'001', '002', '003', '004', '005', '006', '007', '008', '009', '010'};
%{
n=15 low grief
subj_path= {'ISO_04', 'ISO_05', 'ISO_06', 'ISO_09', 'ISO_10', 'ISO_11', 'ISO_14', 'ISO_15', 'ISO_16', 'ISO_17', ......
            'ISO_23', 'ISO_24',  'ISO_27', 'ISO_28', 'ISO_29'};
n=15 high grief
subj_path= { 'ISO_01', 'ISO_02', 'ISO_03', 'ISO_07', 'ISO_08','ISO_12', 'ISO_13', 'ISO_18', 'ISO_19', 'ISO_20', ......            
'ISO_21', 'ISO_22', 'ISO_25', 'ISO_26', 'ISO_30' };
%}

% condition_label={'bf_sad', 'bf_happy', 'af_sad', 'af_happy'}; % exp3_bf_af
condition_label={'rest'};
subdir = dir(fullfile(datafolder,char(condition_label)));
subj_path = {subdir.name};
subj_path = subj_path(3:end);
% remove the first & last few volumes (unit = TR)
remove_volume_start = 1;
remove_volume_end =1;

%define analysis of interest
hypothesis={
1  % happy-sad
    };

% exp3_bf_af
% hypothesis={
%     [1 0 0 0];  % bf_sad
%     [0 1 0 0];  % bf_happy
%     [0 0 1 0];  % af_sad
%     [0 0 0 1];  % af_happy
%     [1 0 -1 0]; % bf_sad - af_sad
%     [-1 0 1 0]; % af_sad - bf_sad
%     [0 1 0 -1]; % bf_happy - af_happy
%     [0 -1 0 1]; % af_happy - bf_happy  
%     [1 -1 0 0]; % bf_sad - bf_happy
%     [-1 1 0 0]; % bf_happy - bf_sad
%     };


%specify the number of permutation and stats threshold
nPerm = 5000;
p_threshold = 0.05;

% import mask
V_mask = spm_vol(maskdir); 

mask = spm_read_vols(V_mask);
Brainout = find(mask==0);
confound_polynomial_order=2;
confound_sinusoidal_order=3;


%% PART 1: calculate inter-subject correlation for each voxel
for condition_idx=1:length(condition_label) % for every condition

    datapath = fullfile(datafolder,condition_label{condition_idx});
    subject=subj_path;
    stc=[];

    for subj_idx=1:length(subject)
        niipath = dir(fullfile(datapath, [subject{subj_idx}, '*']));
        V_data = spm_vol(fullfile(datapath,niipath.name));  % read nii file attribute(V_data)
        data = spm_read_vols(V_data); % get 4-D data(data)
        nX = size(data,1); nY = size(data,2); nZ = size(data,3); % get image n*m and slice number(nX--n,nY--m,mZ--slice)
        disp(['loading ', fullfile(datapath,niipath.name)])

        % 1 subject = 4D img to 2D => put every subject in a matrix   ex.img 2*2*2    *10  *20 --> stc 8  *10  *20
        %                                                                    n*m*slice*time*subject    img*time*subject
        for tt = 1:size(data,4) 
            thisV = data(:,:,:,tt); % find every time point image
            thisV(Brainout) = NaN; % set vales outside the mask to NaN
            stc(:,tt,subj_idx) = reshape(thisV,nX*nY*nZ,1); % set a image data(3D) to a column(1D) ex. img 2*2*2 --> stc 8*1 (stc--voxel*slice*subject)
        end
    end

    ii=find(isnan(stc(:))); % find NaN
    stc(ii)=randn(size(ii)).*eps; % set NaN equal to a random number(base on normal distribute)*2^-56

    %remove the first and last few time points
    stc=stc(:,remove_volume_start:end-remove_volume_end,:);

    %remove global mean
    for s_idx=1:size(stc,3) % for every subject
        tmp=squeeze(stc(:,:,s_idx)); % find every subject image (allvoxel*timepoint)
        ga=mean(tmp,1); % find every img mean(ga)

        tmp=tmp-tmp*ga'*inv(ga*ga')*ga; %?? remove mean

        stc(:,:,s_idx)=tmp; % set new img in same space
    end

    timeVec=[1:size(stc,2)]'; % time vector(1:number of slice)
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

    fprintf('\n');

    for v_idx=1:size(stc,1) % find every voxel 
        if(mod(v_idx,100)==0)
            fprintf('[%1.1f%%]...\r',v_idx/size(stc,1)*100); % print process rate 
        end
        if(~isempty(D_prep))
            data=squeeze(stc(v_idx,:,:))-D_prep*squeeze(stc(v_idx,:,:)); % data(voxel value) -- slice*subject(for every voxel)
        else
            data=squeeze(stc(v_idx,:,:));
        end
        ii=find(abs(data)<eps);
        data(ii)=randn(size(ii)).*eps; % if voxel value is lower than 2^-52 set it a new value

        %save all correlation coef.
        if((condition_idx==1)&&(v_idx==1))
            % corr_buffer(condition*voxel*subject*subject)
            corr_buffer=zeros(length(condition_label),size(stc,1),size(stc,3),size(stc,3)); % set a buffer to restore the correlation value
        end
        corr_buffer(condition_idx,v_idx,:,:)=corrcoef(data); % find evrey voxel in subject correlation
        %
        % condition*voxel*subject*subject
        %                 (corrlation matrix)

        tmp=tril(corrcoef(data),-1); % lower trangular part of matrix

        rr=tmp(find(tmp)); % tmp -- subject*subject(upper trangular part is zero) set tmp to (C(subject,2)*1) ---2D to 1D
        z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3)); % compute z score
        zm=mean(z);
        zmd=median(z);

        if((condition_idx==1)&&(v_idx==1))
            cc=zeros(length(condition_label),size(stc,1),length(z(:))); %  cc--(set z in here) (condition,voxel,C(subject,2))
        end
        cc(condition_idx,v_idx,:)=z(:)'; % set z
        if(find(isnan(z))) keyboard; end % if z have nan give control to keyboard

        ccm(v_idx)=zm; % store every mean in z
        ccmd(v_idx)=zmd; % store every median in z

    end
    fprintf('\n');
    if ~exist(sprintf('%s/pairs', resultfolder))
        mkdir(sprintf('%s/pairs', resultfolder))
    end

    if ~exist(sprintf('%s/group', resultfolder))
        mkdir(sprintf('%s/group', resultfolder))
    end

    if ~exist(sprintf('%s/pairs/%s', resultfolder,condition_label{condition_idx}))
        mkdir(sprintf('%s/pairs/%s', resultfolder,condition_label{condition_idx}))
    end

    for volume_idx = 1:size(cc,3) % save every pair
        V2save = reshape(cc(condition_idx,:,volume_idx),nX,nY,nZ); % V2save -- save all voxel z value in a .nii (all file number will be C(subject,2))
        V_header = V_data(1);
        V_header.fname = [sprintf('%s/pairs/%s/ISC_%s_pair-%02d.nii',resultfolder, condition_label{condition_idx}, condition_label{condition_idx}, volume_idx)];
        spm_write_vol(V_header,V2save);
    end

    V_header_ccm =  V_data(1);  V_header_ccm.fname = sprintf('%s/group/ISC_%s_mean.nii', resultfolder, condition_label{condition_idx});
    V_header_ccm.dt(1)=16;
    V_header_ccmd =  V_data(1); V_header_ccmd.fname = sprintf('%s/group/ISC_%s_median.nii', resultfolder, condition_label{condition_idx});
    V_header_ccmd.dt(1)=16;

    spm_write_vol(V_header_ccm,reshape(ccm, nX, nY, nZ));% save mean data
    spm_write_vol(V_header_ccmd,reshape(ccmd, nX, nY, nZ)); % save median data
end

%% PART 2: calculate the analysis effect and run a permutation test
if ~exist(sprintf('%s/stats', resultfolder))
    mkdir(sprintf('%s/stats', resultfolder))
end
for hypothesis_idx=1:length(hypothesis)  % set every hypothesis mean and median and save to .nii
    for v_idx=1:size(stc,1)
        tmp=hypothesis{hypothesis_idx}(:)'*squeeze(cc(:,v_idx,:)); %  cc--(set z in here) (condition,voxel,C(subject,2))
        effect_mean(hypothesis_idx,v_idx)=mean(tmp); % set hypothesis mean and median
        effect_median(hypothesis_idx,v_idx)=median(tmp);
    end
    
    % save meanz and medianz value
    V_header_mean =  V_data(1);  V_header_mean.fname = sprintf('%s/stats/ISC-hypothesis%02d_mean_z.nii', resultfolder, hypothesis_idx);
    V_header_mean.dt(1)=16;
    V_header_median =  V_data(1); V_header_median.fname = sprintf('%s/stats/ISC-hypothesis%02d_median_z.nii', resultfolder, hypothesis_idx);
    V_header_median.dt(1)=16;
    
    spm_write_vol(V_header_mean,reshape(effect_mean(hypothesis_idx,:), nX, nY, nZ));
    spm_write_vol(V_header_median,reshape(effect_median(hypothesis_idx,:), nX, nY, nZ));

end

for hypothesis_idx=1:length(hypothesis) % ---------hypothesis
    
    outputname=sprintf('%s/stats/ISC-hypothesis%02d',resultfolder,hypothesis_idx);

    zm_perm_p=zeros(size(ccm));
    zmd_perm_p=zeros(size(ccmd));
    n_perm=nPerm;
    tril_idx=find(tril(ones(size(stc,3)),-1)); % find lower triangular part index
    for perm_idx=1:n_perm % ----------permutation 
        for v_idx=1:size(stc,1) % -------------voxel
            if(mod(v_idx,100)==0)
                fprintf('permutation [%03d|%03d]...[%1.1f%%]...\r',perm_idx,n_perm,v_idx/size(stc,1)*100);
            end


            for condition_idx=1:length(condition_label) % ---------------coditional
                n_subj=size(data,2); % subject number
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

                tmp_perm=tril(squeeze(corr_buffer(condition_idx,v_idx,:,:)).*perm_swb_mat,-1); % let some correlation value to negative and get lower triangular area
                rr_perm=tmp_perm(tril_idx); % lower triangular to 1D array
                z_perm=0.5.*log((1+rr_perm)./(1-rr_perm))./(1/sqrt(size(data,1)/2.34-3));
                zm_perm=mean(z_perm);
                zmd_perm=median(z_perm);

                %cc(v_idx,:)=z(:)';
                cc_perm(condition_idx,:)=z_perm(:)'; %cc_perm -- save all condition data(condition*C(subject,2))
                ccm_perm(condition_idx)=zm_perm;
                ccmd_perm(condition_idx)=zmd_perm;
            end

            tmp=hypothesis{hypothesis_idx}(:)'*cc_perm(:,:); 
            null_mean(v_idx,perm_idx)=mean(tmp); % store all null mean and null median -- num_voxel*num_permutation(271633*5000)
            null_median(v_idx,perm_idx)=median(tmp);
        end
    end
    
    % the ratio of how many permutation mean(null mean) greater than
    % original mean(effect_mean) (one-side right-tail)
    for v_idx=1:size(stc,1)
        zm_perm_p(v_idx)=length(find(null_mean(v_idx,:)>effect_mean(hypothesis_idx,v_idx)))./n_perm;
        zmd_perm_p(v_idx)=length(find(null_median(v_idx,:)>effect_median(hypothesis_idx,v_idx)))./n_perm;
    end


    fprintf('\n');

    
    % save meanP and medianP 
    V_header_zm =  V_data(1);  V_header_zm.fname = sprintf('%s_mean_p.nii', outputname);
    V_header_zm.dt(1)=16;
    V_header_zmd =  V_data(1); V_header_zmd.fname = sprintf('%s_median_p.nii', outputname);
    V_header_zmd.dt(1)=16;
    
    spm_write_vol(V_header_zm,reshape(zm_perm_p, nX, nY, nZ));
    spm_write_vol(V_header_zmd,reshape(zmd_perm_p, nX, nY, nZ));
    
    
    % save 1-meanP and 1-medianP    
    V_header_zm =  V_data(1);  V_header_zm.fname = sprintf('%s_mean_lmp.nii', outputname);
    V_header_zm.dt(1)=16;
    V_header_zmd =  V_data(1); V_header_zmd.fname = sprintf('%s_median_lmp.nii', outputname);
    V_header_zmd.dt(1)=16;
    
    
    spm_write_vol(V_header_zm,reshape(1-zm_perm_p, nX, nY, nZ));
    spm_write_vol(V_header_zmd,reshape(1-zmd_perm_p, nX, nY, nZ));
       
    
    
end


if ~exist(sprintf('%s/thresholded', resultfolder))
    mkdir(sprintf('%s/thresholded', resultfolder))
end


for hypothesis_idx=1:length(hypothesis)

    hypothesis_z= sprintf('%s/stats/ISC-hypothesis%02d_mean_z.nii', resultfolder, hypothesis_idx);
    hypothesis_p= sprintf('%s/stats/ISC-hypothesis%02d_mean_p.nii', resultfolder, hypothesis_idx);
    
    % read z value and p value
    V = spm_vol(hypothesis_z);
    thisMask = spm_read_vols(spm_vol(maskdir));
    thisZ = spm_read_vols(spm_vol(hypothesis_z));
    thisP = spm_read_vols(spm_vol(hypothesis_p));

    % use p value to correct the null hypothesis
    q=0.05;    
    h=fdr(thisP,q); % FDR
    
    % correct the z value
    thisZ(h == 0) = 0; % FDR
    thisZ(thisP>=p_threshold) = 0;
    thisZ(thisMask == 0) = 0;
    
    % save mean z value
    p_display = split(string(p_threshold),'.');
    newV =  V;  newV.fname = sprintf('%s/thresholded/ISC-hypothesis%02d_mean_p_%s.nii', resultfolder, hypothesis_idx,string(p_threshold));
    spm_write_vol(newV,thisZ);
end



