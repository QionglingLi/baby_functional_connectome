%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FCS analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load subjects information
paras                                               = readtable('.\basicparas930.csv');
N_sub                                               = size(paras,1);

% load mask
gm                                                  = spm_vol('.\atlas\UNC-BCP_4D_Infant_Brain_Volumetric_Atlas_v1\6Month\BCP-06M-GM_2mm.nii');
[gm_mask,CoordinateMatrix]                          = spm_read_vols(gm);
[x, y, z]                                           = size(gm_mask);
gimg                                                = reshape(gm_mask,x*y*z,1);
ind                                                 = find(gimg>0.5);

% compute distance between paired voxels
D                                                   = single(pdist(CoordinateMatrix(:,ind)'));
clear gimg gm_mask CoordinateMatrix
ZD                                                  = single(squareform(D));
clear D

ZD_ind                                              = zeros(size(ind,1));
ZD_ind(ZD>10 & ZD<=30)                              = 1;
ZD_ind(ZD>30 & ZD<=50)                              = 2;
ZD_ind(ZD>50)                                       = 3;

%% 1.compute FCS for each subject
for i = 1:N_sub
    subject                                         = char(paras.sub_id(i));
    ses                                             = char(paras.ses_id(i));
    if contains(paras.sub_id(i),'sub')
        dpath                                       = '.\qlli\dHCPdata';
        fmripath                                    = strcat(dpath,'\',subject,'\',ses,'\funtemplate\',subject,'_',ses,'_bold_scrubbed.nii.gz');
    else
        dpath                                       = '.\qlli\BCPdata\BCP_func_03yr';
        fmripath                                    = strcat(dpath,'\',subject,'\',subject,'_',ses,'\func\funtemplate\',subject,'_',ses,'_fMRI_AP_scrubbed.nii.gz');
    end
    fmrinii                                         = load_nii(fmripath);
    
    t                                               = size(fmrinii.img,4);
    fimg                                            = reshape(fmrinii.img,x*y*z,t);   
    fmri                                            = fimg(ind,:)';
    
    clear fimg fmrinii    
    disp                                            (strcat(num2str(N_sub),'/',num2str(i)))
    FC                                              = corr(fmri);
    FC(FC<0)                                        = 0;
    FC(isinf(FC)|isnan(FC))                         = 0;   
    clear fmri
    
    FCs(i,:)                                        = sum(FC.*(ZD_ind~=0),1)./sum((ZD_ind~=0),1);    
end
save                                                ('.\output_distance\FCs.mat', 'FCs')
%% 2. cluster for sig FCS
Fmap_mask                                           = load_nii('~\output_distance\FCS\FCs_FMap_GRFcorrected_mask.nii');
Fmap_mask_img                                       = reshape(Fmap_mask.img,x*y*z,1);
Fmap_mask_img_5w                                    = Fmap_mask_img(ind);
FCs_sig                                             = FCs(logical(Fmap_mask_img_5w),:);

k = 4;
[idx, C]                                            = kmeans(FCs_sig, k);
save                                                ('~\output_distance\cluster\Cluster_idx.mat','idx')
for  i = 1:k
    Cluster_FCS(i,:)                                 = mean(FCs_sig(idx==i,:));
end
writematrix                                         (Cluster_FCS, strcat('~\output_distance\cluster\Cluster_FCS.csv'))

L1_5w                                               = zeros(numel(ind),1);
L1_5w(logical(Fmap_mask_img_5w),1)                  = idx;

Lmap                                                = zeros(size(gimg,1),1);
Lmap(ind)                                           = L1_5w;
Lmap                                                = reshape(Lmap, x, y, z);
dimg                                                = gm_mask;
dimg.img                                            = Lmap;
save_nii                                            (dimg, '~\output_distance\cluster\sigFCs_cluster_v2.nii');

%% 3. compute cluster based FC
for i = 1:N_sub
    disp                                            (strcat(num2str(N_sub),'/',num2str(i)))
    
    subject                                         = char(paras.sub_id(i));
    ses                                             = char(paras.ses_id(i));
    
    if contains(paras.sub_id(i),'sub')
        dpath                                       = '~\dHCPdata';
        fmripath                                    = strcat(dpath,'\',subject,'\',ses,'\funtemplate\',subject,'_',ses,'_bold_scrubbed.nii.gz');
    else
        dpath                                       = '~\BCPdata\BCP_func_03yr';
        fmripath                                    = strcat(dpath,'\',subject,'\',subject,'_',ses,'\func\funtemplate\',subject,'_',ses,'_fMRI_AP_scrubbed.nii.gz');
    end
    fmrinii                                         = load_nii(fmripath);
    
    t                                               = size(fmrinii.img,4);
    fimg                                            = reshape(fmrinii.img,x*y*z,t);
    fmri                                            = fimg(ind,:)';
    
    for p = 1:4
        seed_ind_sig                                = find(idx==p);
        seed_ind_5w_tmp                             = find(Fmap_mask_img_5w~=0);
        seed_ind_5w                                 = seed_ind_5w_tmp(seed_ind_sig);
        
        seed_fmri                                   = mean(fmri(:,seed_ind_5w),2);
        seed_FC{p,i}                                = corr(fmri,seed_fmri);
    end
    
    clear fimg fmrinii
    clear fmri
end
seed_zFC                                            = cellfun(@atanh,seed_FC,'UniformOutput',false);

for p = 1:4
    seed_FC1                                        = cell2mat(seed_zFC(p,:));
    seed_FC1(seed_FC1<0)                            = 0;
    seed_FC1(isinf(seed_FC1)|isnan(seed_FC1))       = 0;
    writematrix                                     (seed_FC1,strcat('~\output_distance\cluster\seed\seed_FC',num2str(p),'.csv'));
end

%% 4. ratio of F-values in each system
yeo7_nii                                                = load_nii('~\atlas\UNC-Infant-yeo\Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_BCP24_to_BCP06_2mm.nii.gz');
yeo7                                                    = yeo7_nii.img(ind);

seedname                                                 = {'Cluster1','Cluster2','Cluster3','Cluster4'};
for seed = 1:4
    Fmapnii                                              = load_nii(['~\output_distance\cluster\seed\seedFC',num2str(seed),'_FMap_GRFcorrected.nii']);
    Fimg                                                 = reshape(Fmapnii.img,x*y*z,1);
    Fmap                                                 = Fimg(ind);
    
    for ye = 1:7
        seed_F_values(ye,seed)                            = sum(Fmap & (yeo7==ye));
    end
end

SFvalue                                                  = sum(seed_F_values,1);
seed_F_percent                                           = seed_F_values./SFvalue;

f=figure; %
De_mica_spider                                          (seed_F_percent, 'Fvalues',[0 0.5], ...
                                                        {'Visual','Somatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontoparietal','Default'},...
                                                        {'Cluster1','Cluster2','Cluster3','Cluster4'},color(1:64:256,:),gca);

%% 5. hub changes
% identify hubs
MEAN                                                 = mean(s_pred_FCs,1);
STD                                                  = std(s_pred_FCs,0,1);
hub_voxels                                           = s_pred_FCs>(MEAN+1.5*STD);

for y = 1:7
    for s = 1:930
        hub_voxel_num(y,s)                           = sum(hub_voxels(:,s) & (yeo7==y));
    end
end
hub_voxel_percent                                    = hub_voxel_num./sum(hub_voxel_num,1);
writematrix                                          (hub_voxel_num,'.\output_distance\FCS_k=3\hub\hub_size_yeo7.csv');
writematrix                                          (hub_voxel_percent,'.\output_distance\FCS_k=3\hub\hub_percent_yeo7.csv');

% plot hub voxel bar
figure; 
bar(hub_voxel_num(:,[24 54 70 215 360 430 645 930])','stacked')

%% 6. Compute individual FCS within each bin
for i = 1:N_sub
    subject                                         = char(paras.sub_id(i));
    ses                                             = char(paras.ses_id(i));
    if contains(paras.sub_id(i),'sub')
        dpath                                       = '.\qlli\dHCPdata';
        fmripath                                    = strcat(dpath,'\',subject,'\',ses,'\funtemplate\',subject,'_',ses,'_bold_scrubbed.nii.gz');
    else
        dpath                                       = '.\qlli\BCPdata\BCP_func_03yr';
        fmripath                                    = strcat(dpath,'\',subject,'\',subject,'_',ses,'\func\funtemplate\',subject,'_',ses,'_fMRI_AP_scrubbed.nii.gz');
    end
    fmrinii                                         = load_nii(fmripath);
    
    t                                               = size(fmrinii.img,4);
    fimg                                            = reshape(fmrinii.img,x*y*z,t);   
    fmri                                            = fimg(ind,:)';
    
    clear fimg fmrinii    
    disp                                            (strcat(num2str(N_sub),'/',num2str(i)))
    FC                                              = corr(fmri);
    FC(FC<0)                                        = 0;
    FC(isinf(FC)|isnan(FC))                         = 0;   
    clear fmri
    
    for b = 1:3
       FCs_bin{b}(i,:)                              = sum(FC.*(ZD_ind==b),1)./sum((ZD_ind==b),1); 
    end
    
    clear FC
end

save                                                ('~\output_distance\FCsbins\FCs_bin.mat', 'FCs_bin')
