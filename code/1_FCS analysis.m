%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FCS analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load subjects information
paras                                               = readtable('~\basicparas930.csv');
N_sub                                               = size(paras,1);

% load mask
gm                                                  = spm_vol('~\UNC-BCP_4D_Infant_Brain_Volumetric_Atlas_v1\6Month\BCP-06M-GM_2mm.nii');
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
        dpath                                       = '~\qlli\dHCPdata';
        fmripath                                    = strcat(dpath,'\',subject,'\',ses,'\funtemplate\',subject,'_',ses,'_bold_scrubbed.nii.gz');
    else
        dpath                                       = '~\qlli\BCPdata\BCP_func_03yr';
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
save                                                ('~\output_distance\FCs.mat', 'FCs')
%% 2. extract peak FCs
pcor                                                = [-16 -60 18; -42 -14 42; -8 32 20; -18 28 24; -50 -2 -12; -12 10 -2]; 
CoordinateMatrix                                    = round(CoordinateMatrix');
    
for p = 1:size(pcor,1)
    peak_coor                                       = pcor(p,:);
    [~,peak_ind]                                    = ismember(peak_coor,CoordinateMatrix,'rows');
    mask_ind                                        = find(peak_ind==ind); % find peak location in gm mask
    peakFCs(p,:)                                    = FCs(mask_ind,:);        
end

writematrix                                         (peakFCs,'~\output_distance\FCS_k=3\peakFCs.csv')

%% 3.show all predicted seed FC for real age
for seed=1:6
    pred_table                                      = readtable(strcat('~\output_distance\FCS_k=3\seedFC\seedFC',num2str(seed),'\seedFC_predvalue.txt'));
    pred_table                                      = sortrows(pred_table,'Var1','ascend');
    pred_FCs                                        = table2array(pred_table(:,2:end));
    
    scan_age_pred                                   = load('~\output_distance\FCs\scan_age.txt');
    [a_scan_age,aind]                               = sort(scan_age_pred);
    
    s_pred_FCs                                      = pred_FCs(:,aind);
    
    for i = [21:26, 51:55, 68:73, 212:217, 356:360, 427:431, 641:646, 928:930]%1:930
        FCsmap                                      = zeros(size(gimg,1),1);
        FCsmap(ind(pred_table.Var1))                = s_pred_FCs(:,i);%
        FCsmap                                      = reshape(FCsmap, x, y, z);
        dimg                                        = gm_mask;
        dimg.img                                    = FCsmap;
        save_nii                                    (dimg, ['F:\OneDrive - 北京师范大学\project2\output_distance\FCS_k=3\seedFC\seedFC',num2str(seed),...
                                                             '\allsub\seedFC',num2str(seed),'_pred_',num2str(i),'_',num2str(a_scan_age(i)),'w.nii']);
        
        BrainVolume                                 = strcat('F:\OneDrive - 北京师范大学\project2\output_distance\FCS_k=3\seedFC\seedFC',num2str(seed),...
                                                             '\allsub\seedFC',num2str(seed),'_pred_',num2str(i),'_',num2str(a_scan_age(i)),'w.nii');
        H_BrainNet                                  = BrainNet_MapCfg(Surf,BrainVolume,['F:\OneDrive - 北京师范大学\project2\output_distance\FCS_k=3\seedFC\seedFC',num2str(seed),'\Cfg_FC',num2str(seed),'.mat']);
        colormap(cmap)
        Gmap                                        = strcat('F:\OneDrive - 北京师范大学\project2\output_distance\FCS_k=3\seedFC\seedFC',num2str(seed),...
                                                             '\allsub\seedFC',num2str(seed),'_pred_',num2str(i),'_',num2str(a_scan_age(i)),'w.tif');
        saveas                                      (H_BrainNet,Gmap)
        clear BrainVolume H_BrainNet
        close all
    end
end


%% 4. compute seed based FC
for i = N_sub:-1:1
    disp                                            (strcat(num2str(N_sub),'/',num2str(i)))
    
    subject                                         = char(paras.sub_id(i));
    ses                                             = char(paras.ses_id(i));
    
    if contains(paras.sub_id(i),'sub')
        dpath                                       = '~\qlli\dHCPdata';
        fmripath                                    = strcat(dpath,'\',subject,'\',ses,'\funtemplate\',subject,'_',ses,'_bold_scrubbed.nii.gz');
    else
        dpath                                       = '~\qlli\BCPdata\BCP_func_03yr';
        fmripath                                    = strcat(dpath,'\',subject,'\',subject,'_',ses,'\func\funtemplate\',subject,'_',ses,'_fMRI_AP_scrubbed.nii.gz');
    end
    fmrinii                                         = load_nii(fmripath);
    
    t                                               = size(fmrinii.img,4);
    fimg                                            = reshape(fmrinii.img,x*y*z,t);
    fmri                                            = fimg(ind,:)';
    
    for p = 1:size(pcor,1)
        peak_coor                                   = pcor(p,:);
        [~,peak_ind]                                = ismember(peak_coor,CoordinateMatrix,'rows');
        mask_ind                                    = find(peak_ind==ind); % find peak location in gm mask        
        seed_fmri                                   = fmri(:,mask_ind);
        seed_FC{p,i}                                = corr(fmri,seed_fmri);
    end
    
    clear fimg fmrinii
    clear fmri
end
seed_zFC                                            = cellfun(@atanh,seed_FC,'UniformOutput',false);

for p = 1:size(pcor,1)
    seed_FC1                                        = cell2mat(seed_zFC(p,:));
    seed_FC1(seed_FC1<0)                            = 0;
    seed_FC1(isinf(seed_FC1)|isnan(seed_FC1))       = 0;
    writematrix                                     (seed_FC1,strcat('~\output_distance\FCS_k=3\seedFC\seed_FC',num2str(p),'.csv'));
end

%% 5. hub changes
% identify hubs
MEAN                                                    = mean(s_pred_FCs,1);
STD                                                     = std(s_pred_FCs,0,1);
hub_voxels                                              = s_pred_FCs>(MEAN+1.5*STD);

% read yeo7 networks
yeo7_nii                                                = load_nii('~\UNC-Infant-yeo\Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_BCP24_to_BCP06_2mm.nii.gz');
yeo7                                                    = yeo7_nii.img(ind);

for y = 1:7
    for s = 1:930
        hub_voxel_num(y,s)                               = sum(hub_voxels(:,s) & (yeo7==y));
    end
end
hub_voxel_percent                                       = hub_voxel_num./sum(hub_voxel_num,1);
writematrix                                             (hub_voxel_num,'~\output_distance\FCS_k=3\hub\hub_size_yeo7.csv');
writematrix                                             (hub_voxel_percent,'~\output_distance\FCS_k=3\hub\hub_percent_yeo7.csv');

% plot hub voxel bar
figure; bar(hub_voxel_num(:,[24 54 70 215 360 430 645 930])','stacked')

%% 6. Compute individual FCS within each bin
for i = 1:N_sub
    subject                                         = char(paras.sub_id(i));
    ses                                             = char(paras.ses_id(i));
    if contains(paras.sub_id(i),'sub')
        dpath                                       = '~\qlli\dHCPdata';
        fmripath                                    = strcat(dpath,'\',subject,'\',ses,'\funtemplate\',subject,'_',ses,'_bold_scrubbed.nii.gz');
    else
        dpath                                       = '~\qlli\BCPdata\BCP_func_03yr';
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
