
%% modularity analyis
modir                                                   = string(ls('~\metric_all\unweighted\*modu.nii'));
modularitydir                                           = string(ls('~\metric_all\unweighted\*modu.txt'));

cd                                                      '~\metric_all\unweighted'
ModularityParas                                         = paras;
for sub = 1:N_sub
    subject                                             = cell2mat(subjects(sub));
    ses                                                 = cell2mat(sess(sub));
    disp                                                (strcat(num2str(N_sub),'/',num2str(sub)))
    
    fildID                                              = fopen(modularitydir(contains(modularitydir,[subject '_' ses])~=0));
    moduletxt                                           = textscan(fildID,'%s');
    mind                                                = find(contains(moduletxt{1},'Final')~=0);
    Mind                                                = mind(1)+3;
    Nind                                                = mind(1)+8;
    if isnan(str2double(moduletxt{1}{Mind,1}))
        Mind                                            = Mind+1;
        Nind                                            = Nind+2;
    end
    
    ModularityParas.modularity(sub)                     = str2double(moduletxt{1}{Mind,1});
    ModularityParas.moduleNum(sub)                      = str2double(moduletxt{1}{Nind,1});
    
    fclose(fildID);
end
writetable                                              (ModularityParas,'~\paras_module.csv')

%% within and between FC

for sub = 1:N_sub
    subject                                             = cell2mat(subjects(sub));
    ses                                                 = cell2mat(sess(sub));
    disp                                                (strcat(num2str(N_sub),'/',num2str(sub)))
    
    % load module nii
    Ind                                                 = ~cellfun(@isempty,(regexp(modir,[subject '_' ses])));
    Indmodir                                            = strtrim(char(modir(Ind)));
    
    modulenii                                           = load_nii(['~\metric_all\unweighted\' Indmodir]);
    moduleImg                                           = modulenii.img;
    moduleInd                                           = moduleImg(ind);
    
    % load fmri nii
    if contains(subject,'sub')
        dpath                                           = '~\qlli\dHCPdata';
        fmripath                                        = strcat(dpath,'\',subject,'\',ses,'\funtemplate\',subject,'_',ses,'_bold_scrubbed.nii.gz');
    else
        dpath                                           = '~\qlli\BCPdata\BCP_func_03yr';
        fmripath                                        = strcat(dpath,'\',subject,'\',subject,'_',ses,'\func\funtemplate\',subject,'_',ses,'_fMRI_AP_scrubbed.nii.gz');
    end
    fmrinii                                             = load_nii(fmripath);
    
    t                                                   = size(fmrinii.img,4);
    fimg                                                = reshape(fmrinii.img,x*y*z,t);
    fmri                                                = fimg(ind,:)';
    
    % within FC
    uniqueMind                                          = unique(moduleInd);
    for mo = 1:size(uniqueMind,1)
        mofmri                                          = fmri(:,moduleInd==uniqueMind(mo));
        withFC                                          = corr(mofmri);
        ZwithFC                                         = atanh(withFC);        
        ZwithFC(ZwithFC<0)                              = 0;
        ZwithFC(isinf(ZwithFC)|isnan(ZwithFC))          = 0;
        all_withFC(mo)                                  = mean(mean(ZwithFC));
    end
    sub_WBFC(sub,1)                                     = mean(all_withFC(all_withFC~=0));
    
    % between FC
    for mo = 1:size(uniqueMind,1)
        mofmri                                          = fmri(:,moduleInd==uniqueMind(mo));
        other_mofmri                                    = fmri(:,moduleInd~=uniqueMind(mo));
        betweenFC                                       = corr(mofmri,other_mofmri);
        ZbetweenFC                                      = atanh(betweenFC);        
        ZbetweenFC(ZbetweenFC<0)                        = 0;
        ZbetweenFC(isinf(ZbetweenFC)|isnan(ZbetweenFC)) = 0;
        all_betweenFC(mo)                               = mean(mean(ZbetweenFC));
    end
    sub_WBFC(sub,2)                                     = mean(all_betweenFC);
    
end
writematrix                                             (sub_WBFC,'~\within_between_FC_desktop.csv')

%% pc analysis
pcdir                                                   = string(ls('~\metric_all\unweighted\*pc.nii'));

for sub = 1:N_sub
    subject                                             = cell2mat(subjects(sub));
    ses                                                 = cell2mat(sess(sub));
    disp                                                (strcat(num2str(N_sub),'/',num2str(sub)))
    
    pcnii                                               = load_nii(char(strrep(pcdir(contains(pcdir,[subject '_' ses])~=0),' ','')));
    rpc                                                 = reshape(pcnii.img,x*y*z,1);
    Npc                                                 = rpc(ind);
    pc(:,sub)                                           = Npc;
end
Rind                                                    = all(pc == 0,1);
pcparas                                                 = paras(~Rind,:);
writetable                                              (pcparas, '~\pcparas.csv')
writematrix                                             (pc,'~\pc_allvoxels.csv')

% pca for orignal pc
npc                                                     = readmatrix('~\pc_allvoxels.csv');
[coeff,score,latent,tsquared,explained,mu]              = pca(npc');
L1                                                      = coeff(:,1);

Lmap                                                    = zeros(size(gimg,1),1);
Lmap(ind)                                               = L1;
Lmap                                                    = reshape(Lmap, x, y, z);
dimg                                                    = gm_mask;
dimg.img                                                = Lmap;
save_nii                                                (dimg, '~\pc_loading.nii');

nbin                                                    = 10;
[Y,E]                                                   = discretize(L1,nbin);
for nb = 1:nbin
    bpc(nb,:)                                           = mean(npc(Y==nb,:),'omitnan');
end
writematrix                                             (bpc,'~\pc_bin10.csv');