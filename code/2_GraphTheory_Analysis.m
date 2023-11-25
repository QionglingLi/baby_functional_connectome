

paras                                               	= readtable('~\basicparas930.csv');
N_sub                                                   = size(paras,1);
subjects                                                = paras.sub_id;
sess                                                    = paras.ses_id;

% load mask for generating new image
gm_mask                                                 = load_nii('~\atlas\UNC-BCP_4D_Infant_Brain_Volumetric_Atlas_v1\6Month\BCP-06M-GM_2mm.nii.gz');
[x, y, z]                                               = size(gm_mask.img);
gimg                                                    = reshape(gm_mask.img,x*y*z,1);
ind                                                     = find(gimg>0.5);
%% organize graph theoretical measures
cd                                                      '~\metric_all\unweighted'
Cpdir                                                   = string(ls('~\metric_all\unweighted\*cp.txt'));
effdir                                                  = string(ls('~\metric_all\unweighted\*eff.nii'));
Cpniidir                                                = string(ls('~\metric_all\unweighted\*cp.nii'));

Sparas                                                  = paras;
for sub = 1:N_sub
    subject                                             = cell2mat(subjects(sub));
    ses                                                 = cell2mat(sess(sub));
    disp                                                (strcat(num2str(N_sub),'/',num2str(sub)))
    
    Cp                                                  = load(Cpdir(contains(Cpdir,[subject '_' ses])~=0));
    cpnii                                               = load_nii(char(strrep(Cpniidir(contains(Cpniidir,[subject '_' ses])~=0),' ','')));
    rcp                                                 = reshape(cpnii.img,x*y*z,1);
    Ncp                                                 = rcp(ind);
    cp_allvoxels(:,sub)                                 = Ncp;
    
    effnii                                              = load_nii(char(strrep(effdir(contains(effdir,[subject '_' ses])~=0),' ','')));
    reff                                                = reshape(effnii.img,x*y*z,1);
    Ne                                                  = reff(ind);
    eff_allvoxels(:,sub)                                = Ne;
    Geff                                                = mean(Ne);
    
    Sparas.Cp(sub)                                      = Cp(1);
    Sparas.GlobalEff(sub)                               = Geff;
    
end


writetable                                              (Sparas,'~\Sparas.csv')
writematrix                                             (cp_allvoxels,'~\cp_allvoxels.csv')
writematrix                                             (eff_allvoxels,'~\eff_allvoxels.csv')

%% pca for principal axis
ncp                                                     = readmatrix('~\cp_allvoxels.csv');
[coeff,score,latent,tsquared,explained,mu]              = pca(ncp');
L1                                                      = coeff(:,1);

Lmap                                                    = zeros(size(gimg,1),1);
Lmap(ind)                                               = L1;
Lmap                                                    = reshape(Lmap, x, y, z);
dimg                                                    = gm_mask;
dimg.img                                                = Lmap;
save_nii                                                (dimg, '~\cp_loading.nii');

% divide measures into 10 bins
nbin                                                    = 10;
[Y,E]                                                   = discretize(L1,nbin);
for nb = 1:nbin
    bcp(nb,:)                                           = mean(ncp(Y==nb,:),'omitnan');
end
writematrix                                             (bcp,'~\cp_bin10.csv');

% resorted ncp
[B,I]                                                   = sort(L1,'descend');
resorted_cp                                             = ncp(I,:);
figure;imagesc(resorted_cp')


