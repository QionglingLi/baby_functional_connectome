
%% import brainspan atlas data
expression_matrix                               = load('~\BrainSpanAtlas\downloaded\genes_matrix_csv\expression_matrix.csv');
% remove all zeros
rowsWithZeros                                   = all(expression_matrix(:, 2:end) == 0, 2);
zeroRows                                        = find(rowsWithZeros);

rows_metadata                                   = readtable('~\BrainSpanAtlas\downloaded\genes_matrix_csv\rows_metadata.csv');
columns_metadata                                = readtable('~\BrainSpanAtlas\downloaded\genes_matrix_csv\columns_metadata.csv');
GeneSymbol                                      = rows_metadata.gene_symbol;

expression_matrix(zeroRows,:)                   = [];
rows_metadata(zeroRows,:)                       = [];
GeneSymbol(zeroRows,:)                          = [];

% % check duplicate genes
% [unique_strings,uniqueSIndex]           = unique(GeneSymbol,'stable');
% duplicate_positions = {};
% 
% for i = 1:length(unique_strings)
%     positions = find(strcmp(GeneSymbol, unique_strings{i}));
%     
%     if numel(positions) > 1
%         cleanStr = matlab.lang.makeValidName(unique_strings{i});
%         duplicate_positions.(cleanStr) = positions;
%     end
% end
% postion1 = expression_matrix(duplicate_positions.AC008069_1(1),:);
% postion2 = expression_matrix(duplicate_positions.AC008069_1(2),:);
% rows_metadata(duplicate_positions.AC008069_1(1),:)
% rows_metadata(duplicate_positions.AC008069_1(2),:)
% 
% expression_matrix                       = expression_matrix(uniqueSIndex,:);
% rows_metadata                           = rows_metadata(uniqueSIndex,:);

age                                             = columns_metadata.week;
gender                                          = columns_metadata.gender;
Gender                                          = zeros(size(gender,1),1);
Gender(ismember(gender,'M'))                    = 1;

structure_id                                    = columns_metadata.structure_id;
load                                            ('~\Genes\Sig.txt')
load                                            ('~\Genes\NonSig.txt')

SigIndex                                        = ismember(structure_id, Sig);
NonSigIndex                                     = ismember(structure_id, NonSig);

%% statistical analysis between sig and nonsig
SigIndex                                        = ismember(structure_id, Sig);
NonSigIndex                                     = ismember(structure_id, NonSig);

SigAgeMatchedIndex                              = logical(SigIndex.*(age>=28&age<=184)); % keep samples aged between 28 pmw and 3 yrs
NonSigAgeMatchedIndex                           = logical(NonSigIndex.*(age>=28&age<=184));

SigAge                                          = age(SigAgeMatchedIndex, 1);
NonSigAge                                       = age(NonSigAgeMatchedIndex, 1);

ExpressionData_Sig                              = expression_matrix( :, [ 0; SigAgeMatchedIndex ] > 0.5 )';
ExpressionData_NonSig                           = expression_matrix( :, [ 0; NonSigAgeMatchedIndex ] > 0.5 )';

% extract removed genes' id
removeid                                        = zeros(size(ExpressionData_Sig,2),1);
for i = 1:size(ExpressionData_Sig,2)
    tmp_sig                                     = unique(SigAge(ExpressionData_Sig(:,i)~=0));% check the number of age samples 
    tmp_nonsig                                  = unique(NonSigAge(ExpressionData_NonSig(:,i)~=0));
    if numel(tmp_sig)<6 | numel(tmp_nonsig)<6
    removeid(i)                                 = 1;
    end
end
removeid                                        = logical(removeid);
writecell                                       (unique(rows_metadata.gene_symbol(~removeid)),'~\Genes\enrich\AllGenes_6timepoint.csv')

% % plot gene expression trajectories
% x1                                              = SigAge;
% y1                                              = ExpressionData_Sig(:,102);
% 
% x2                                              = NonSigAge;
% y2                                              = ExpressionData_NonSig(:,102);
% 
% span                                            = 0.3;  % Span can be in the range 0 to 1
% smoothed_y1                                     = smooth(x1, y1, span, 'rlowess');  % 'rloess' could also be 'rlowess'
% smoothed_y2                                     = smooth(x2, y2, span, 'rlowess'); 

% % Plot the data
% figure;
% plot                                            (x1, y1, 'o', 'MarkerFaceColor', 'b');
% hold on;
% plot                                            (x1, smoothed_y1, 'b-', 'LineWidth', 2);
% hold on;
% plot                                            (x2, y2, 'o', 'MarkerFaceColor', 'g');
% hold on;
% plot                                            (x2, smoothed_y2, 'g-', 'LineWidth', 2);

% hold off;
% title('1D LOWESS Smoothing');
% legend('Data', 'Smoothed Curve');

% for i = 1:size(ExpressionData_Sig,2)
%     p(i)                                = ranksum(ExpressionData_Sig(:,i), ExpressionData_NonSig(:,i));
% end
% 
% q                                       = mafdr(p);
% selectedGenesIndex                      = find(q<0.05);
% selectedGenesExpression_sig             = ExpressionData_Sig(:,selectedGenesIndex);
% selectedgenes                           = rows_metadata.gene_symbol(selectedGenesIndex); 
% [~,sort_ind]                            = sort(q(q<0.05));
% orderedgenes                            = selectedgenes(sort_ind);
% 
% writecell(orderedgenes,'F:\OneDrive - 北京师范大学\project2\Genes\OrderedGenes1.csv')
% writecell(selectedgenes,'F:\OneDrive - 北京师范大学\project2\Genes\SelectedGenes1.csv')
% writecell(unique(rows_metadata.gene_symbol),'F:\OneDrive - 北京师范大学\project2\Genes\AllGenes.csv')

% removed genes with less than 6 time point in samples of donors
ReExpressionData_Sig                            = ExpressionData_Sig;
ReExpressionData_Sig(:,removeid)                = [];
ReExpressionData_NonSig                         = ExpressionData_NonSig;
ReExpressionData_NonSig(:,removeid)             = [];

ZExpressionData                                 = zscore([ReExpressionData_Sig; ReExpressionData_NonSig]);
ZExpressionData_Sig                             = ZExpressionData(1:61,:);
ZExpressionData_NonSig                          = ZExpressionData(62:end,:);
figure; imagesc(ZExpressionData); colormap(cmap)
saveas                                          (gcf,'~\Genes\figures\expression_matrix.tiff')
print                                           (gcf, '~\Genes\figures\expression_matrix.pdf','-dpdf', ['-r' '600'],'-bestfit')

% Wilcoxon rank-sum test
for i = 1:size(ZExpressionData_Sig,2)
    [p(i),~,stat]                               = ranksum(ZExpressionData_Sig(:,i), ZExpressionData_NonSig(:,i));
    z(i)                                        = stat.zval;
end

[sort_z,sort_ind]                               = sort(z);
orderedallgenes                                 = rows_metadata.gene_symbol(sort_ind);
Weighted                                        = sort_z'; 
orderedtable                                    = table(orderedallgenes,Weighted);
writetable                                      (orderedtable,'~\Genes\enrich\OrderedAllGenes_s6.csv')

q                                               = mafdr(p);
tmpIndex                                        = find(q<0.05);
indices_of_zeros                                = find(removeid == 0);
nth_sig_index                                   = indices_of_zeros(tmpIndex);
selectedgenes                                   = rows_metadata.gene_symbol(nth_sig_index);

[sort_z,sort_ind]                               = sort(z(tmpIndex));
orderedgenes                                    = selectedgenes(sort_ind);
Weighted                                        = sort_z'; 
orderedtable                                    = table(orderedgenes,Weighted);
writetable                                      (orderedtable,'~\Genes\enrich\OrderedSigGenes_s6.csv')

