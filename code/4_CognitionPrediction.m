paras_cognition                                         = readtable('~\paras_cognition.csv');

xVariableNames                                          = {'meanFCS', 'Cp', 'efficiency','modularity','withinFC'}; 
beha_names                                              = {'Cognitive', 'Language', 'Motor'};
nPermutations                                           = 10000;

permutedRvalues                                         = zeros(numel(xVariableNames), numel(beha_names), nPermutations);

for i = 1:numel(xVariableNames)
    xVariableName                                       = xVariableNames(i);
    X                                                   = paras_cognition{:, {xVariableName{1}, 'scan_age', 'sex', 'mFD','deltage','Cogage'}};
    for b = 1:3
        yname                                           = beha_names(b);
        y                                               = paras_cognition{:, yname{1}};
        
        svmModel                                        = fitrsvm(X, y, 'KernelFunction', 'linear', 'Standardize', true);
        cvSVMModel                                      = crossval(svmModel, 'KFold', 10);
        y_pred                                          = kfoldPredict(cvSVMModel);
        
        % actual prediction accuracy        
        [rvalue,~]                                      = corr(y,y_pred);
        rvalues(i,b)                                    = rvalue;
        
        % permutation test
        parfor p = 1:nPermutations
            %             disp(strcat('i=',num2str(i),';b=',num2str(b)))
            y_permuted                                  = y(randperm(length(y)));
            
            [rvalue_perm,pvalue_perm]                   = corr(y_permuted,y_pred);
            permutedRvalues(i, b, p)                    = rvalue_perm;
        end
%         parfor p = 1:nPermutations
%             y_permuted                                  = y(randperm(length(y)));
%             svmModelPermuted                            = fitrsvm(X, y_permuted, 'KernelFunction', 'linear', 'Standardize', true);
%             cvSVMModelPermuted                          = crossval(svmModelPermuted, 'KFold', 10);
%             y_pred_permuted                             = kfoldPredict(cvSVMModelPermuted);
%             
%             [rvalue_perm,pvalue_perm]                   = corr(y_permuted,y_pred_permuted);
%             SS_res_perm                                 = sum((y_permuted - y_pred_permuted).^2);
%             SS_tot_perm                                 = sum((y_permuted - mean(y_permuted)).^2);
%             R_squared_perm                              = 1 - (SS_res_perm/SS_tot_perm);
%             
%             permutedR2(i, b, p)                         = R_squared_perm;
%             permutedRvalues(i, b, p)                    = rvalue_perm;
%         end
        
        % plot scatters
        lm                                              = fitlm(y,y_pred);
        y_fit                                           = predict(lm, y);
        
        f                                               = figure;
        scatter                                         (y, y_pred,'filled');
        hold on;
        plot                                            (y, y_fit, 'r', 'LineWidth', 2);
        
        title                                           (strcat('Predict'," ",yname{1}," ",'using'," ",xVariableName{1}));
        xlabel                                          ('Actual values');
        ylabel                                          ('Predicted values');
        hold off;
        
        saveas                                          (f,strcat('~\cognition\GlobalVarsPermutation\',xVariableName{1},'_',yname{1},'_scatter.tiff'))
        print                                           (strcat('~\cognition\GlobalVarsPermutation\',xVariableName{1},'_',yname{1},'_scatter.pdf'), '-dpdf', ['-r' '600'],'-bestfit')
        
    end
end

% compute pvalue using permutation tests
pValue_r                                                = zeros(numel(xVariableNames), numel(beha_names));

for i = 1:numel(xVariableNames)
    xVariableName                                       = xVariableNames(i);
    for b = 1:3
        yname                                           = beha_names(b);
        
        observed_r                                      = rvalues(i, b);
        permuted_r                                      = squeeze(permutedRvalues(i, b, :));
        pValue_r(i, b)                                  = mean(permuted_r > observed_r);
        
        % plot hist
        f                                               = figure;
        histfit                                         (permuted_r);
        line                                            ([observed_r observed_r],[0 300],'linewidth',1.5);
%         ylim                                        ([0 90]);
        saveas                                          (f,strcat('~\cognition\GlobalVarsPermutation\',xVariableName{1},'_',yname{1},'_hist.tiff'))
        print                                           (strcat('~\cognition\GlobalVarsPermutation\',xVariableName{1},'_',yname{1},'_hist.pdf'), '-dpdf', ['-r' '600'],'-bestfit')
        
    end
end
save                                                    ('~\cognition\pValue_r_GlobalVars.mat','pValue_r')
save                                                    ('~\cognition\r_GlobalVars.mat','rvalues')

fprintf('Correlation coefficient r p-values:\n');
disp(pValue_r);
disp(rvalues)

% plot r and p matrix
figure;
imagesc(rvalues);
colorbar; 
title('Observed r and p-values');
colormap(cmap); 
caxis([-0.15, 0.15]);

% text p label
[numRows, numCols] = size(rvalues);
for row = 1:numRows
    for col = 1:numCols
        pValText                                        = sprintf('p=%.3f', pValue_r(row, col));
        text                                            (col, row, pValText, 'HorizontalAlignment', 'center', ...
                                                            'VerticalAlignment', 'middle', 'Color', 'k', 'FontSize', 8);
    end
end

set                                                     (gca, 'XTick', 1:numel(beha_names), 'XTickLabel', beha_names);
set                                                     (gca, 'YTick', 1:numel(xVariableNames), 'YTickLabel', xVariableNames);

saveas                                                  (gcf,'~\cognition\rvalues.tiff')
print                                                   ('~\cognition\rvalues.pdf', '-dpdf', ['-r' '600'],'-bestfit')
