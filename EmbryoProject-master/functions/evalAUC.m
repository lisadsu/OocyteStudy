

function AUC = evalAUC(params, m, sigma)

% want to maximize the AUC_ROC, NOT minimize the k-fold loss
% maybe modify to maximize AUC_PR instead?

currClassifier = fitcsvm(params, m, 'KernelFunction', 'rbf', ...
    'KernelScale', sigma, 'CrossVal', 'on');
[~, decDist] = kfoldPredict(currClassifier);

[~, ~, ~, AUC] = perfcurve(m, decDist(:,2), 1);
AUC = -1*AUC;



