

function AUC = evalKfoldLoss(params, m, sigma)

% want to maximize the AUC_ROC, NOT minimize the k-fold loss

currClassifier = fitcsvm(params, m, 'KernelFunction', 'rbf', ...
    'KernelScale', sigma);
crossValModel = crossval(currClassifier, 'Kfold', 10);

[~, ~, ~, AUC] = perfcurve(m, -decDistTest, 1);
AUC = kfoldLoss(crossValModel);



