function fval = gridSearch(mTrain, pTrain, vals)

model = svmtrain(mTrain', pTrain, ['-c ' num2str(vals(1)) ' -g ' num2str(vals(2)) ' -b 1']);

[predict_train, accuracyTrain] = ...
    svmpredict(mTrain', pTrain, model, '-b 1');

predict_train = predict_train';
true_pos = length(mTrain(mTrain == 1 & predict_train == 1));
false_pos = length(mTrain(mTrain == -1 & predict_train == 1));
true_neg = length(mTrain(mTrain == -1 & predict_train == -1));
false_neg = length(mTrain(mTrain == 1 & predict_train == -1));

sensitivity = true_pos / (true_pos + false_neg)
specificity = true_neg / (true_neg + false_pos)
ppv = true_pos / length(predict_train(predict_train == 1))

% fval = -(sensitivity^2 + specificity^2 -(sensitivity - specificity)^2);
fval = (specificity - sensitivity)^2;