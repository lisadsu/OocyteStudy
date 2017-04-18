% find best RBF sigma to minimize cross-validation error


function bestSigma = findBestSVMParams(params, m)


initClassifier = fitcsvm(params, m, 'KernelFunction', 'rbf', ...
    'KernelScale', 'auto');
startSigma = initClassifier.KernelParameters.Scale;

opts = optimset('display', 'iter', 'TolFun', 1e-4, 'TolX', .001, ...
    'MaxFunEvals', 100, 'Algorithm', 'sqp');
[bestSigma, fval, exitflag] = fmincon(@(sigma) evalAUC(params,m,sigma), ...
    startSigma, [], [], [], [], 0.001, 2, [], opts)













