% given predicted survival values and distances from classifier decision
% boundary, 

function survivalLikelihoods = findSurvivalLikelihoods(decDistTest, ...
    decDistTrain, mOut, methodToUse)

mTrain = mOut(mOut < 2); %(decDistTrain > 0);
mTest = (decDistTest > 0);

% decDistTest corresponds to mTest
% decDistTrain corresponds to mTrain

if methodToUse == 1
    
    % ========================================================================
    % Method 1: Estimate based on how many viable and nonviable embryos have a
    % similar distance to the decision boundary
    % ========================================================================
    
    % for a given embryo in the test set, calculate its likelihood of survival
    % by looking at the relative numbers of viable/non-viable embryos in the
    % training set with decision boundary distances within +/- decInterval
    
    % weighting envelope is 1/10th of total decDistTrain range
    envelopeWidth = (max(decDistTrain) - min(decDistTrain))/10;
    survivalLikelihoods = zeros(1,length(decDistTest));
    
    for i = 1:length(mTest)
        
        % define gaussian envelope centered at curr decDist
        % then compute locally weighted viability average
        currDist = decDistTest(i);
        gaussCurrDist = normpdf(decDistTrain, currDist, envelopeWidth);
        survivalLikelihoods(i) = sum(gaussCurrDist .* mTrain') / ...
            sum(gaussCurrDist);
        
    end
    
elseif methodToUse == 2
    
    % ========================================================================
    % Method 2: Estimate based on Platt's method, fit to logistic curve
    % ========================================================================    
    
        
%         minFun = inline(['sum(mTest.*(AB(1)*decDistTest + AB(2)) ' ...
%         '+ log(1 + exp(-AB(1)*decDistTest - AB(2))))'], ...
%         'mTest', 'decDistTest', 'AB');
    

    options = optimset('display', 'off');
    ABout = fminsearch(@(AB) findSigmoidParams(AB, mTrain', decDistTrain), [0 0], options);
    
    survivalLikelihoods = ((1 + exp(ABout(1)*decDistTest + ABout(2))).^(-1))';
    
end
