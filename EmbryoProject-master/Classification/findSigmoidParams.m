function fval = findSigmoidParams(AB, mTest, decDistTest)

Pab = (1 + exp(AB(1)*decDistTest + AB(2))) .^ (-1);

fval = -1*sum(mTest.*log(Pab) + (1-mTest).*log(1-Pab));