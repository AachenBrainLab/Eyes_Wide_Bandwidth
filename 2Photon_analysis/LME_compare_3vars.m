function [pVals, tStats, fullmodel, modelCompare1, modelCompare2, modelCompare3 ] = LME_compare_3vars(targData, predict1, predict2, predict3, randomVar)
% function to compare to role of two variables predict1 and predict2 to explain 
% 'targData' a mixed effect model, while controlling for the impact of the
% random variable 'randomVar'. The pVals and tStats are quantifying how predictive 
% predict1 or predict2 are of targData, when controlling for the role of randomVar in a reduced model.
% SM 09.01.2023


% combine into one table and create model
tbl = table(targData, predict1, predict2, predict3, randomVar, 'VariableNames',{'y','X1','X2','X3', 'randomVar'});
fullmodel = fitlme(tbl,'y ~ X1 + X2 + X3 +(1 |randomVar)');
pVals = fullmodel.Coefficients.pValue(2:4); %return p-Value for both regressors in the full model
tStats = fullmodel.Coefficients.tStat(2:4); %return t-statistics for both regressors in the full model

if nargout > 4
    nullmodel1 = fitlme(tbl,'y ~ X1 + (1|randomVar)');
    nullmodel2 = fitlme(tbl,'y ~ X2 + (1|randomVar)');
    nullmodel3 = fitlme(tbl,'y ~ X3 + (1|randomVar)');
    modelCompare1 = compare(nullmodel1, fullmodel, 'CheckNesting',true);
    modelCompare2 = compare(nullmodel2, fullmodel, 'CheckNesting',true);
    modelCompare3 = compare(nullmodel3, fullmodel, 'CheckNesting',true);
end