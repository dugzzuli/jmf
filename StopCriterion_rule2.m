function retVal=StopCriterion_rule2(X,gradX)
% Stopping criterion 2
pGrad=gradX(gradX<0|X>0);
retVal=norm(pGrad);

