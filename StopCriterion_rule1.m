function [retVal,diff1]=StopCriterion_rule1(VtV,WtW,W,sumXHt,sumHHt,diff2,diff3,diff4,diff5)
% Compute the first topping criterion
WtXHt = W'*sumXHt;
WtWHHt = WtW*sumHHt;
diff1 = trace(VtV) - 2*trace(WtXHt) + trace(WtWHHt);
retVal = diff1+diff2+diff3+diff4+diff5;
        
  