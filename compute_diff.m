function [diff2,diff3,diff4,diff5] = compute_diff(WtW,H_record,HHt_record,Hsumtheta_record,sumHRt_record,L1,L2,r1,r2,K)
% compute error of each part other than the first part
% <Inputs>:
% WtW:              W'*W
% H_record:         Cell data with low rank basis matrices
% HHt_record:       Cell data with Hi*Hi'as its element
% Hsumtheta_record: Cell data with Hi*sumtheta_i as its element
% sumHRt_record:    Cell data with sum of Hi*Rj_i' as its element and j is
% not equal to i
% L1:               The parameter weigh the must link constraints in theta
% L2:               The parameter weigh the must link constraints in R
% r1;               The parameter limit the growth of W
% r2:               The parameter limit the growth of Hi
% K:                Target low rank
% <Outputs>:
% diff2,diff3,diff4,diff5:  The value of the part of objective function
% other than the first part
numN = length(H_record);
D2 = 0;
for i = 1:numN
    Hsumtheta = Hsumtheta_record{1,i};
    HtheHt = Hsumtheta*H_record{1,i}';
    D2 = D2 + trace(HtheHt);
end
diff2 = -L1*D2;
D3 = 0;
for i = 1:numN
    sumHRt = sumHRt_record{1,i};
    HRHt = H_record{1,i}*sumHRt';
    D3 = D3 + trace(HRHt);
end
diff3 = -L2/2*D3;
diff4 = r1*trace(WtW);
D5 = 0;
for i = 1:numN
    HHt = HHt_record{1,i}';
    D5 = D5+ ones(1,K)*HHt*ones(K,1);
end
diff5 = r2*D5;