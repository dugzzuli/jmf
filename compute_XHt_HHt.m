function [sumXHt,sumHHt,XHt_record,HHt_record] = compute_XHt_HHt(X_record,H_record,K)
% compute the product of some matrices used in the main update rule
% <Inputs>:
% X_record :   Cell data with each nonnegative input matrix as its element
% H_record :   Cell data with each coefficient matrix as its element
% K:           Target low rank
% <Outputs>;
% sumXHt:      The sum matrix of element in XHt_record
% sumHHt:      The sum matrix of element in HHt_record
% XHt_record:  Cell data with the product of Xi and Hi as its element
% HHt_record:  Cell data with the product of Hi and Hi' as its element
numN = length(X_record);M = size(X_record{1,1},1);
XHt_record = cell(1,numN); HHt_record = cell(1,numN);
sumXHt = zeros(M,K); sumHHt = zeros(K);
for i = 1:numN
    XHt_record{1,i} = X_record{1,i}*H_record{1,i}'; HHt_record{1,i} = H_record{1,i}*H_record{1,i}'; 
    sumXHt = sumXHt + XHt_record{1,i}; sumHHt = sumHHt + HHt_record{1,i};
end

