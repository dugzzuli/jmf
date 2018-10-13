function sumHRt = compute_sumHRt(H_record,R_record,index)
% Compute the product of Hi and Ri(index)' and i is not equal to index
% <Inputs>:
% H_record:  Cell data with low rank basis matrices
% R_record:  Cell data with relationship constrained matrix between each pair of data as its element
% index:     The index of the idex-th view data
% <Outputs>:
% sumHRt:    The product between Hi and Ri(index)'
numN = length(H_record); H = H_record{1,index}; [K,N_index] = size(H);
sumHRt = zeros(K,N_index);
for i = 1:numN
    if i ~= index && ~isempty(R_record{index,i})
        sumHRt = sumHRt + H_record{1,i}*R_record{index,i}';
    end
end
  

