function sumtheta_record = compute_sumtheta_record(H_record,theta_record)
% Compute the the sum of constraint network matrices on each view data
% <Inputs>
% H_record:        Cell data with low rank basis matrices
% theta_record:    Cell data with network constrained matrix on each data as its element
% <Outputs>
% sumtheta_record: Cell data with the sum of constraint network matrices on
% each view data as its element
numN = length(H_record); sumtheta_record = cell(1,numN);
for index = 1:numN
    H = H_record{1,index}; [~,N] = size(H); 
    theta_element = theta_record{1,index};
    %compute sumtheta
    if isempty(theta_element)
       sumtheta = zeros(N);
    elseif iscell(theta_element)
        num_theta = length(theta_element); sumtheta = zeros(N);
        for j1 = 1:num_theta
            sumtheta = sumtheta + theta_element{1,j1}+theta_element{1,j1}';
        end
    else
        sumtheta = theta_element+theta_element';
    end
    sumtheta_record{1,index} = sumtheta;
end