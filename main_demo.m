clear;
%% simulate data 
mu = 2;
[W,H_record,X_record,theta_record,R_record] = JMF_synthetic_dataset1(mu);
data.W = W; data.H_record = H_record; data.X_record = X_record; data.theta_record = theta_record; data.R_record = R_record;
cd('./simulation data/synthetic dataset 1')
save('data1.mat','data')
%% To see the influence of each constraint term on the performance of JMF
load('data1.mat')
% back to the main code file
cd ..
cd ..
X_record = data.X_record; theta_record = data.theta_record; R_record = data.R_record;
K = 4; maxiter = 1000; maxtime = 10^(6); tol = 10^(-7);
% we take the influence of the first part as an example
L1s = [0.001,0.01,0.1,1,10,100,1000]; L2 =0; r1 = 0; r2 = 0;
result = cell(1,length(L1s));
repeat = 10;
for i = 1:length(L1s) 
    L1 = L1s(i);
    re = cell(1,repeat);
    for j = 1:repeat
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run JMF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [W_result, H_result, niter, ~, stop_control] = JMF(X_record, theta_record, R_record, L1, L2, r1, r2, K,maxiter, maxtime, tol,'TYPE','MUR','STOP_RULE','rule 1');
        % can be updated by other update rules
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run JMF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        re{1,j}.W = W_result{1,niter}; re{1,j}.H_record = H_result{1,niter}; re{1,j}.stop_control = stop_control;
    end
    result{1,i} = re;
end    
% To make a new folder to store the results.
cd('./simulation data/synthetic dataset 1')
influence_file = 'influence';
if ~isdir(influence_file)
    mkdir(influence_file);
end
cd('./influence')
save('MUR_result_L1.mat','result')


