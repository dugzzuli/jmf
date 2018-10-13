function H_record = JMF_prediction_R(X_record,W,theta_record, R_record, L1, L2, r2, K, tol, varargin)
% update by four widely used methods (MUR, PG, Ne, PANLS)
% <Inputs>:
% X_record:      Cell data with each nonnegative input matrix as its element
% W:             Learned basis matrix
% theta_record:  Cell data with network constrained matrix on each data as its element
% R_record:      Cell data with relationship constrained matrix between each pair of data as its element
% L1:            The parameter weigh the must link constraints in theta
% L2:            The parameter weigh the must link constraints in R
% r2:            The parameter limit the growth of Hi
% K;             Target low rank
% tol;           Predefined precision

%        (Below are optional arguments: can be set by providing name-value pairs)
%        ITER_MAX:      The maximum iteration in each subproblem
%        TYPE:          Update rule
%        H_RECORD_INIT: Initial value for low rank coefficients matrices 
        
% <Outputs>:
% H_record:         Cell data with low rank basis matrices
%% Default configuration
iterMax = 100; type = 'PANLS';
% initialize random W, Hi
numN = length(X_record);
vecN = zeros(1,numN);
for i = 1:numN
    vecN(i) = size(X_record{1,i},2);
end
H_record_init = cell(1,numN);
for i = 1:numN
    H_record_init{1,i} = rand(K,vecN(i));
end
% Read optional parameters
if (rem(length(varargin),2) == 1)
    error('Error:Optional parameters should always go by pairs');
else
    for i = 1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'TYPE',                type = varargin{i+1};
            case 'ITER_MAX',            iterMax = varargin{i+1};
            case 'H_RECORD_INIT',       H_record_init = varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end
%%%%%%%%%%%%%%%%%paramter configuration end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test for negative values in X
flag = 0;
for i = 1:numN
    X = X_record{1,i};
    if min(min(X)) < 0
        flag = 1;
    end
    if flag == 1
        error('Error:Input matrix elements can not be negative');
    end
end
%% test for same rows in each X
n1 = size(X_record{1,1},1);
for i = 2:numN
    X = X_record{1,i}; n2 = size(X,1);
    if n1 ~= n2
        error('Error:Input matrices should have the same rows');
    end
end
%% compute the initial error of each part
WtW = W'*W;
H_record = H_record_init;
sumtheta_record = compute_sumtheta_record(H_record,theta_record);
WtX_record = cell(1,numN); sumHRt_record = cell(1,numN);
for i = 1:numN
    WtX_record{1,i} = W'*X_record{1,i};
    sumHRt_record{1,i} = compute_sumHRt(H_record,R_record,i);
end
if ~strcmp(type, 'MUR')
    tolH_record = cell(1,numN); Hsumtheta_record = cell(1,numN);
    for i = 1:numN
        Hsumtheta_record{1,i} = H_record{1,i}*sumtheta_record{1,i};
        [Gradh,~] = get_gradH(WtW,WtX_record{1,i},H_record{1,i},Hsumtheta_record{1,i},sumHRt_record{1,i},L1,L2,r2,K,i);
        delta_init = StopCriterion_rule2(H_record{1,i},Gradh);
        tolH_record{1,i} = max(0.001,tol)*delta_init;
    end
end
%% main update Hi part
for i2 =1:numN
    sumHRt = compute_sumHRt(H_record,R_record,i2);
    switch type
        case 'MUR'
            XitXi = X_record{1,i2}'*X_record{1,i2};
            h = MUR_updateH(XitXi,H_record{1,i2},WtX_record{1,i2},WtW,sumtheta_record{1,i2},sumHRt,L1,L2,r2,iterMax,tol,K);
        case 'PG'
            [h,~,~] = PG_updateH(H_record{1,i2},WtX_record{1,i2},WtW,sumtheta_record{1,i2},Hsumtheta_record{1,i2},sumHRt,L1,L2,r2,iterMax,tolH_record{1,i2});
        case 'Ne'
            L0 = 2*norm(WtW+r2*ones(K));% part Lipschitz constant
            L = L0 + L1*norm(sumtheta_record{1,i2});% Lipschitz constant
            iterMin = 1;
            [h,~,~] = Ne_updateH(H_record{1,i2},WtX_record{1,i2},WtW,sumtheta_record{1,i2},Hsumtheta_record{1,i2},sumHRt,L1,L2,r2,L,iterMin,iterMax,tolH_record{1,i2},K);
        case 'PANLS'
            [h,~] = PANLS_updateH(H_record{1,i2},WtX_record{1,i2},WtW,sumtheta_record{1,i2},Hsumtheta_record{1,i2},sumHRt,L1,L2,r2,iterMax,tolH_record{1,i2});
    end
    H_record{1,i2} = h;
end
