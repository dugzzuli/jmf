function W = JMF_prediction_L(X_record, H_record, r1, K, tol, varargin)
% prediction based on left basis matrix W and update by four widely used methods (MUR, PG, Ne, PANLS)
% <Inputs>:
% X_record:      Cell data with each nonnegative input matrix as its element
% H_record:      Cell data with low rank basis matrices
% r1;            The parameter limit the growth of W
% K;             Target low rank
% tol;           Small data control stopping value

%        (Below are optional arguments: can be set by providing name-value pairs)
%        ITER_MAX:   The maximum iteration used in each subproblem
%        TYPE:       Update rule
%        W_INIT:     Initial value for basis matrix
        
% <Outputs>:
% W:         Low rank basis matrix
%% Default configuration
iterMax = 100; type = 'PANLS'; 
% initialize random W,
M = size(X_record{1,1},1); 
numN = length(X_record);
Winit = rand(M,K);  
% Read optional parameters
if (rem(length(varargin),2) == 1)
    error('Error:Optional parameters should always go by pairs');
else
    for i = 1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'TYPE',                type = varargin{i+1};
            case 'ITER_MAX',            iterMax = varargin{i+1};
            case 'W_INIT',              Winit = varargin{i+1};
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
W = Winit;
[sumXHt, sumHHt, ~,~] = compute_XHt_HHt(X_record,H_record,K);
if ~strcmp(type, 'MUR')
    GradW = get_gradW(W,sumXHt,sumHHt,r1);
    delta_init = StopCriterion_rule2(W,GradW);
    tolW = max(0.001,tol)*delta_init;
end
%% main update W part
switch type
    case 'MUR'
        V = [];
        for i = 1:numN
            V = [V,X_record{1,i}];
        end
        VtV = V'*V;
        W = MUR_updateW(VtV,W,sumXHt,sumHHt,r1,iterMax,tol);
    case 'PG'
        [W,~,~] = PG_updateW(W,sumXHt,sumHHt,r1,iterMax,tolW);
    case 'Ne'
        iterMin = 1;
        [W,~,~] = Ne_updateW(W,sumXHt,sumHHt,r1,iterMin,iterMax,tolW);
    case 'PANLS'
        [W,~] = PANLS_updateW(W,sumXHt,sumHHt,r1,iterMax,tolW);
end
