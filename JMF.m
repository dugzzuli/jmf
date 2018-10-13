function [W_result, H_result, niter, elapse_record, stop_control] = JMF(X_record, theta_record, R_record, L1, L2, r1, r2, K,maxiter, maxtime, tol, varargin)
% update by four widely used methods (MUR, PG, Ne, PANLS)
% <Inputs>:
% X_record:      Cell data with each nonnegative input matrix as its element
% theta_record:  Cell data with network constrained matrix on each data as its element
% R_record:      Cell data with relationship constrained matrix between each pair of data as its element
% L1:            The parameter weigh the must link constraints in theta
% L2:            The parameter weigh the must link constraints in R
% r1;            The parameter limit the growth of W
% r2:            The parameter limit the growth of Hi
% K;             Target low rank
% maxiter:       The maximum number of iterations to run
% maxtime:       The maxmuim time limition
% tol;           Small data control stopping value

%        (Below are optional arguments: can be set by providing name-value pairs)
%        ITER_MAX:      The maximum iteration in each subproblem
%        TYPE:          Update rule
%        STOP_RULE:     The stop criterion
%        W_INIT:        Initial value for basis matrix
%        H_RECORD_INIT: Initial value for low rank coefficients matrices 
        
% <Outputs>:
% W_result:         Cell data with low rank basis matrix on each iteration
% H_result:         Cell data with low rank coefficients matrices on each iteration
% niter:            The number of iterations
% elapse_record:    Record cpu time of each iteration
% stop_control:     Record stop value by rule1 on each iteration
%% Default configuration
iterMax = 100; type = 'PANLS'; stop_rule = 'rule 2';
% initialize random W, Hi
M = size(X_record{1,1},1); numN = length(X_record); vecN = zeros(1,numN);
for i = 1:numN
    vecN(i) = size(X_record{1,i},2);
end
Winit = rand(M,K); H_record_init = cell(1,numN);
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
            case 'STOP_RULE',           stop_rule = varargin{i+1};
            case 'ITER_MAX',            iterMax = varargin{i+1};
            case 'W_INIT',              Winit = varargin{i+1};
            case 'H_RECORD_INIT',       H_record_init = varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end
%%%%%%%%%%%%%%%%%paramter configuration end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MUR only stopped by rule 1
if strcmp(type,'MUR') && strcmp(stop_rule, 'rule 2')
    error('Error:MUR is stopped by rule 1')
end
%% record the columns of each data matrix
vec_sumN = zeros(1,3); vec_sumN(1) = vecN(1);
for j = 2:numN
    sumn = 0;
    for s=1:j
        sumn = sumn+vecN(s);
    end
    vec_sumN(j) = sumn;
end
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
%% normalization
H_init = H_record_init{1,1};
for i = 2:numN
    H_init = [H_init, H_record_init{1,i}];
end
value = 1;
[W,H] = normalize_WH_rowH(Winit,H_init,value);
%% compute the initial error of each part
WtW = W'*W;
H_record = cell(1,numN); H_record{1,1} = H(:,1:vec_sumN(1));
for i = 2:numN
    H_record{1,i} = H(:,vec_sumN(i-1)+1:vec_sumN(i));
end
sumtheta_record = compute_sumtheta_record(H_record,theta_record);
Hsumtheta_record = cell(1,numN); sumHRt_record = cell(1,numN);
for i = 1:numN
    Hsumtheta_record{1,i} = H_record{1,i}*sumtheta_record{1,i};
    sumHRt_record{1,i} = compute_sumHRt(H_record,R_record,i);
end
[sumXHt, sumHHt, ~,HHt_record] = compute_XHt_HHt(X_record,H_record,K);
if strcmp(stop_rule,'rule 1')
    V = [];
    for i = 1:numN
        V = [V, X_record{1,i}];
    end
    VtV = V'*V;
    [diff2,diff3,diff4,diff5] = compute_diff(WtW,H_record,HHt_record,Hsumtheta_record,sumHRt_record,L1,L2,r1,r2,K);
    [delta1_init,~] = StopCriterion_rule1(VtV,WtW,W,sumXHt,sumHHt,diff2,diff3,diff4,diff5);
end
if ~strcmp(type, 'MUR')
    WtX_record = cell(1,numN);
    for i = 1:numN
        WtX_record{1,i} = W'*X_record{1,i};
    end
    GradW = get_gradW(W,sumXHt,sumHHt,r1);
    [~, GradH] = get_gradH(WtW,WtX_record,H_record,Hsumtheta_record,sumHRt_record,L1,L2,r2,K,1,'NUM','multiple');
    delta_init = StopCriterion_rule2([W',H],[GradW',GradH]);
    tolW = max(0.001,tol)*delta_init; tolH = repmat(tolW,1,numN);
end
elapse_record = zeros(1,maxiter); W_result = cell(1,maxiter); H_result = cell(1,maxiter); stop_control = zeros(1,maxiter);
%% main part
for iter = 1:maxiter
    elapse = cputime;
    % update W with any Hi fixed
    switch type
        case 'MUR'
            sumWHHt = W*sumHHt;
            W0 = W.* sumXHt./(sumWHHt+r1*W+eps);
            WtW = W0'*W0;
        case 'PG'
            [W0,~,iterW] = PG_updateW(W,sumXHt,sumHHt,r1,iterMax,tolW);
            WtW = W0'*W0;
            if iterW == 1,
                tolW = 0.1 * tolW;
            end
        case 'Ne'
            iterMin = 1;
            [W0,iterW,~] = Ne_updateW(W,sumXHt,sumHHt,r1,iterMin,iterMax,tolW);
            WtW = W0'*W0;
            L0 = 2*norm(WtW+r2*ones(K));% part Lipschitz constant
            if iterW == 1,
                tolW = 0.1 * tolW;
            end
        case 'PANLS'
            [W0,iterW] = PANLS_updateW(W,sumXHt,sumHHt,r1,iterMax,tolW);
            WtW = W0'*W0;
            if iterW == 1,
                tolW = 0.1 * tolW;
            end
    end    
    % update Hi with W fixed
    WtX_record = cell(1,numN);
    for i1 = 1:numN
        WtX_record{1,i1} = W0'*X_record{1,i1};
    end
    
    for i2 =1:numN
        WtX = WtX_record{1,i2};
        Hsumtheta = Hsumtheta_record{1,i2};
        sumHRt = compute_sumHRt(H_record,R_record,i2);
        switch type
            case 'MUR'
                h = H_record{1,i2}.*(WtX + L1/2*Hsumtheta + L2/2*sumHRt)./((WtW+r2*ones(K))*H_record{1,i2}+eps);
            case 'PG'
                [h,~,iterh] = PG_updateH(H_record{1,i2},WtX,WtW,sumtheta_record{1,i2},Hsumtheta,sumHRt,L1,L2,r2,iterMax,tolH(i2));
                if iterh == 1,
                    tolH(i2) = 0.1 * tolH(i2);
                end
                
            case 'Ne'
                L = L0 + L1*norm(sumtheta_record{1,i2});% Lipschitz constant
                iterMin = 1;
                [h,iterh,~]=Ne_updateH(H_record{1,i2},WtX,WtW,sumtheta_record{1,i2},Hsumtheta,sumHRt,L1,L2,r2,L,iterMin,iterMax,tolH(i2),K);
                if iterh == 1,
                    tolH(i2) = 0.1 * tolH(i2);
                end
                
            case 'PANLS'
                [h,iterh] = PANLS_updateH(H_record{1,i2},WtX,WtW,sumtheta_record{1,i2},Hsumtheta,sumHRt,L1,L2,r2,iterMax,tolH(i2));
                if iterh == 1,
                    tolH(i2) = 0.1 * tolH(i2);
                end
        end
        H_record{1,i2} = h;
    end
    % normalization
    H0 = [];
    for j = 1:numN
        H0 = [H0, H_record{1,j}];
    end
    [W,H] = normalize_WH_rowH(W0,H0,value);
    H_record = cell(1,numN); H_record{1,1} = H(:,1:vec_sumN(1));
    for i3 = 2:numN
        H_record{1,i3} = H(:,vec_sumN(i3-1)+1:vec_sumN(i3));
    end
    W_result{1,iter} = W; H_result{1,iter} = H_record;
    % stop criterion
    WtW = W'*W;
    Hsumtheta_record = cell(1,numN); sumHRt_record = cell(1,numN);
    for i4 = 1:numN
        Hsumtheta_record{1,i4} = H_record{1,i4}*sumtheta_record{1,i4};
        sumHRt_record{1,i4} = compute_sumHRt(H_record,R_record,i4);
    end
    [sumXHt, sumHHt, ~,HHt_record] = compute_XHt_HHt(X_record,H_record,K);
    switch stop_rule
        case 'rule 1'
            [diff2,diff3,diff4,diff5] = compute_diff(WtW,H_record,HHt_record,Hsumtheta_record,sumHRt_record,L1,L2,r1,r2,K);
            [delta,~] = StopCriterion_rule1(VtV,WtW,W,sumXHt,sumHHt,diff2,diff3,diff4,diff5);
            if iter == 1
                delta_old = delta1_init;
            end
            stop_control(iter) = abs((delta_old - delta)/(delta1_init - delta));
            delta_old  = delta;          
        case 'rule 2'
            GradW = get_gradW(W,sumXHt,sumHHt,r1);
            [~, GradH] = get_gradH(WtW,WtX_record,H_record,Hsumtheta_record,sumHRt_record,L1,L2,r2,K,1,'NUM','multiple');
            delta = StopCriterion_rule2([W',H],[GradW',GradH]);
            stop_control(iter) = delta/delta_init;
            if iter > 100
                Dstop = stop_control(iter) - stop_control(iter-10);
            else
                Dstop = 1;
            end
            if abs(Dstop) < 10^(-3)*tol, break, end
    end
    if stop_control(iter) < tol && iter > 10
        break,
    end
    elapse_record(1,iter) = cputime-elapse;
    if sum(elapse_record) > maxtime,break,end
    % return the inappropriate parameters information
    if ~isempty(find(isnan(H)))
        disp('The parameters is inappropriate')
        break,
    end
end
niter = iter;
elapse_record = elapse_record(1:niter); stop_control = stop_control(1:niter); W_result = W_result(1:niter); H_result = H_result(1:niter);
