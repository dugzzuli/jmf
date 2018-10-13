function [Gradh, GradH] = get_gradH(WtW,WtX_record,H_record,Hsumtheta_record,sumHRt_record,L1,L2,r2,K,index,varargin)
% compute the gradient of Hi or H
% <Inputs>:
% WtW:                W'*W         
% WtX_record:         Cell data with W'*Xi as its element, can be matrix data
% H_record:           Cell data with low rank basis matrices, can be matrix data
% Hsumtheta_record:   Cell data with Hi*sumtheta_i as its element, can be matrix data
% sumHRt_record:      Cell data with sum of Hi*Rj_i' as its element and j is not equal to i, can be matrix data
% L1:                 The parameter weigh the must link constraints in theta
% L2:                 The parameter weigh the must link constraints in R
% r1;                 The parameter limit the growth of W
% r2:                 The parameter limit the growth of Hi
% K;                  Target low rank
%        (Below are optional arguments: can be set by providing name-value pairs)
%        NUM:          Compute gradient of Hi if NUM is 'single', and H if
%        NUM is 'multiple'

% <Outputs>:
% Gradh:  The gradient of Hi
% GradH:  The gradient of H
%% Default configuration
num = 'single';
% Read optional parameters
if (rem(length(varargin),2) == 1)
    error('Error:Optional parameters should always go by pairs');
else
    for i = 1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'NUM',                num = varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end
if strcmp(num,'multiple')
    numN = length(H_record);
    Gradh_record = cell(1,numN); GradH = [];
    for i = 1:numN
        sumHtheta = Hsumtheta_record{1,i};
        Gradh_record{1,i} = -2*WtX_record{1,i} + 2*WtW*H_record{1,i} - L1*sumHtheta - L2*sumHRt_record{1,i} + 2*r2*ones(K)*H_record{1,i};
        GradH = [GradH,Gradh_record{1,i}];
    end
    Gradh = [];
elseif iscell(WtX_record)
    WtX = WtX_record{1,index}; h = H_record{1,index}; Hsumtheta = Hsumtheta_record{1,index}; sumHRt = sumHRt_record{1,index};
    Gradh = -2*WtX + 2*WtW*h - L1*Hsumtheta - L2*sumHRt + 2*r2*ones(K)*h;
    GradH = [];
else
    WtX = WtX_record; h = H_record; Hsumtheta = Hsumtheta_record; sumHRt = sumHRt_record; 
    Gradh = -2*WtX + 2*WtW*h - L1*Hsumtheta - L2*sumHRt + 2*r2*ones(K)*h;
    GradH = [];   
end
    




