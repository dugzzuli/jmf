function [W,iter,Grad] = Ne_updateW(Z,sumXHt,sumHHt,r1,iterMin,iterMax,tol)
% optimize W by solving the subproblem
K = size(Z,2);
l = 2*norm(sumHHt+r1*eye(K));    % Lipschitz constant
W = Z;    % Initialization
Grad = 2*Z*sumHHt-2*sumXHt+2*r1*Z;     % Gradient
alpha1 = 1;
for iter = 1:iterMax,
    W0 = W;
    W = max(Z-Grad/l,eps);    % Calculate sequence 'Y'
    alpha2 = 0.5*(1+sqrt(1+4*alpha1^2));
    Z = W+((alpha1-1)/alpha2)*(W-W0);
    alpha1 = alpha2;
    Grad = 2*Z*sumHHt-2*sumXHt+2*r1*Z; 
    
    % Stopping criteria
    % Lin's stopping condition
    if iter >= iterMin
        pgn = StopCriterion_rule2(Z,Grad);
        if pgn <= tol,
            break;
        end
    end
end
