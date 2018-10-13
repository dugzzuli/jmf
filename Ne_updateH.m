function [H,iter,Grad]=Ne_updateH(Z,WtX,WtW,sumtheta,Hsumtheta,sumHRt,L1,L2,r2,L,iterMin,iterMax,tol,K)
% optimize Hi by solving the subproblem
H=Z;    % Initialization
Grad = -2*WtX + 2*WtW*H - L1*Hsumtheta - L2*sumHRt + 2*r2*ones(K)*H;  % Gradient
alpha1=1;
for iter=1:iterMax,
    H0=H;
    H=max(Z-Grad/L,eps);    % Calculate sequence 'Y'
    alpha2=0.5*(1+sqrt(1+4*alpha1^2));
    Z=H+((alpha1-1)/alpha2)*(H-H0);
    alpha1=alpha2;
    Grad = -2*WtX + 2*WtW*Z - L1*Z*sumtheta - L2*sumHRt + 2*r2*ones(K)*Z;  % Gradient;
    
    % Stopping criteria
    % Lin's stopping condition
    if iter >= iterMin
        pgn=StopCriterion_rule2(Z,Grad);
        if pgn<=tol,
            break;
        end
    end
end
