function [W,grad,iter] = PG_updateW(W,sumXHt,sumHHt,r1,iterMax,tol)
% optimize W by solving the subproblem
alpha = 1; beta = 0.1; [~,K] = size(W); 
for iter = 1:iterMax
    grad = get_gradW(W,sumXHt,sumHHt,r1);
    projgrad = norm(grad(grad < 0 | W >0));
    if projgrad < tol
        break
    end
    
    % search step size
    for inner_iter = 1:20
        Wn = max(W-alpha*grad,eps);d = Wn-W;
        Q = 2*sumHHt + 2*r1*eye(K); 
        gradd = sum(sum(grad.*d)); 
        dQd = sum(sum((d*Q).*d));
        suff_decr = 0.99*gradd + 0.5*dQd < 0;
        if inner_iter == 1,
            decr_alpha = ~suff_decr; Wp = W;
        end
        if decr_alpha,
            if suff_decr,
                W = Wn; break;
            else
                alpha = alpha * beta;
            end
        else
            if ~suff_decr || isequal(Wp, Wn),
                W = Wp; break;
            else
                alpha = alpha/beta; Wp = Wn;
            end
        end
    end
end
