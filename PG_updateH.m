function [H,grad,iter] = PG_updateH(H,WtX,WtW,sumtheta,Hsumtheta,sumHRt,L1,L2,r2,iterMax,tol)
% optimize Hi by solving the subproblem
beta = 0.1; alpha = 1;
[K,~] = size(H);
for iter = 1:iterMax
    [grad, ~] = get_gradH(WtW,WtX,H,Hsumtheta,sumHRt,L1,L2,r2,K,1);
    projgrad = norm(grad(grad < 0 | H >0));
    if projgrad < tol
        break
    end
    
    % search step size
    for inner_iter = 1:20
        Hn = max(H - alpha*grad, eps); d = Hn-H;
        %Q = 2*kron(eye(M),WtW) - L1*kron(sumtheta_record{1,index},eye(K)) + 2*r2*kron(eye(M),ones(K)); 
         Q1 = 2*WtW+2*r2*ones(K); Q2 = -L1*sumtheta;
        gradd = sum(sum(grad.*d)); 
        %dQd = (vect(d))'*Q*vect(d);
        dQd = sum(sum((Q1*d).*d)) + sum(sum((d*Q2).*d));
        suff_decr = 0.99*gradd + 0.5*dQd < 0;
        if inner_iter == 1,
            decr_alpha = ~suff_decr; Hp = Hn;
        end
        if decr_alpha
            if suff_decr
                H = Hn; break;
            else
                alpha = alpha * beta;
            end
        else
            if ~suff_decr || isequal(Hp,Hn),
                H = Hp; break;
            else
                alpha = alpha/beta; Hp = Hn;
            end
        end
    end
end