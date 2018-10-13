function W = MUR_updateW(VtV,W,sumXHt,sumHHt,r1,iterMax,tol)
% solve the subproblem: argminF(W) =||X1-WH1||_F^2+...+||XN-WHN||_F^2+r1||W||_F^2 
WtW = W'*W;
WtXHt = W'*sumXHt;
WtWHHt = WtW*sumHHt;
delta_init = trace(VtV) - 2*trace(WtXHt) + trace(WtWHHt) + r1*norm(W,'fro')^2 ;
for iter = 1 : iterMax
    sumWHHt = W*sumHHt;
    W = W.* sumXHt./(sumWHHt+r1*W+eps); %update W
    % compute F(W)
    WtW = W'*W;
    WtXHt = W'*sumXHt;
    WtWHHt = WtW*sumHHt;
    delta = trace(VtV) - 2*trace(WtXHt) + trace(WtWHHt) + r1*norm(W,'fro')^2;
    if iter == 1
        delta_old = delta_init;
    end
    stop_control = (delta_old - delta)/(delta_init - delta);
    %     stop_control2(iter) = delta2/delta2_init;
    if abs(stop_control) < tol && iter >10
        break,
    end
    delta_old  = delta;
end
    
    