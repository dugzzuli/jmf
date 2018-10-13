function h = MUR_updateH(XitXi,h,WtX,WtW,sumtheta,sumHRt,L1,L2,r2,iterMax,tol,K)
% optimize Hi by MUR
HHt = h*h';
WtXHt = WtX*h'; HHtWtW = HHt*WtW;
Hsumtheta = h*sumtheta; HsumthetaHt = Hsumtheta*h'; sumHRtHt = sumHRt*h';
oHHto = ones(1,K)*HHt*ones(K,1);
delta_init = trace(XitXi)- 2*trace(WtXHt) + trace(HHtWtW) - L1*trace(HsumthetaHt) - L2*trace(sumHRtHt) + r2*trace(oHHto);
for iter = 1 : iterMax
    h = h.*(WtX + L1/2*Hsumtheta + L2/2*sumHRt)./((WtW+r2*ones(K))*h+eps);
    % compute stop value
    HHt = h*h';
    WtXHt = WtX*h'; HHtWtW = HHt*WtW;
    Hsumtheta = h*sumtheta; HsumthetaHt = Hsumtheta*h'; sumHRtHt = sumHRt*h';
    oHHto = ones(1,K)*HHt*ones(K,1);
    delta = trace(XitXi)- 2*trace(WtXHt) + trace(HHtWtW) - L1*trace(HsumthetaHt) - L2*trace(sumHRtHt) + r2*trace(oHHto);
    if iter == 1
        delta_old = delta_init;
    end
    stop_control = abs((delta_old - delta)/(delta_init - delta));
    if stop_control < tol && iter >10
        break,
    end
end
