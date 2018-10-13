function [H,iterH] = PANLS_updateH(Hinit,WtX,WtW,sumtheta,Hsumtheta,sumHRt,L1,L2,r2,iterMax,tol)
% <INPUT>
%        H, grad: output solution and gradient
%        iter: #iterations used
%        V, W: constant matrices
%        Hinit: initial solution
%        tol: stopping tolerance
%        maxiter: limit of iterations
% <OUTPUT>
%       H
%       iterW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global  alpha beta eta rho eps0 n_1 n_2  r
alpha = 1; beta = 0.1; eta = 0.1; rho = 0.5; n_1 = 2; n_2 = 1;
eps0 = 1; r = 2; tau = 1e-7;
[K,M] = size(Hinit);
H = Hinit;
Q1 = 2*(WtW+r2*ones(K)+tau*eye(K)); Q2= -L1*sumtheta;
[grad, ~] = get_gradH(WtW,WtX,H,Hsumtheta,sumHRt,L1,L2,r2,K,1);
grad = grad + tau*(H-Hinit);
projnorm  = norm(grad(grad < 0 | H >0));
iterH = 0;
while  iterH <  iterMax
    
    if projnorm <  tol
        break
    end
    %%% PG loop
    inner_iter_PG = 0;
    while projnorm >= tol && iterH < iterMax
        
        inner_iter_PG = inner_iter_PG + 1;
        
        iterH = iterH+1;
        for inner_iter = 1:20,
            
            Hn = max(H - alpha*grad, 0); d = Hn-H;
            gradd = sum(sum(grad.*d)); dQd = sum(sum((Q1*d).*d))+sum(sum((d*Q2).*d));
            suff_decr = 0.99*gradd + 0.5*dQd < 0;
            if inner_iter==1,
                decr_alpha = ~suff_decr; Hp = H;
            end
            if decr_alpha,
                if suff_decr,
                    H = Hn; break;
                else
                    alpha = alpha * beta;
                end
            else
                if ~suff_decr || isequal(Hp, Hn)
                    H = Hp; break;
                else
                    alpha = alpha/beta; Hp = Hn;
                end
            end
        end
        
        [grad, ~] = get_gradH(WtW,WtX,H,Hsumtheta,sumHRt,L1,L2,r2,K,1); grad = grad + tau*(H-Hinit);
        projnorm = norm(grad(grad < 0 | H >0));
        
        IA = find(H > 0);
        Z = H(IA);
        g_red = grad(IA);
        ng_I = norm(g_red);
        
        if ng_I < eta*projnorm
            eta=rho*eta;
        else
            if inner_iter_PG >=  n_1
                break;
            end
        end
    end
% 	if isempty(IA),break,end
    %%%
    %%% CG loop
    d_CG = -g_red;

    while projnorm >= tol && iterH < iterMax

        iterH = iterH + 1;

        d_ex = zeros(K,M);
        
        %  for inner_iter = 1:5
        for inner_iter = 1:5
            Zn = max(Z + alpha*d_CG, 0); s0 = Zn-Z;
            gradd = g_red'*s0;
            % d_ex = zeros(nr,nc);
            d_ex(IA) = s0;
            y0_ex = Q1*d_ex + d_ex*Q2;
            y0 = y0_ex(IA);
            dQd = y0'*s0;
            suff_decr = 0.9*gradd + 0.5*dQd < 0;
            if inner_iter == 1,
                decr_alpha = ~suff_decr; Zp = Z;
            end
            if decr_alpha,
                if suff_decr,
                    Z = Zn; break;
                else
                    alpha = alpha * beta;
                end
            else
                if ~suff_decr || isequal(Zp,Zn)
                    Z = Zp; break;
                else
                    alpha = alpha/beta;
                    Zp = Zn;
                end
            end
        end
        H(IA) = Z;
        [grad, ~] = get_gradH(WtW,WtX,H,Hsumtheta,sumHRt,L1,L2,r2,K,1); grad = grad + tau*(H-Hinit);
        %%
        g_red = grad(IA);
        %%
        projnorm = norm(g_red(g_red < 0 | Z >0));
        
        ng_I = norm(g_red);
        
        if ng_I < eta * projnorm
            break
        else
            a1 = length(find(Z == 0));
            
            if a1 > n_2

               IA = find(H>0);
               g_red = grad(IA);
               d_CG = -g_red;
               Z = H(IA);
            else
                if a1 <= n_2 && a1 > 0 && ~ isempty( find( abs(g_red) >= projnorm^(1/2) & Z >= projnorm^(3/2) ,1))
                    break
                else
                    d0 = d_CG;
                    t1 = eps0 * ng_I^r ;
                    z0 = y0 + t1 * s0;
                    d0z0 = d0'*z0;
                    g1d0 = g_red'*d0;
                    if d0z0 < 1.0e-8
                    beta_1 = 0;
                    theta_1 = 0;
                    else
                    beta_1 = (g_red'*z0)/d0z0 - 2 * (z0'*z0 * g1d0)/(d0z0^2);
                    theta_1 = g1d0/d0z0;
                    end
                    d_CG = -g_red + beta_1 * d0 + theta_1 * z0;  
                end
            end
        end
    end
end

