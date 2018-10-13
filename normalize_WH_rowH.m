function [W,H] = normalize_WH_rowH(W0,H0,value)
% normalization after each step
% <Inputs>:
% W0:    Low rank basis matrix after one step
% H0:    Stacked low rank coefficients matrices after one step
% value: The sum of each row of H0 and often set to 1
% <Outputs>:
% W:    Low rank normalized basis matrix
% H:    Stacked low rank normalized coefficients matrices

sh = sum(H0,2);  sh(sh == 0) = value; psh = value./sh; 
W=W0*diag(sh);H=diag(psh)*H0;
    