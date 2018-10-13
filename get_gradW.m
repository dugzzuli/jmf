function GradW = get_gradW(W,sumXHt,sumHHt,r1)
% compute the gradient of W
GradW = 2*(W*sumHHt - sumXHt + r1*W); 

