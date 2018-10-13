function [W,H_record,X_record,theta_record,R_record] = JMF_synthetic_dataset4(mu)
% simulation synthetic dataset 4 
% <Input>:
% mu:    The noise level
% <Outputs>:
% W:             The golden standard of low rank basis matrix
% H_record:      The golden standard of low rank coefficients matrices
% X_record:      The golden standard of original data matrices
% theta_record:  Cell data with network constrained matrix on each data as its element
% R_record:      Cell data with relationship constrained matrix between each pair of data as its element
w = ones(1,100); h1 = ones(1,240); h2 = ones(1,260); h3 = ones(1,400);
W = zeros(500,5);    
H1 = zeros(5,1200); H2 = zeros(5,1300); H3 = zeros(5,2000); 
for i = 1:5
    W(100*(i-1)+1:100*i,i) = w';
end
W0 = W;
for j = 1:5
    H1(j,240*(j-1)+1:240*j) = h1;
    H2(j,260*(j-1)+1:260*j) = h2;
    H3(j,400*(j-1)+1:400*j) = h3;
end    
H0_record = cell(1,3); H0_record{1,1} = H1; H0_record{1,2} = H2; H0_record{1,3} = H3;
% Original matrix
X1 = W*H1; X2 = W*H2; X3 = W*H3;

% add noise
XX1 = X1 + mu*randn(500,1200);  RX1 = XX1;
XX2 = X2 + mu*randn(500,1300);  RX2 = XX2;
XX3 = X3 + mu*randn(500,2000);  RX3 = XX3;

I1 = find(XX1<0);  XX1(I1) = 0;
I2 = find(XX2<0);  XX2(I2) = 0;
I3 = find(XX3<0);  XX3(I3) = 0;
%% theta_record
theta_record = cell(1,3); 
A = zeros(1200,1200);% on the first view data
for i = 1:5
    hi = H1(i,:); si = find(hi == 1);
    ai = zeros(1200,1200);
    for j = 1:length(si)
        for k = 1:length(si)
            ai(si(j),si(k)) = 1;
        end
    end
    A = A + ai;
end

B = zeros(1300,1300);% on the second view data
for i = 1:5
    hi = H2(i,:); si = find(hi == 1);
    bi = zeros(1300,1300);
    for j = 1:length(si)
        for k = 1:length(si)
            bi(si(j),si(k)) = 1;
        end
    end
    B = B + bi;
end


C = zeros(2000,2000);% on the third view data
for i = 1:5
    hi = H3(i,:); si = find(hi == 1);
    ci = zeros(2000,2000);
    for j = 1:length(si)
        for k = 1:length(si)
            ci(si(j),si(k)) = 1;
        end
    end
    C = C + ci;
end

theta_record{1,1} = A; theta_record{1,2} = B; theta_record{1,3} = C;
%% R_record
n_vec = [1200,1300,2000]; 
R_record = cell(3,3);
for i = 1:3
    for j = 1:3
        if j ~= i
            R = zeros(n_vec(i),n_vec(j)); 
            Hr = H0_record{1,i}; Hc = H0_record{1,j};
            for k = 1:5
                r = zeros(n_vec(i),n_vec(j));
                hr = Hr(k,:); hc = Hc(k,:);
                sr = find(hr == 1); sc = find(hc == 1);
                for tr = 1:length(sr)
                    for tc = 1:length(sc)
                        r(sr(tr),sc(tc)) = 1;
                    end
                end
                R = R +r;
            end
             R_record{i,j} = R;    
        end
    end
end
% random permutation
p = randperm(500); p1 = randperm(1200); p2 = randperm(1300); p3 = randperm(2000);

for i = 1:500
    RX1(i,:) = XX1(p(i),:);  RXX1 = RX1; 
    RX2(i,:) = XX2(p(i),:);  RXX2 = RX2; 
    RX3(i,:) = XX3(p(i),:);  RXX3 = RX3; 
    WW(i,:) = W(p(i),:);
end

for i = 1:1200
    RXX1(:,i) = RX1(:,p1(i));
    RA(i,:) = A(p1(i),:); RAA = RA; RAA(:,i) = RA(:,p1(i)); 
    R12 = R_record{1,2}; R13 = R_record{1,3};
    R121(i,:) = R12(p1(i),:); 
    R131(i,:) = R13(p1(i),:); 
    RH1(:,i) = H1(:,p1(i));
end


for i = 1:1300
    RXX2(:,i) = RX2(:,p2(i));
    RB(i,:) = B(p2(i),:); RBB = RB; RBB(:,i) = RB(:,p2(i));
    R23 = R_record{2,3};
    R1211(:,i) = R121(:,p2(i));
    R231(i,:) = R23(p2(i),:); 
    RH2(:,i) = H2(:,p2(i));
end

for i = 1:2000
    RXX3(:,i) = RX3(:,p3(i));
    RC(i,:) = C(p3(i),:); RCC = RC; RCC(:,i) = RC(:,p3(i));
    R1311(:,i) = R131(:,p3(i));
    R2311(:,i) = R231(:,p3(i));  
    RH3(:,i) = H3(:,p3(i));
end

H_record = cell(1,3); H_record{1,1} = RH1; H_record{1,2} = RH2; H_record{1,3} = RH3;
theta_record1 = max(RAA + 0.1*randn(1200,1200),0); theta_record2 = max(RBB + 0.1*randn(1300,1300),0); theta_record3 = max(RCC + 0.1*randn(2000,2000),0);
theta_record{1,1} = (theta_record1 + theta_record1')/2; theta_record{1,2} = (theta_record2 + theta_record2')/2; theta_record{1,3} = (theta_record3 + theta_record3')/2;
R_record{1,2} = max(R1211 + 0.1*randn(1200,1300),0); R_record{1,3} = max(R1311 + 0.1*randn(1200,2000),0); R_record{2,3} = max(R2311 + 0.1*randn(1300,2000),0);
R_record{2,1} = R_record{1,2}'; R_record{3,1} = R_record{1,3}'; R_record{3,2} = R_record{2,3}';
W = WW;
X_record = cell(1,3);
X_record{1,1} = RXX1; X_record{1,2} = RXX2; X_record{1,3} = RXX3;
