function [W,H_record,X_record,theta_record, R_record] = JMF_synthetic_dataset3(mu)
% simulation synthetic dataset 3
% <Input>:
% mu:    The noise level
% <Outputs>:
% W:             The golden standard of low rank basis matrix
% H_record:      The golden standard of low rank coefficients matrices
% X_record:      The golden standard of original data matrices
% theta_record:  Cell data with network constrained matrix on each data as its element
% R_record:      Cell data with relationship constrained matrix between each pair of data as its element
index = 20; coph = 15; p1 =1/20; p2 = 1/20; p3 = 1/20;
n1 = 200; n2 = 150; n3 = 300;
W0 = zeros(2000,index);
for j = 1:2000
    for k = 1:index
        if j>= 1+(k-1)*(100-coph) && j <= 100+(k-1)*(100-coph)
            W0(j,k) = 1;
        end
    end
end
H10 = random('Binomial',1,p1,index,n1); H20 = random('Binomial',1,p2,index,n2); H30 = random('Binomial',1,p3,index,n3);
H_record = cell(1,3); H_record{1,1} = H10; H_record{1,2} = H20; H_record{1,3} = H30;
X1 = W0*H10 + mu*randn(2000,n1);
X2 = W0*H20 + mu*randn(2000,n2);
X3 = W0*H30 + mu*randn(2000,n3);
I1 = find(X1<0);  X1(I1) = 0;
I2 = find(X2<0);  X2(I2) = 0;
I3 = find(X3<0);  X3(I3) = 0;
%% theta_record
theta_record = cell(1,3); 
A = zeros(n1,n1);% on the first view data
for i = 1:index
    hi = H10(i,:); si = find(hi == 1);
    ai = zeros(n1,n1);
    for j = 1:length(si)
        for k = 1:length(si)
            ai(si(j),si(k)) = 1;
        end
    end
    A = A + ai;
end

A = max(A + 0.1*randn(n1,n1),0); A = (A + A')/2;

B = zeros(n2,n2);% on the second view data
for i = 1:index
    hi = H20(i,:); si = find(hi == 1);
    bi = zeros(n2,n2);
    for j = 1:length(si)
        for k = 1:length(si)
            bi(si(j),si(k)) = 1;
        end
    end
    B = B + bi;
end
B = max(B + 0.1*randn(n2,n2),0); B = (B + B')/2;

C = zeros(n3,n3);% on the third view data
for i = 1:index
    hi = H30(i,:); si = find(hi == 1);
    ci = zeros(n3,n3);
    for j = 1:length(si)
        for k = 1:length(si)
            ci(si(j),si(k)) = 1;
        end
    end
    C = C + ci;
end

C = max(C + 0.1*randn(n3,n3),0); C = (C + C')/2;
theta_record{1,1} = A; theta_record{1,2} = B; theta_record{1,3} = C;
%% R_record
n_vec = [n1,n2,n3]; 
R_record = cell(3,3);
for i = 1:2
    for j = i+1:3
        if j ~= i
            R = zeros(n_vec(i),n_vec(j)); 
            Hr = H_record{1,i}; Hc = H_record{1,j};
            for k = 1:index
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
             R_record{i,j} = max(R + 0.1*randn(n_vec(i),n_vec(j)),0);    
        end
    end
end

for i = 2:3
    for j = 1:i-1
        R_record{i,j} = R_record{j,i}';
    end
end

%% permutation rows
p = randperm(2000);
for i=1:2000
    RX1(i,:)=X1(p(i),:); 
    RX2(i,:)=X2(p(i),:); 
    RX3(i,:)=X3(p(i),:);
    W(i,:)=W0(p(i),:);
end
X_record= cell(1,3); X_record{1,1} = RX1; X_record{1,2} = RX2; X_record{1,3} = RX3;
