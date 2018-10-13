function [W,H_record, X_record,theta_record, R_record] = JMF_synthetic_dataset2(mu)
% simulation synthetic dataset 2
% <Input>:
% mu:    The noise level
% <Outputs>:
% W:             The golden standard of low rank basis matrix
% H_record:      The golden standard of low rank coefficients matrices
% X_record:      The golden standard of original data matrices
% theta_record:  Cell data with network constrained matrix on each data as its element
% R_record:      Cell data with relationship constrained matrix between each pair of data as its element

p = 1/10; coph1 = 0; coph2 = 5; coph3 = 10;
W0 = random('Binomial',1,p,1000,10);
H01 = zeros(10,200);
for j = 1:10
    for k = 1:200
        if k >= 1+(j-1)*(20-coph1) && k <= 20+(j-1)*(20-coph1)
            H01(j,k) = 1;
        end
    end
end

H02 = zeros(10,300);
for j = 1:10
    for k = 1:300
        if k >= 1+(j-1)*(30-coph2) && k <= 30+(j-1)*(30-coph2)
            H02(j,k) = 1;
        end
    end
end

H03 = zeros(10,500);
for j = 1:10
    for k = 1:500
        if k >= 1+(j-1)*(50-coph3) && k <= 50+(j-1)*(50-coph3)
            H03(j,k) = 1;
        end
    end
end
H0_record = cell(1,3);H0_record{1,1} = H01; H0_record{1,2} = H02; H0_record{1,3} = H03;
X1 = W0*H01 + mu*randn(1000,200);
X2 = W0*H02 + mu*randn(1000,300);
X3 = W0*H03 + mu*randn(1000,500);
I1 = find(X1<0);  X1(I1) = 0;
I2 = find(X2<0);  X2(I2) = 0;
I3 = find(X3<0);  X3(I3) = 0;
%% theta_record
theta_record = cell(1,3); 
A = zeros(200);% on the first view
for i = 1:10
    hi = H01(i,:); si = find(hi == 1);
    ai = zeros(200);
    for j = 1:length(si)
        for k = 1:length(si)
            ai(si(j),si(k)) = 1;
        end
    end
    A = A + ai;
end

A = max(A + 0.1*randn(200),0); A = (A + A')/2;

B = zeros(300);% on the second view
for i = 1:10
    hi = H02(i,:); si = find(hi == 1);
    bi = zeros(300);
    for j = 1:length(si)
        for k = 1:length(si)
            bi(si(j),si(k)) = 1;
        end
    end
    B = B + bi;
end
B = max(B + 0.1*randn(300),0); B = (B + B')/2;

C = zeros(500); % on the third view
for i = 1:10
    hi = H03(i,:); si = find(hi == 1);
    ci = zeros(500);
    for j = 1:length(si)
        for k = 1:length(si)
            ci(si(j),si(k)) = 1;
        end
    end
    C = C + ci;
end

C = max(C + 0.1*randn(500),0); C = (C + C')/2;

%% R_record
n_vec = [200,300,500]; 
R_record = cell(3,3);
for i = 1:2
    for j = i+1:3
        if j ~= i
            R = zeros(n_vec(i),n_vec(j)); 
            Hr = H0_record{1,i}; Hc = H0_record{1,j};
            for k = 1:10
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

%% permutation columns
p1 = randperm(200); p2 = randperm(300); p3 = randperm(500);
for i = 1:200
    RX1(:,i) = X1(:,p1(i));
    RA(i,:) = A(p1(i),:);
    R12 = R_record{1,2}; R13 = R_record{1,3};
    R121(i,:) = R12(p1(i),:); 
    R131(i,:) = R13(p1(i),:); 
    RH1(:,i) = H01(:,p1(i));
end
for j = 1:200
    RAA(:,j) = RA(:,p1(j)); 
end

for i = 1:300
    RX2(:,i) = X2(:,p2(i));
    RB(i,:) = B(p2(i),:); 
    R23 = R_record{2,3};
    R1211(:,i) = R121(:,p2(i));
    R231(i,:) = R23(p2(i),:); 
    RH2(:,i) = H02(:,p2(i));
end

for j = 1:300
    RBB(:,j) = RB(:,p2(j));
end


for i=1:500
    RX3(:,i) = X3(:,p3(i));
    RC(i,:) = C(p3(i),:); 
    R1311(:,i) = R131(:,p3(i));
    R2311(:,i) = R231(:,p3(i));  
    RH3(:,i) = H03(:,p3(i));
end

for j = 1:500
   RCC(:,j) = RC(:,p3(j));
end
W = W0;
X_record= cell(1,3); X_record{1,1} = RX1; X_record{1,2} = RX2; X_record{1,3} = RX3;
H_record = cell(1,3); H_record{1,1} = RH1; H_record{1,2} = RH2; H_record{1,3} = RH3;
theta_record{1,1} = RAA; theta_record{1,2} = RBB; theta_record{1,3} = RCC;
R_record{1,2} = R1211; R_record{1,3} = R1311; R_record{2,3} = R2311;
R_record{2,1} = R_record{1,2}'; R_record{3,1} = R_record{1,3}'; R_record{3,2} = R_record{2,3}';
