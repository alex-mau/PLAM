close all;  clear all;  % clc;
addpath('../solvers','../subfun','./data','./spg_real')

train_data=importdata('leukemiatrain.txt');
test_data=importdata('leukemiatest.txt');

X_train=train_data(:,1:3051);
y_train=train_data(:,3052);
X_test=test_data(:,1:3051);
y_test=test_data(:,3052);
n_test=length(y_test);
n_train=length(y_train);
[n p]=size(X_train);
sigma=0.15;
delta=sqrt(2*log(p))*sigma;
y_train=y_train-sigma*randn(n_train,1);
y_test=y_test-sigma*randn(n_test,1);
X_train=zscore(X_train);
X_test=zscore(X_test);
y_train=zscore(y_train);
y_test=zscore(y_test);

eps=[ 1e-5 1e-6];

for i= 1:length(eps)
    fprintf('The epsilon= %f\n', eps(i));
    %%================== CPPA-PD================================
                para1.s = 100;
            rhos = 2.2203e+005;
        fprintf('Running CPPA-PD \n');
    fixp.rule = 'SRII';    fixp.eps = eps(i);    fixp.MAX = 20000;  fixp.detail = 0;
    X=X_train;  y=y_train;
    D=ones(p,1);
    for k = 1:p
        X(:,k)=X(:,k)/norm(X(:,k));
    end 
    t1=1.01;
    Mtype = 'PDM';    para1.tau = 1.2;  para1.r  = t1*rhos/para1.s;
    out1 = CPPA(X,D,y,delta,Mtype,para1,fixp);
    pred1=sign(X_test*out1.beta);
    real=sign(y_test);
    dif=pred1-real;
    fprintf('The accruarcy is %.2f, running time is %.2f  iteration is %d\n', length(find(dif==0))/n_test,out1.time,out1.iter)

    %=================== P-PLAM ==================================
    fprintf('Running P-PLAM \n');
    X=X_train;  y=y_train;
    [U,S,V]=svd(X,'econ');
    F=U*diag(1./diag(S))*U';
    X=F*X;
    D=ones(p,1)*sigma;
    para2.gamma =10;  para2.mu = 1;
    fixp.eps = eps(i);    fixp.MAX = 60000;  fixp.detail = 0;
    fixp.rule = 'SRII'; 
    out2 = PLAM(X,D,y,delta,para2,fixp);
    pred2=sign(X_test*out2.solution);
    real=sign(y_test);
    dif=pred2-real;
    fprintf('The accruarcy is %.2f, running time is %.2f iteration is %d\n', length(find(dif==0))/n_test,out2.time,out2.iter);
    
        %=================== P-LADM ======================================
    fprintf('Running P-LADM \n');
    fixp.rule = 'SRII';    fixp.eps = eps(i);    fixp.MAX = 20000;  fixp.detail = 0;
    X=X_train;  y=y_train;
    D=ones(p,1);
    for k = 1:p
        X(:,k)=X(:,k)/norm(X(:,k));
    end 
    para3.gamma =0.01;
      t1 = 2.1;
       rhos = 2.2203e+005;
    para3.mu = t1*para3.gamma*rhos;
    Mtype = 'original'; 
    out3 = PLADM(X,D,y,delta,Mtype,para3,fixp);
        pred3=sign(X_test*out3.beta);
    real=sign(y_test);
    dif=pred3-real;
    fprintf('The accruarcy is %.2f, running time is %.2f  iteration is %d\n', length(find(dif==0))/n_test,out3.time,out3.iter);
        fprintf('=============================== \n');
end
