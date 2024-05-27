close all;  clear all;  % clc;
addpath('../solvers','../subfun','./data','./spg_real')
data=importdata('scaledbctrain.txt');
[n p]=size(data);
p=p-1;
sigma=0.05;
delta=sqrt(2*log(p))*sigma;

times=10;
n_test=20;

time=zeros(2,3);
acc=zeros(2,3);
iters=zeros(2,3);

eps=[1e-5 1e-6];
for i=1:length(eps)
    for j=1:times
        r=randperm( size(data,1)); 
        data=data(r, :);
        breast=data(:,1:p);
        breast_y=data(:,p+1);
        
        y=breast_y-sigma*randn(n,1);
        X=zscore(breast);
        y=zscore(y);
        
        n_train=n-n_test;
        X_train=X(1:n_train,:);
        y_train=y(1:n_train);
        X_test=X(n_train+1:n,:);
        y_test=y(n_train+1:n);
        fprintf('Data shuffled......\n');
%%================== CPPA-PD================================
            para1.s = 200;
            rhos = 4.5627e+005;
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
    time(i,1)=time(i,1)+out1.time;
    acc(i,1)=acc(i,1)+length(find(dif==0))/n_test;
    iters(i,1)=iters(i,1)+out1.iter;
    fprintf('The accruarcy is %.2f, running time is %.2f  iteration is %d\n', length(find(dif==0))/n_test,out1.time,out1.iter)

    %=================== P-PLAM ==================================
    fprintf('Running P-PLAM \n');
    X=X_train;  y=y_train;
    [U,S,V]=svd(X,'econ');
    F=U*diag(1./diag(S))*U';
    X=F*X;
    D=ones(p,1)*sigma;
    para2.gamma =10;  para2.mu = 1;
    fixp.eps = eps(i);    fixp.MAX = 20000;  fixp.detail = 0;
    fixp.rule = 'SRII'; 
    out2 = PLAM(X,D,y,delta,para2,fixp);
    % fprintf('i= %d Iter %d  Time %2.2f  Obj %2.2f \r',...
    %     i,out2.iter,out2.time,out2.obj);
    pred2=sign(X_test*out2.solution);
    real=sign(y_test);
    dif=pred2-real;
    time(i,2)=time(i,2)+out2.time;
    acc(i,2)=acc(i,2)+length(find(dif==0))/n_test;
    iters(i,2)=iters(i,2)+out2.iter;
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
       rhos = 4.5627e+005;
    para3.mu = t1*para3.gamma*rhos;
    Mtype = 'original'; 
    out3 = PLADM(X,D,y,delta,Mtype,para3,fixp);
        pred3=sign(X_test*out3.beta);
    real=sign(y_test);
    dif=pred3-real;
    time(i,3)=time(i,3)+out3.time;
    acc(i,3)=acc(i,3)+length(find(dif==0))/n_test;
    iters(i,3)=iters(i,3)+out3.iter;
    fprintf('The accruarcy is %.2f, running time is %.2f iteration %d \n', length(find(dif==0))/n_test,out3.time,out3.iter);
        fprintf('=============================== \n');
    end

end

time=time/times;
iters=iters/10;
acc=acc/10;

