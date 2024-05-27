close all;  clear all;  % clc;
addpath('../solvers','../subfun')
cases = 'I';  datatype = 'state';
fprintf('Results for Unit Case %s with %s data.....\n',cases,datatype);

SIG = [0.05];

for j = 1:length(SIG)
    sigma = SIG(j);
    fprintf( '======== Results of sigma = %2.3f ====== \r\n', sigma);
    for i = 1:1
        
        p = i*5120; n = i*1440;  T = i*160;
        [X,D,y,betas,delta] = DATA_DS_Unit(n,p,T,sigma,cases,datatype);
        fprintf(' %d  & %2.3f \r',i,norm(betas,1));
        
        fixp.eps = 2*10^-4;    fixp.MAX = 1000;  peps = 0.2;
        fixp.rule = 'SRII';    rhos = 70; fixp.detail=0;
        %%=================== P-LADM ======================================
        %         para1.gamma = 0.2;  % para1.mu = 15;
        para1.gamma = 3/(log(p));
        %         para1.gamma = 20/(sqrt(p)+log(p));
        Mtype = 'original'; t1 = 2.0;
        para1.mu = t1*para1.gamma*rhos;
        out1 = PLADM(X,D,y,delta,Mtype,para1,fixp);
        [xm1, rhoo1, rhop1] = PostRho(out1.beta,X,y,betas,sigma,peps);
        fprintf('  P-LADM: %d & %2.2f & %2.2f & %2.3f (%2.3f)\r',...
            out1.iter,out1.time,out1.obj,rhoo1,rhop1);
        
        %=================== P-PLAM ==================================
        [U,S,V]=svd(X);
        F=U*inv(S(:,1:n))*U';
        X=F*X;
        y=F*y;
        D=zeros(p,1);
        for ki=1:p
            D(ki)=norm(X(:,ki));
        end
        D=D.*sigma;  
        para6.gamma =35;  para6.mu = 1; fixp.detail = 0;      
        out6 = PLAM(X,D,y,delta,para6,fixp);
        [xm6, rhoo6, rhop6] = PostRho(out6.solution,X,y,betas,sigma,peps);
        fprintf('  P-PLAM: %d & %2.2f & %2.2f & %2.3f (%2.3f)\r',...
            out6.iter,out6.time,out6.obj,rhoo6,rhop6);
       
        
    end
    
    fprintf('\n');
end
figure(1)
subplot(2,1,1);
axis equal
stem(betas,'ks');
hold on
stem(out6.solution,'ro');
hold on
stem(out1.beta,'b+');


set(gca,'XLim',[0 p]);

subplot(2,1,2);
axis equal
stem(betas,'ks');
hold on
stem(xm6,'ro');
hold on
stem(xm1,'b+');
set(gca,'XLim',[0 p]);
fprintf('Running is completed!\r\n')

