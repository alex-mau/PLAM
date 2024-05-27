addpath('../solvers','../subfun')
close all;  clear all;
cases = 'II';  datatype = 'state';
fprintf('Results for Orth Case %s with %s data.....\r\n',cases,datatype);
SIG = [0.05];
for j = 1:length(SIG)
    sigma = SIG(j);
    fprintf( '======== Results of sigma = %2.3f ====== \n', sigma);
    for i = 1:1
        p = i*5120; n = i*1440;  T = i*160;
        [X,D,y,betas,delta] = DATA_DS_Orth(n,p,T,sigma,cases,datatype);
         betas=sign(betas).*(1+rand(p,1));  y=X*betas+sigma*randn(n,1);
        fprintf(' %d  & %2.3f \r',i,norm(betas,1));
        
        fixp.eps = 2*10^-4;    fixp.MAX = 2000;  peps = 0.2;
        fixp.rule = 'SRII';   fixp.detail = 0 ;
        %%=================== P-LADM ======================================
        para1.gamma = 10;  para1.mu = 18; Mtype = 'original';
        out1 = PLADM(X,D,y,delta,Mtype,para1,fixp);
         [xm1, rhoo1, rhop1] = PostRho(out1.beta,X,y,betas,sigma,peps);
        fprintf('  P-LADM: %d & %2.2f & %2.2f & %2.3f (%2.3f) \\\\ \r',...
            out1.iter,out1.time,out1.obj,rhoo1,rhop1);
        
        %%=================== P-PLAM ======================================

        D=ones(size(D)).*sigma;
        para3.gamma = 15;  para3.mu = 1; fixp.detail = 0; fixp.eps = 10^-4;
        out3 = PLAM(X,D,y,delta,para3,fixp);
         [xm3, rhoo3, rhop3] = PostRho(out3.solution,X,y,betas,sigma,peps);
        fprintf('P-PLAM: %d & %2.2f & %2.2f & %2.3f (%2.3f) \\\\ \r',...
            out3.iter,out3.time,out3.obj,rhoo3,rhop3);
        
    end

figure(1)
subplot(2,1,1);
axis equal
stem(betas,'ks');
hold on
stem(out3.solution,'ro');
hold on
stem(out1.beta,'b+');


set(gca,'XLim',[0 p]);

subplot(2,1,2);
axis equal
stem(betas,'ks');
hold on
stem(xm3,'ro');
hold on
stem(xm1,'b+');
set(gca,'XLim',[0 p]);


fprintf('\n');
end
fprintf('Running is completed!\n')

