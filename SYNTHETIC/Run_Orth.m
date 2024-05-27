close all;  clear all;
addpath('../solvers','../subfun')
cases = 'I';  datatype = 'state';
fprintf('Results for Orth Case %s with %s data.....\r\n',cases,datatype);
T_Orth=zeros(24,4);
m=10;  %m is the times of repeated experiments, I is uper bound of i.
sigma=0.05;

for  j=1:m
    for i =1:8
        fprintf( '======== Results of sigma = %2.3f ====== \n', sigma);
        p = i*2560; n = i*720;  T = i*80;
        [X,D,y,betas,delta] = DATA_DS_Orth(n,p,T,sigma,cases,datatype);
        fprintf(' %d  & %2.3f \r',i,norm(betas,1));

        fixp.eps = 2*10^-4;    fixp.MAX = 1000;  peps = 0.2;
        fixp.rule = 'SRII';    fixp.detail=0;
        %=================== P-LADM ======================================
        para1.gamma = 10;  para1.mu = 18; Mtype = 'original';
        out1 = PLADM(X,D,y,delta,Mtype,para1,fixp);
        [xm1, rhoo1, rhop1] = PostRho(out1.beta,X,y,betas,sigma,peps);
        fprintf('  P-LADM: %d & %2.2f & %2.2f & %2.3f (%2.3f) \\\\ \r',...
            out1.iter,out1.time,out1.obj,rhoo1,rhop1);
        T_Orth(i*3-2,:)=T_Orth(i*3-2,:)+[out1.iter,out1.time,rhoo1,rhop1];

        %=================== CPPA-PD ======================================
        para2.r = 11;  para2.s = 0.1; para2.tau = 1.2; Mtype = 'PDM';
        out2 = CPPA(X,D,y,delta,Mtype,para2,fixp);
        [xm4, rhoo2, rhop2] = PostRho(out2.beta,X,y,betas,sigma,peps);
        fprintf('CPPA-PD: %d & %2.2f & %2.2f & %2.3f (%2.3f)\\\\ \r',...
            out2.iter,out2.time,out2.obj,rhoo2,rhop2);
        T_Orth(i*3-1,:)=T_Orth(i*3-1,:)+[out2.iter,out2.time,rhoo2,rhop2];

        %=================== P-PLAM ======================================
        D=D*sigma;
        para3.gamma = 15;  para3.mu = 1; 
        out3 = PLAM(X,D,y,delta,para3,fixp);
        [xm3, rhoo3, rhop3] = PostRho(out3.solution,X,y,betas,sigma,peps);
        fprintf('P-PLAM: %d & %2.2f & %2.2f & %2.3f (%2.3f) \\\\ \r',...
            out3.iter,out3.time,out3.obj,rhoo3,rhop3);
        T_Orth(i*3,:)= T_Orth(i*3,:)+[out3.iter,out3.time,rhoo3,rhop3]; 
            
    end
        fprintf('\n\n');
        save([num2str(j) '_orth.mat'],'T_Orth');
end


fprintf('\n Running is completed!\n')