function Plot_Line(time1,std_time1,time2,std_time2,time3,std_time3,...
    error1,std_error1,error2,std_error2,error3,std_error3,...
    iters1,std_iters1,iters2,std_iters2,iters3,std_iters3)

figure(1)
hold('on');
plot(1:8,time1,'--b','LineWidth',2);
e1=errorbar(1:8,time1,std_time1,std_time1,'--b','LineWidth',1.5);
plot(1:8,time2,'--k','LineWidth',2);
e2=errorbar(1:8,time2,std_time2,std_time2,'--k','LineWidth',1.5);
plot(1:8,time3,'--r','LineWidth',2);
e3=errorbar(1:8,time3,std_time3,std_time3,'--r','LineWidth',1.5);
xlabel('i');
ylabel('Running time');
legend([e1 e2 e3],'P-LADM','CPPA-PD','P-PLAM');

figure(2)
hold('on');
plot(1:8,error1,'--b','LineWidth',2);
e1=errorbar(1:8,error1,std_error1,std_error1,'--b','LineWidth',1.5);
plot(1:8,error2,'--k','LineWidth',2);
e2=errorbar(1:8,error2,std_error2,std_error2,'--k','LineWidth',1.5);
plot(1:8,error3,'--r','LineWidth',2);
e3=errorbar(1:8,error3,std_error3,std_error3,'--r','LineWidth',1.5);
xlabel('i');
ylabel('Quality of solutions defined by \rho^2');
legend([e1 e2 e3],'P-LADM','CPPA-PD','P-PLAM');

figure(3)
hold('on');
plot(1:8,iters1,'--b','LineWidth',2);
e1=errorbar(1:8,iters1,std_iters1,std_iters1,'--b','LineWidth',1.5);
plot(1:8,iters2,'--k','LineWidth',2);
e2=errorbar(1:8,iters2,std_iters2,std_iters2,'--k','LineWidth',1.5);
plot(1:8,iters3,'--r','LineWidth',2);
e3=errorbar(1:8,iters3,std_iters3,std_iters3,'--r','LineWidth',1.5);
xlabel('i');
ylabel('Iterations');
legend([e1 e2 e3],'P-LADM','CPPA-PD','P-PLAM');

end
