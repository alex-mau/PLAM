function y=Bar(iters1,time1,acc1,iters2,time2,acc2)
%The input are CPPA-PD P-PLAM PLADM in order%

% 
% iters=diag([4079 3831 2513]);
% time=diag([1.15 0.83 0.78]);
% acc=diag([0.97 0.97 0.97]);
%sigma=0.05;

map=[0 0 0; 1 0 0; 0 0 1];
colormap(map)

subplot(2,3,1)
b1=bar(diag(iters1),'stacked');
set(b1,'edgecolor','none');
ylabel('Iterations')
set(gca, 'xticklabels', {'CPPA-PD', 'P-PLAM', 'P-LADM'})

subplot(2,3,2)
b2=bar(diag(time1),'stacked');
set(b2,'edgecolor','none');
ylabel('Computing time')
set(gca, 'xticklabels', {'CPPA-PD', 'P-PLAM', 'P-LADM'})

subplot(2,3,3)
b3=bar(diag(acc1),'stacked');
set(b3,'edgecolor','none');
ylabel('Accuracy')
set(gca, 'xticklabels', {'CPPA-PD', 'P-PLAM', 'P-LADM'})

subplot(2,3,4)
b4=bar(diag(iters2),'stacked');
set(b4,'edgecolor','none');
ylabel('Iterations')
set(gca, 'xticklabels', {'CPPA-PD', 'P-PLAM', 'P-LADM'})

subplot(2,3,5)
b5=bar(diag(time2),'stacked');
set(b5,'edgecolor','none');
ylabel('Computing time')
set(gca, 'xticklabels', {'CPPA-PD', 'P-PLAM', 'P-LADM'})

subplot(2,3,6)
b6=bar(diag(acc2),'stacked');
set(b6,'edgecolor','none');
ylabel('Accuracy')
set(gca, 'xticklabels', {'CPPA-PD', 'P-PLAM', 'P-LADM'})

end