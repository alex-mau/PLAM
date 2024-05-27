% leukemia dataset
% results when sigma=0.05
figure(1)
iters1=[4079 3831 2513];
time1=[1.15 0.83 0.78];
acc1=[0.97 0.97 0.97];

iters2=[20000 6822 20000];
time2=[4.93 1.42 6.11];
acc2=[0.97 0.97 0.97];
Bar(iters1,time1,acc1,iters2,time2,acc2)

% results when sigma=0.10
figure(2)
iters1=[4839 2990 4806];
time1=[1.32 0.63 1.48];
acc1=[0.94 0.97 0.94];

iters2=[20000 6006 20000];
time2=[4.96 1.33 6.14];
acc2=[0.94 0.97 0.94];
Bar(iters1,time1,acc1,iters2,time2,acc2)

% results when sigma=0.15
figure(3)
iters1=[4800 3546 4489];
time1=[1.31 0.74 1.38];
acc1=[0.94 0.97 0.94];

iters2=[20000 5821 20000];
time2=[4.88 1.20 6.09];
acc2=[0.97 0.97 0.94];
Bar(iters1,time1,acc1,iters2,time2,acc2)

%breast dataset
% results when sigma=0.05
figure(1)

iters1=[6395.8 3225 2327.1];
time1=[4.48 1.93 2.15];
acc1=[0.57 0.68 0.64];

iters2=[20000 7844.3 20000];
time2=[13.34 4.71 18.24];
acc2=[0.61 0.66 0.63];
Bar(iters1,time1,acc1,iters2,time2,acc2)

% results when sigma=0.10
figure(2)

iters1=[4373.2 2329.1 3081.7];
time1=[3.15 1.38 2.82];
acc1=[0.66 0.67 0.63];

iters2=[20000 6715.9 20000];
time2=[13.11 3.99 18.03];
acc2=[0.65 0.69 0.64];
Bar(iters1,time1,acc1,iters2,time2,acc2)

% results when sigma=0.15
figure(3)

iters1=[4299.6 2668.6 3943.6];
time1=[3.13 1.59 3.62];
acc1=[0.64 0.68 0.63];

iters2=[20000 2780.6 20000];
time2=[13.42 1.69 18.42];
acc2=[0.64 0.63 0.61];
Bar(iters1,time1,acc1,iters2,time2,acc2)
