close all;  clear all;  % clc;
addpath('../solvers','../subfun','./data','./spg_real')

data=[importdata('leukemiatrain.txt');importdata('leukemiatest.txt')];
leukemia=data(:,1:3051);
leukemia_y=data(:,3052);

data=importdata('lymph.dat');
lymph=data(:,1:4514);
lymph_y=data(:,4515);

data=importdata('scaledbctrain.txt');
breast=data(:,1:4919);
breast_y=data(:,4920);

colon=importdata('I2000.txt')';
colon_y=importdata('I2000b.txt');