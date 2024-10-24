clc 
clear all;
%addpath('../MY/')
%addpath('../MY/clusteringmeasure/');
%addpath('../MY/dataset2/dataset/');
%addpath('../MY/dataset/')
addpath('C:\Users\86180\Desktop\MY2\clusteringmeasure')
addpath('C:\Users\86180\Desktop\MY2\dataset2\dataset');
addpath('C:\Users\86180\Desktop\MY2\dataset_0.15');
dataname = 'BBC_685';
load([dataname,'.mat']);

num_class = length(unique(truth));%返回聚类数
num_views = length(X);

for v= 1:num_views
    X{v} = normalize(X{v},'L2');
end

Result = [];

Lambda_c=0.01;
n=0.00032;
Lambda_d=10000;
mu=0.06;         
epsilon=10^-5;
max_time=500;           
gamma=0.1;

for lci = 1 :length(Lambda_c)
    for ldi = 1 : length(Lambda_d)
        for ln = 1 : length(n)
            
            lambda_c = Lambda_c(lci);
            lambda_d = Lambda_d(ldi);
            nn=n(ln);
            Noise_ratios=[1,1,1,1,1,1]*nn;
            tic;
            [Z,t,W]=WWYWAY(X,lambda_c,lambda_d,mu,epsilon,max_time,gamma,Noise_ratios);
            whole_time=toc;
            disp(whole_time);
            [result] = CluEst_Spec2(Z,num_class,truth);
            Result = [Result,result]; 
        end
    end
end

Result
t


