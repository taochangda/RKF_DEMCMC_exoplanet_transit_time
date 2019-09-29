% Tao Changda, June, 2019, All rights reserved.

%% initialization and pre-process 初始化，预处理

clear all
close all
clc
% % clear Nbody_C_7Np_RKF78.mexw64
% % delete Nbody_C_7Np_RKF78.mexw64
% % mex Nbody_C_7Np_RKF78.c


%-----begin------basic parameter 程序运行基础参数---------
source_file_path='C:\works\Nbody\'; %path of datas 数据源所在路径
startup_file=[source_file_path,'for_startup20190301.mat'];
ncores=700;%Number of distributed cores for run matlab  能调用的分布式计算的总核(或者线程)数，每核(线程)会运行一个分布式matlab进程
reboot_count=1200;
%-----end------basic parameter 程序运行基础参数---------

%-----begin-----parameters of exoplanet system 星系系统参数-----------
load([source_file_path,'transit_data_7Np_XXX.mat']) %including transit time oberseved, range of fitting track parameter 观测值及基础参数包括观测的凌星时刻，包括观测的凌星时刻，参数取值范围
m=35;%number of fitting parameter 拟合的参数个数
Nplanet=7; %number of planet 行星个数
mc=0.0802; %mass of star（times of sun）恒星质量(单位：太阳质量倍数）
Rs=0.0005382; %radius of star 恒星半径
err=1e-12; %precision of RKF and so on. 1e-12 for RKF78 and 1e-10 for RKF45?  I am not sure.  RKF算法精度
%the more little of err, longer time will be consumed for running. And err could not be set to too big for the correct result 
time_precise=1e-5; %minimal time interval for searching transit time 寻找凌星时刻方法中对时间步长的精度控制，如果不控制，可能会出现时间步长很小，迭代步数过多，消耗过多的计算时间。详见相关程序说明
n=0; %number of transit time oberseved 所有行星总凌星时刻个数
%如下：算出最长的观测时间max_time，积分计算的时间tstop比max_time稍大
max_time=0;
for i=1:Nplanet
    n=n+size(planet_data{i},1);
    if (max_time<max(planet_data{i}(:,2 )))
        max_time=max(planet_data{i}(:,2));
    end
end
n=size(planet_data{1},1)+size(planet_data{2},1)+size(planet_data{3},1)+size(planet_data{4},1)+size(planet_data{5},1)+size(planet_data{6},1)+size(planet_data{7},1);
tstop=1.04*max_time/365.2422*2*pi; % tstop：integration time（transform to circular measure）,a little bigger than max_time, and 1.04 is a value of experience. It will be consumed  long time if tstop is set to be big
%积分总时间(转化成弧度)，要比总观测时间max_time稍放大，1.04只是个经验值，不要过大，太大的话，耗时过多
%----end-------parameters of exoplanet system 星系系统参数-----------


%----begin----parameter of Differential Evolution（DE）algorithm---DE算法参数-----------
n_link=4*ncores;  % n_link：1) population of DE algorithm, 2) number of Markov Chains of DE-MCMC, n_link is reset to a small value burn_phase and conve_phase. 
%3),4 times of ncores is a empirical value. what is the value of n_link,
%I am not sure, please choose the value base on your compute power. 
%1)DE算法中,种群规模数，2)DE-MCMC算法中算法马尔可夫链数，此处值是设置给DE算法的，在阶段2和3(burn_phase和conve_phase会重置)
%3）4倍ncores不是固定取值 ，请根据实际的计算能力情况，n_link不能过小，可能取值比较大会好，但这只是经验，并不是确定的。
conve_N=20; %number of links selected for DE-MCMC(burn_phase and conve_phase) after DE phase. I am not sure what is the right value. but it seems that it should not to set to big.
% Another method, set different values for  conve_N to 10, 20, 30, 40,...
% and caculate respectively in different machines
%DE算法结束后，选出conve_N组参数，在阶段2和3(burn_phase和conve_phase),用于DE-MCMC算法进行收敛计算。
%conve_N值设置多少我并不是很确定，但好像太大了不好收敛。有一种方法，可以在DE算法结束后，分别取conve_N为10，20，30，40。。用多台机器分别计算，看如一组收敛
N_DE=1*1e4; % minum step of algorithm DE, DE算法最少迭代步数，根据实际情况设置；
N_var=500; %maximum step of mutation.
%DE算法每次迭代的变异步骤中，因为参数都有取值范围，不能保证变异操作后的参数值都取值范围内，所以要进行多次变异操作。需要设置一个最大的变异操作次数。如果超过此次数，则保留原值
gamma_ini=2.38/sqrt(2*m);
gamma=gamma_ini*(1+1e-4*randn); %initailize gamma
max_dist=0.05; %control parameter of DE algorithm: Subtraction difference between maximum reduced-chi and minimun reduced-chi
%DE算法停止控制参数之一，表示最大和最小re_chi之差
CR=0.3;%cross rate 交叉率
aim_acc_p=0.1; %one of the parameter for updating gamma
%目标接受百分比，这里只是个经验值，因为n_link比较大，所以这里不适合设置过大。
lower_acc_ratio=0.8; %one of the parameter for updating gamma
%gamma每次迭代变化的控制参数之一，aim_acc_p*lower_acc_ratio;
upper_acc_ratio=1.2; %one of the parameter for updating gamma
%gamma每次迭代变化的控制参数之一，aim_acc_p*upper_acc_ratio;
% pre_acc_p=ones(1,gamma_reset_count);
pre_acc_p=1;
acc_p=1; %初始化
gamma_reset_count=100; %reset gamma after gamma_reset_count step for preventing very small gamma 必须是偶数，为避免迭代过程中gamma过小，经过一定步骤后reset gamma值为1
half_gamma_reset_count=gamma_reset_count/2;
%-----end----parameter of Differential Evolution（DE）algorithm--DE算法参数-----------


%-----begin------initialize fitting parameter of planets初始化拟合参数-----------
loop_out=1;
temp_para=zeros(n_link,m);
temp_var1=rand(n_link,m);
p_para=range_lower+(range_upper-range_lower).*temp_var1;%初始化m个拟合的参数值，在取值范围内取随机数


delete(gcp('nocreate'));
pool=parpool('myMJS', ncores);%---start parallel job, refer matlab guide for more information 启动并行任务
addAttachedFiles(pool,[source_file_path,'Nbody_C_7Np_RKF78.mexw64']);


%parallel computing, and the speed base on the distrubited  compute
%resource and networks
%并行计算，需要调用大量分布式计算资源，计算快慢取决和计算资源和网络情况
parfor loop_1=1:n_link 
    [tr,error_flag]=Nbody_C_7Np_RKF78(p_para(loop_1,:),mc,Rs,tstop,time_precise,err); %CALL funtion Nbody_C_7Np_RKF78  调用函数Nbody_C_7Np_RKF78进行凌星时间求解
    chi_i(loop_1)=sum(((tr(planet_data{1}(:,1)+1,1)-planet_data{1}(:,2))./planet_data{1}(:,3)).^2)+sum(((tr(planet_data{2}(:,1)+1,2)-planet_data{2}(:,2))./planet_data{2}(:,3)).^2)+sum(((tr(planet_data{3}(:,1)+1,3)-planet_data{3}(:,2))./planet_data{3}(:,3)).^2)+sum(((tr(planet_data{4}(:,1)+1,4)-planet_data{4}(:,2))./planet_data{4}(:,3)).^2)+sum(((tr(planet_data{5}(:,1)+1,5)-planet_data{5}(:,2))./planet_data{5}(:,3)).^2)+sum(((tr(planet_data{6}(:,1)+1,6)-planet_data{6}(:,2))./planet_data{6}(:,3)).^2)+sum(((tr(planet_data{7}(:,1)+1,7)-planet_data{7}(:,2))./planet_data{7}(:,3)).^2);
% chi_i: errors 
end
delete(gcp('nocreate')); %---end parallel job, refer matlab guide for more information 结束并行任务

DE_phase=1;
burn_phase=1;
conve_phase=1;

save(startup_file,'-v7.3');
%---end------initialize fitting parameter of planets初始化拟合参数-----------


