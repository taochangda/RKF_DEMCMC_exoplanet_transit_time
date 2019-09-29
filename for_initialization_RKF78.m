% Tao Changda, June, 2019, All rights reserved.

%% initialization and pre-process ��ʼ����Ԥ����

clear all
close all
clc
% % clear Nbody_C_7Np_RKF78.mexw64
% % delete Nbody_C_7Np_RKF78.mexw64
% % mex Nbody_C_7Np_RKF78.c


%-----begin------basic parameter �������л�������---------
source_file_path='C:\works\Nbody\'; %path of datas ����Դ����·��
startup_file=[source_file_path,'for_startup20190301.mat'];
ncores=700;%Number of distributed cores for run matlab  �ܵ��õķֲ�ʽ������ܺ�(�����߳�)����ÿ��(�߳�)������һ���ֲ�ʽmatlab����
reboot_count=1200;
%-----end------basic parameter �������л�������---------

%-----begin-----parameters of exoplanet system ��ϵϵͳ����-----------
load([source_file_path,'transit_data_7Np_XXX.mat']) %including transit time oberseved, range of fitting track parameter �۲�ֵ���������������۲������ʱ�̣������۲������ʱ�̣�����ȡֵ��Χ
m=35;%number of fitting parameter ��ϵĲ�������
Nplanet=7; %number of planet ���Ǹ���
mc=0.0802; %mass of star��times of sun����������(��λ��̫������������
Rs=0.0005382; %radius of star ���ǰ뾶
err=1e-12; %precision of RKF and so on. 1e-12 for RKF78 and 1e-10 for RKF45?  I am not sure.  RKF�㷨����
%the more little of err, longer time will be consumed for running. And err could not be set to too big for the correct result 
time_precise=1e-5; %minimal time interval for searching transit time Ѱ������ʱ�̷����ж�ʱ�䲽���ľ��ȿ��ƣ���������ƣ����ܻ����ʱ�䲽����С�������������࣬���Ĺ���ļ���ʱ�䡣�����س���˵��
n=0; %number of transit time oberseved ��������������ʱ�̸���
%���£������Ĺ۲�ʱ��max_time�����ּ����ʱ��tstop��max_time�Դ�
max_time=0;
for i=1:Nplanet
    n=n+size(planet_data{i},1);
    if (max_time<max(planet_data{i}(:,2 )))
        max_time=max(planet_data{i}(:,2));
    end
end
n=size(planet_data{1},1)+size(planet_data{2},1)+size(planet_data{3},1)+size(planet_data{4},1)+size(planet_data{5},1)+size(planet_data{6},1)+size(planet_data{7},1);
tstop=1.04*max_time/365.2422*2*pi; % tstop��integration time��transform to circular measure��,a little bigger than max_time, and 1.04 is a value of experience. It will be consumed  long time if tstop is set to be big
%������ʱ��(ת���ɻ���)��Ҫ���ܹ۲�ʱ��max_time�ԷŴ�1.04ֻ�Ǹ�����ֵ����Ҫ����̫��Ļ�����ʱ����
%----end-------parameters of exoplanet system ��ϵϵͳ����-----------


%----begin----parameter of Differential Evolution��DE��algorithm---DE�㷨����-----------
n_link=4*ncores;  % n_link��1) population of DE algorithm, 2) number of Markov Chains of DE-MCMC, n_link is reset to a small value burn_phase and conve_phase. 
%3),4 times of ncores is a empirical value. what is the value of n_link,
%I am not sure, please choose the value base on your compute power. 
%1)DE�㷨��,��Ⱥ��ģ����2)DE-MCMC�㷨���㷨����ɷ��������˴�ֵ�����ø�DE�㷨�ģ��ڽ׶�2��3(burn_phase��conve_phase������)
%3��4��ncores���ǹ̶�ȡֵ �������ʵ�ʵļ������������n_link���ܹ�С������ȡֵ�Ƚϴ��ã�����ֻ�Ǿ��飬������ȷ���ġ�
conve_N=20; %number of links selected for DE-MCMC(burn_phase and conve_phase) after DE phase. I am not sure what is the right value. but it seems that it should not to set to big.
% Another method, set different values for  conve_N to 10, 20, 30, 40,...
% and caculate respectively in different machines
%DE�㷨������ѡ��conve_N��������ڽ׶�2��3(burn_phase��conve_phase),����DE-MCMC�㷨�����������㡣
%conve_Nֵ���ö����Ҳ����Ǻ�ȷ����������̫���˲�����������һ�ַ�����������DE�㷨�����󣬷ֱ�ȡconve_NΪ10��20��30��40�����ö�̨�����ֱ���㣬����һ������
N_DE=1*1e4; % minum step of algorithm DE, DE�㷨���ٵ�������������ʵ��������ã�
N_var=500; %maximum step of mutation.
%DE�㷨ÿ�ε����ı��첽���У���Ϊ��������ȡֵ��Χ�����ܱ�֤���������Ĳ���ֵ��ȡֵ��Χ�ڣ�����Ҫ���ж�α����������Ҫ����һ�����ı��������������������˴���������ԭֵ
gamma_ini=2.38/sqrt(2*m);
gamma=gamma_ini*(1+1e-4*randn); %initailize gamma
max_dist=0.05; %control parameter of DE algorithm: Subtraction difference between maximum reduced-chi and minimun reduced-chi
%DE�㷨ֹͣ���Ʋ���֮һ����ʾ������Сre_chi֮��
CR=0.3;%cross rate ������
aim_acc_p=0.1; %one of the parameter for updating gamma
%Ŀ����ܰٷֱȣ�����ֻ�Ǹ�����ֵ����Ϊn_link�Ƚϴ��������ﲻ�ʺ����ù���
lower_acc_ratio=0.8; %one of the parameter for updating gamma
%gammaÿ�ε����仯�Ŀ��Ʋ���֮һ��aim_acc_p*lower_acc_ratio;
upper_acc_ratio=1.2; %one of the parameter for updating gamma
%gammaÿ�ε����仯�Ŀ��Ʋ���֮һ��aim_acc_p*upper_acc_ratio;
% pre_acc_p=ones(1,gamma_reset_count);
pre_acc_p=1;
acc_p=1; %��ʼ��
gamma_reset_count=100; %reset gamma after gamma_reset_count step for preventing very small gamma ������ż����Ϊ�������������gamma��С������һ�������reset gammaֵΪ1
half_gamma_reset_count=gamma_reset_count/2;
%-----end----parameter of Differential Evolution��DE��algorithm--DE�㷨����-----------


%-----begin------initialize fitting parameter of planets��ʼ����ϲ���-----------
loop_out=1;
temp_para=zeros(n_link,m);
temp_var1=rand(n_link,m);
p_para=range_lower+(range_upper-range_lower).*temp_var1;%��ʼ��m����ϵĲ���ֵ����ȡֵ��Χ��ȡ�����


delete(gcp('nocreate'));
pool=parpool('myMJS', ncores);%---start parallel job, refer matlab guide for more information ������������
addAttachedFiles(pool,[source_file_path,'Nbody_C_7Np_RKF78.mexw64']);


%parallel computing, and the speed base on the distrubited  compute
%resource and networks
%���м��㣬��Ҫ���ô����ֲ�ʽ������Դ���������ȡ���ͼ�����Դ���������
parfor loop_1=1:n_link 
    [tr,error_flag]=Nbody_C_7Np_RKF78(p_para(loop_1,:),mc,Rs,tstop,time_precise,err); %CALL funtion Nbody_C_7Np_RKF78  ���ú���Nbody_C_7Np_RKF78��������ʱ�����
    chi_i(loop_1)=sum(((tr(planet_data{1}(:,1)+1,1)-planet_data{1}(:,2))./planet_data{1}(:,3)).^2)+sum(((tr(planet_data{2}(:,1)+1,2)-planet_data{2}(:,2))./planet_data{2}(:,3)).^2)+sum(((tr(planet_data{3}(:,1)+1,3)-planet_data{3}(:,2))./planet_data{3}(:,3)).^2)+sum(((tr(planet_data{4}(:,1)+1,4)-planet_data{4}(:,2))./planet_data{4}(:,3)).^2)+sum(((tr(planet_data{5}(:,1)+1,5)-planet_data{5}(:,2))./planet_data{5}(:,3)).^2)+sum(((tr(planet_data{6}(:,1)+1,6)-planet_data{6}(:,2))./planet_data{6}(:,3)).^2)+sum(((tr(planet_data{7}(:,1)+1,7)-planet_data{7}(:,2))./planet_data{7}(:,3)).^2);
% chi_i: errors 
end
delete(gcp('nocreate')); %---end parallel job, refer matlab guide for more information ������������

DE_phase=1;
burn_phase=1;
conve_phase=1;

save(startup_file,'-v7.3');
%---end------initialize fitting parameter of planets��ʼ����ϲ���-----------


