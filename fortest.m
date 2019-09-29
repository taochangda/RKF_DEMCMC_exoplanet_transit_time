%%
%预处理
clc
clear all
close all
delete(gcp('nocreate'))


[notfound,warnings]=loadlibrary('C:\works\Nbody\codegen\dll\rkf78_7Np\rkf78_7Np.dll','C:\works\Nbody\codegen\dll\rkf78_7Np\rkf78_7Np.h');
[notfound,warnings]=loadlibrary('C:\works\Nbody\codegen\dll\rkf78_forC\rkf78_forC.dll','C:\works\Nbody\codegen\dll\rkf78_forC\rkf78_forC.h');
[notfound,warnings]=loadlibrary('C:\works\Nbody\codegen\dll\NbodySubProg\NbodySubProg.dll','C:\works\Nbody\codegen\dll\NbodySubProg\NbodySubProg.h');


ncores=feature('numCores');
% ncores=480;
parpool(ncores)

parfor i=1:ncores
    [notfound,warnings]=loadlibrary('C:\works\Nbody\codegen\dll\rkf78_7Np\rkf78_7Np.dll','C:\works\Nbody\codegen\dll\rkf78_7Np\rkf78_7Np.h');
    [notfound,warnings]=loadlibrary('C:\works\Nbody\codegen\dll\rkf78_forC\rkf78_forC.dll','C:\works\Nbody\codegen\dll\rkf78_forC\rkf78_forC.h');
    [notfound,warnings]=loadlibrary('C:\works\Nbody\codegen\dll\NbodySubProg\NbodySubProg.dll','C:\works\Nbody\codegen\dll\NbodySubProg\NbodySubProg.h');
end

%%
%观测值及基础参数
load lingxing_2  % 包括观测的凌星时刻，以及rkf78算法中使用到的常数a（多加一行，为了后续计算不出错，a :14*12）
% load lingxing
% 观测值起点t不等于0,需要减去一个常数，保证每个星的第一凌星时刻都小于其周期(p(1),p(2),...)
planet_1_lx(:,2)=planet_1_lx(:,2)-65;
planet_2_lx(:,2)=planet_2_lx(:,2)-65;

%%基础参数
Nplanet=2;
gamma=2.38/sqrt(2*13);
mc=0.55;
Rs=0.002438;
err=1e-16;
n=size(planet_1_lx,1)+size(planet_2_lx,1);
m=13;%参数个数
% [p(1),p(2),e(1),e(2),ci(1),ci(2),w(1),w(2),cm0(1),cm0(2),mp(1),mp(2),domega]=deal(19.2259140634,39.0154721795,0.0631304497,0.0681090770,91.2178250485,90.0606079171,356.6852434554,169.2203655656,225.3157586887,252.8455786655,44.2383233216,30.4181676605,1.4998923879);
% ap=(mc+mp).^(1/3).*(p./365.2422).^(2/3);
% delta_ci=atan(Rs./ap);
% % delta_ci=mod(atan(Rs./ap),2*pi);
% ci_l=pi/2-delta_ci;
% ci_u=pi/2+delta_ci;

% 各参数范围，统一以转换单位后的范围
%p1:19.1-19.3, p2:38.9-39.1, e1:0-0.0.7, e2:0-0.07, ci1:90-icrit1---90+icrit1, ci2:90-icrit2---90+icrit2,
% w1:350-360, w2:165-170, cm1:0-360, cm2:0-360,mp1:40-46, mp2:28-32, domega:0-20
range_lower=[7.2,10.91,0,0,88*pi/180,88.5*pi/180,0,0,0,0,1*3e-6,1*3e-6,0];
range_upper=[7.21,10.92,0.1,0.1,92*pi/180,91.5*pi/180,2*pi,2*pi,2*pi,2*pi,15*3e-6,15*3e-6,10*pi/180];

% range_lower=[19.1,38.9,0,0,ci_l(1),ci_l(2),350*pi/180,165*pi/180,0,0,40*3e-6,28*3e-6,0];
% range_upper=[19.3,39.1,0.07,0.07,ci_u(1),ci_u(2),360*pi/180,170*pi/180,2*pi,2*pi,46*3e-6,32*3e-6,20*pi/180];


%%
%GA算法参数
loop_step=1e4;%进化代数
N_chr=ncores*6;  %群体大小，总共多少组参数
N_gene=m; %基因数量（每组多少元素）
P_var=0.01;% 变异概率，0.1-0.001之间
N_var=ceil(N_chr*N_gene*P_var);
P_copy=0.03;%复制概率
N_copy=ceil(P_copy*N_chr);
if mod(N_copy,2)==1
    N_copy=N_copy+1;
end
N_diff=N_chr-N_copy;

%要保证N_chr和N_copy都是偶数



%%
%预分配空间
re_chi=zeros(1,N_chr);
ini_tr=zeros(500,Nplanet);
ini_num=zeros(1,Nplanet);
min_chi_array=zeros(1,loop_step);
% Fitness=zeros(1,N_chr);
% prob=zeros(1,N_chr);
% p_cum=zeros(1,N_chr);
% min_re_chi=zeros(1,N_chr);
% min_index=zeros(1,N_chr);
% temp_para=zeros(N_chr,N_gene);
% loca=zeros(1,N_chr/2);
% rand_order=zeros(1,1);
% tr=zeros(200,Nplanet);
% num=zeros(1,Nplanet);



%%
% GA算法开始
temp_var1=rand(N_chr,N_gene);
para_array=range_lower+(range_upper-range_lower).*temp_var1;

tic


for loop_out=1:loop_step
    
    parfor loop_1=1:N_chr
%                 [tr,num]=Nbody(para_array(loop_1,:),Nplanet,mc,Rs,err,a);
        [tr,num]=Nbody_forC(ini_tr,ini_num,para_array(loop_1,:),Nplanet,mc,Rs,err,a);
        re_chi(loop_1)=(sum((tr(planet_1_lx(:,1),1)-planet_1_lx(:,2)).^2./(planet_1_lx(:,3)).^2)+sum((tr(planet_2_lx(:,1),2)-planet_2_lx(:,2)).^2./(planet_2_lx(:,3)).^2))/(n-m);
    end
    
    %     [min_re_chi,min_index]=min(re_chi);
    %     min_chi_array(loop_out)=min_re_chi;
    %     loop_out,min_re_chi
    %     if abs(min_re_chi-1)<1e-3
    %         break;
    %     end
    
    [min_re_chi,min_index]=sort(re_chi);
    min_chi_array(loop_out)=min_re_chi(1);
    loop_out,min_re_chi(1)
    %         for loop_5=1:N_copy
    %             min_para(loop_5,:)=para_array(min_index(loop_5),:);
    %         end
    min_para=para_array(min_index(1:N_copy),:);
    
    if abs(min_re_chi(1)-1)<1e-3
        break;
    end
    
    
    
    
    %选择
    Fitness=1+max(re_chi)-re_chi; %适应度,目标 re_chi接近于1
    %         Fitness=1.0001*max(re_chi)-re_chi; %适应度,目标 re_chi接近于1
    prob=Fitness./sum(Fitness); %概率
    p_cum=cumsum(prob);%累积概率
    temp_para=para_array;
    %     parfor loop_2=1:N_chr
    %         [temp1,temp_index]=sort([rand,p_cum]);%找出rand位于哪个数之间
    %         sel_in=temp_index(find(temp_index==1)+1)-1;
    %         para_array(loop_2,:)=temp_para(sel_in,:);
    %     end
    [temp1,temp_index]=sort([rand(1,N_diff);repmat(p_cum',1,N_diff)]);%只选择N_diff行，找出rand位于哪个数之间
    sel_in=temp_index(find(temp_index==1)+1)-1;%sel_in有N_diff个数，代表被选中的哪组参数
    para_array=temp_para(sel_in,:);
    
    %交叉,两两配对交叉
    loca=randi(N_gene-1,1,N_diff/2);%交换的位置
    rand_order=randperm(N_diff); %随机配对
    temp_para=para_array;
    for loop_3=1:N_diff/2
        para_array(rand_order(loop_3),1:loca(loop_3))=temp_para(rand_order(loop_3+N_diff/2),1:loca(loop_3));
        para_array(rand_order(loop_3+N_diff/2),1:loca(loop_3))=temp_para(rand_order(loop_3),1:loca(loop_3));
    end
    
    %突变
    for loop_4=1:N_var
        rand_col=randi(N_gene);
        para_array(randi(N_diff),rand_col)=range_lower(rand_col)+(range_upper(rand_col)-range_lower(rand_col))*rand;
    end
    %     rand_col=randi(N_gene,1,N_var);
    %     rand_row=randi(N_diff,1,N_var);
    %     temp_var=range_lower(rand_col)+(range_upper(rand_col)-range_lower(rand_col)).*rand(1,N_var);
    %     para_array((rand_col-1)*N_diff+rand_row)=temp_var;
    
    %复制
    para_array=[para_array;min_para]; %每次保留最优值在群体里
    para_array=para_array(randperm(N_chr),:); %随机排序
    
    if (mod(loop_out,3000)==0)  %%暂停休息12分钟，要不然可能会变得很慢
        save temp_file
        clear all
        delete(gcp('nocreate'))
        pause(12*60)
        clear
        clear all
        
        load temp_file
        parpool(ncores)
        parfor i=1:ncores
            [notfound,warnings]=loadlibrary('C:\works\Nbody\codegen\dll\rkf78_7Np\rkf78_7Np.dll','C:\works\Nbody\codegen\dll\rkf78_7Np\rkf78_7Np.h');
            [notfound,warnings]=loadlibrary('C:\works\Nbody\codegen\dll\rkf78_forC\rkf78_forC.dll','C:\works\Nbody\codegen\dll\rkf78_forC\rkf78_forC.h');
            [notfound,warnings]=loadlibrary('C:\works\Nbody\codegen\dll\NbodySubProg\NbodySubProg.dll','C:\works\Nbody\codegen\dll\NbodySubProg\NbodySubProg.h');
        end
        clc
    end
    
    
end



times=toc


% save 1e4_6_480_1008



