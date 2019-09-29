% Tao Changda, June, 2019, All rights reserved.

%% import observed data, set range of fitting parameter  ��������ʱ�����ݣ�����m������ȡֵ��Χ
clear all
loop_index=1;
n=0;

%range_lower and range_upper:  range of fitting parameters is suggestibly not set to be so wide��
%�ֱ��ʾ����ȡֵ��Χ�����޺�����ֵ��ȡֵ��Χ��Ҫ�����ܵ�С����Ȼ���ܼ�������ܴ�����ݾ������������ci�̶�ȡֵpi/2,
%omg�̶�ȡֵΪ0�����������������������������ơ�����mp����3e-6Ϊ���ǽ����������ͺ���������ͳһ
range_lower=[0*ones(1,7),0*ones(1,7),0*ones(1,7),0*ones(1,7),3e-6*0.1*ones(1,7)]; %[p1-p7,e1-e7,w1-w7,cm1-cm7,mp1-mp7],
range_upper=[0*ones(1,7),0.1*ones(1,7),2*pi*ones(1,7),2*pi*ones(1,7),3e-6*2*ones(1,7)]; %[[p1-p7,e1-e7,w1-w7,cm1-cm7,mp1-mp7]


source_file_path='C:\works\Nbody\'; %path of source data ����Դ����·��
data_file_path=[source_file_path,'correct20190531\7Np_data\284']; %path of source data  ����Դ����·��
output_file_path=[source_file_path,'transit_data_7Np_XXX.mat']; %path of output file����ļ���·���Լ��ļ��������Ը���ʵ���������

% % ��Ҫ���۲�ֵ���t������0,��Ҫ��ȥһ������delta_t,��֤��0���۲��(����еĻ�)ʱ�̶�С��������(p(1),p(2),...)
% time correct, if the start time of observed data is not "zeros", transit
% time minus constant value delta_t
delta_t=0;

% for i='b':'h'
for i='1':'7'
    filename=[data_file_path,'\t',i,'.txt'];
    [data1,data2,data3]=textread(filename,'%n%n%n');
    data2=data2-delta_t;
    planet_data{loop_index}=[data1,data2,data3];
    loop_index=loop_index+1;
    n=n+length(data1);
end
% planet_data: transit time and errors of all planets.
%planet_dataΪ������������ʱ�����ݣ�planet_data{1}Ϊ��1���������ݣ�Ϊ��ά���飬��һ��Ϊ��ţ��ڶ���Ϊ����ʱ��ֵ��������Ϊ��Ӧʱ����
%planet_data{2}��planet_data{3}������������



% estimate the range of the period of the planet:
%�������Ϊ���������p1-p7�����±߽磬ͨ��Դ���ݣ����p�����ֵ ����Сֵ ������ٷŴ�һ�·�Χ
%ʵ�ַ������ҵ�Դ����������֮���ʱ��������Ų�(delta_time/delta_No)����������Сֵ
for k=1:7
    data=planet_data{k};
    NN=size(data,1);
    ma=0;
    mi=1e10;
    for i=1:NN-1
        for j=i+1:NN
            temp_var=(data(j,2)-data(i,2))/(data(j,1)-data(i,1)); %(delta_time/delta_No)
            if (temp_var>ma)
                ma=temp_var;
            end
            
            if (temp_var<mi)
                mi=temp_var;
            end
            
        end
    end
    
    range_lower(k)=mi-0.01; %0.01 is just estimate or empirical value, not   standard value. 0.01���ﲻ��һ����׼��ֵ���������һ�����ƣ�����ʵ���������
    range_upper(k)=ma+0.01; %0.01 is just estimate or empirical value, not   standard value.0.01���ﲻ��һ ����׼��ֵ���������һ�����ƣ�����ʵ���������
    
end

save(output_file_path,'planet_data','range_lower','range_upper','-v7.3');



