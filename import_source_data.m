% Tao Changda, June, 2019, All rights reserved.

%% import observed data, set range of fitting parameter  导入凌星时刻数据，设置m个参数取值范围
clear all
loop_index=1;
n=0;

%range_lower and range_upper:  range of fitting parameters is suggestibly not set to be so wide。
%分别表示参数取值范围的下限和上限值，取值范围需要尽可能的小，不然可能计算量会很大，请根据具体情况而定。ci固定取值pi/2,
%omg固定取值为0，所以这两个参数不在这里做限制。质量mp乘以3e-6为的是将行星质量和恒星质量作统一
range_lower=[0*ones(1,7),0*ones(1,7),0*ones(1,7),0*ones(1,7),3e-6*0.1*ones(1,7)]; %[p1-p7,e1-e7,w1-w7,cm1-cm7,mp1-mp7],
range_upper=[0*ones(1,7),0.1*ones(1,7),2*pi*ones(1,7),2*pi*ones(1,7),3e-6*2*ones(1,7)]; %[[p1-p7,e1-e7,w1-w7,cm1-cm7,mp1-mp7]


source_file_path='C:\works\Nbody\'; %path of source data 数据源所在路径
data_file_path=[source_file_path,'correct20190531\7Np_data\284']; %path of source data  数据源所在路径
output_file_path=[source_file_path,'transit_data_7Np_XXX.mat']; %path of output file输出文件的路径以及文件名，可以根据实际情况更改

% % 重要：观测值起点t不等于0,需要减去一个常数delta_t,保证第0个观测点(如果有的话)时刻都小于其周期(p(1),p(2),...)
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
%planet_data为所有行星凌星时刻数据，planet_data{1}为第1个行星数据：为三维数组，第一列为序号，第二列为凌星时刻值，第三列为对应时刻误差。
%planet_data{2}，planet_data{3}，……，类似



% estimate the range of the period of the planet:
%下面程序为：大体估算p1-p7的上下边界，通过源数据，求出p的最大值 和最小值 ，最后再放大一下范围
%实现方法是找到源数据两两组之间的时间差除以序号差(delta_time/delta_No)，求最大和最小值
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
    
    range_lower(k)=mi-0.01; %0.01 is just estimate or empirical value, not   standard value. 0.01这里不是一个标准的值，这里仅是一个估计，根据实际情况更改
    range_upper(k)=ma+0.01; %0.01 is just estimate or empirical value, not   standard value.0.01这里不是一 个标准的值，这里仅是一个估计，根据实际情况更改
    
end

save(output_file_path,'planet_data','range_lower','range_upper','-v7.3');



