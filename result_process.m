% Tao Changda, June, 2019, All rights reserved.

%% result process and plot histogram �������ת������ͼ
clear all
close all
clc
load('C:\works\Nbody\XXX\FinalResult_XXX.mat') %load  result data ���������ݣ�����ʵ���������

% transform 3-dimenstion array to 2-dimenstion array ����ά����תΪ��ά����
para_all=para_all(:,:,1:loop_out);
[dim_1,dim_2,dim_3]=size(para_all);
temp_para=reshape(permute(para_all,[2 1 3]),dim_2,dim_1*dim_3);
temp_para=temp_para';

% close all
% for j=1:m
%     figure
%     histogram(temp_para(:,j));
% end


%delete repeat data (Is it neccessary? I am not sure) ȥ�ظ���ȥ���صĲ���ȡֵ
for j=1:m
    [temp_var1,row_index]=unique(temp_para(:,j),'stable');
    temp_para=temp_para(row_index,:);
    %         chi_i_post=chi_i_post(row_index);
end

%plot histogram ��ֱ��ͼ
close all
for j=1:m
    figure
    histogram(temp_para(:,j));
end



