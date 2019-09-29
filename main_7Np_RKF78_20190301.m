% This is the matlab script for fitting parameters of planets, 
% using Differential Evolution (DE) algorithm,DE-MCMC(Markov Chain Monte Carlo) algorithm
% the method of integration is RKF78(Runge-Kutta-Fehlberg7(8) )
%Version 1.0
% Tao Changda, June, 2019, All rights reserved.


clear all
file_path='C:\works\Nbody\'; %path of files �ļ�����·��
startup_file=[file_path,'for_startup20190301.mat']; %path of source data ����ʵ�����
load(startup_file);

%% Phase one: Differential Evolution��DE��algorithm, finding many groups optimal fitting parameter
%DE�㷨���ҳ�����ֵ
if (DE_phase==1)
    
    delete(gcp('nocreate'));
    pool=parpool('myMJS', ncores);
    addAttachedFiles(pool,[file_path,'Nbody_C_7Np_RKF78.mexw64']);
    
    while(1)
        
        re_chi=chi_i./(n-m);
        %         if (loop_out>N_DE&&(max(re_chi)-min(re_chi))<max_dist&&sum(pre_acc_p)==0)
        %         if (loop_out>N_DE&&(max(re_chi)-min(re_chi))<max_dist&&(pre_acc_p+acc_p)==0)
        if (loop_out>N_DE&&(max(re_chi)-min(re_chi))<max_dist)
            % condition of ending DE algorithm. and I am not sure for this condition.  value of max_dist is just for the reference.
            % If value of max_dist is too small, all the group fitting parameters will converge to almost the same value.
            % If value of max_dist is big, DE algorithm will end early.
            
            %DE�㷨����ֹ������Ŀǰ�Ҳ�û���ҵ�ȷ�е���ֹ�����������Ҫ�ٽ�һ��̽����Ŀǰ��ȡֵmax_distΪ0.05������0.001��������ǹ̶��ı�׼��
            %1)�����һֱ���㵽re_chi�����ټ�СȻ������ֹ���㣬����max_dist��С(����1e-12������������ν��ܸ���Ϊ0(��������100��Ϊ0��sum(pre_acc_p)==0)�����п������Ĳ�����������ͬһ��ֵ����Щ�޷����������ֲ�
            %2)�������max_dist������Ҳ�п���DE�㷨�������ֹ��re_chi��ֵ��ʱ����ܻ���Ƚϴ�
            break;
        end
        
        
        temp_N1=zeros(1,n_link); %DE�㷨�У��������û�õ���
        temp_N2=zeros(1,n_link); %record steps of mutation ��¼ÿ�ε����ı��첽����
        pre_para=p_para; %record the last group fitting parameter �ⲽ�Ǳ���ģ������ϴβ���ֵ
        
        % import half of groups random values after certain  iteration
        % step. I am not sure it is useful. what I wanted is to control the DE algorithm not to converge the local best
        % �����ⲿ���������һ���Ҳ��Ǻ�ȷ���Ƿ����ã��ҵı����Ǿ���һ�����������һ��(n_link/2)�������������DE�㷨��������ֲ����Ž�
        if (mod(loop_out,gamma_reset_count)==half_gamma_reset_count)
            pre_para=[p_para(1:n_link/2,:);range_lower+(range_upper-range_lower).*rand(n_link/2,m)];
        end
        
        %reset gamma to aviod a too small gamma value ����gamma��Ϊ����gamma��С
        if (mod(loop_out,gamma_reset_count)==0)
            save([file_path,'DE_data20190301.mat'],'loop_out','re_chi','pre_acc_p','acc_p','gamma'); %Ϊ����۲���
            %             gamma=gamma_ini*(1+1e-4*randn);
            %             gamma=rand;
            gamma=1.0;
        end
        
        
        %--begin--step 1 of DE algorithm: MUTATIAON------DE�㷨 ��1��������-------------------
        parfor loop_1=1:n_link % parallel compute. refer to matlab offical guide for more help   ���м��㣬��Ҫ���ô����ֲ�ʽ������Դ���������ȡ���ͼ�����Դ���������
            
            while (1) %ȷ������Чֵ��Χ
                %R2,R3: random value in [1:n_link], and R2~=R3~=loop_1 ��1��n_link�����ȡֵ����R2~=R3~=loop_1
                temp_var2=[1:loop_1-1,loop_1+1:n_link];
                temp_var2=temp_var2(randperm(n_link-1));
                %                 R1=temp_var2(1);
                R2=temp_var2(2);
                R3=temp_var2(3);
                %                 temp_para(loop_1,:)=pre_para(R1,:)+gamma*(pre_para(R2,:)-pre_para(R3,:));
                temp_para(loop_1,:)=pre_para(loop_1,:)+gamma*(pre_para(R2,:)-pre_para(R3,:));
                temp_N2(loop_1)=temp_N2(loop_1)+1;
                
                %all(temp_para(loop_1,:)-pre_para) �����ظ�
                %                 if ((temp_para(loop_1,:)>=range_lower)&(temp_para(loop_1,:)<=range_upper)&all(temp_para(loop_1,:)-pre_para))
                if ((temp_para(loop_1,:)>=range_lower)&(temp_para(loop_1,:)<=range_upper)) %ensure valid value of the fitting parameter Ҫȷ�����в���������Чֵ��Χ�ڣ�����������ıȽ�
                    break;
                end
                if temp_N2(loop_1)>N_var %equal to last value if exceed  max steps �����������ֱ�ӵ���ԭ����ֵ
                    temp_para(loop_1,:)=pre_para(loop_1,:);
                    break;
                end
            end%ȷ������Чֵ��Χ
            
        end
        %--end--step 1 of DE algorithm: MUTATIAON------DE�㷨 ��1��������-------------------
        
        %--begin--step 2 of DE algorithm: CROSSOVER--DE�㷨 ��2��������---------------
        for loop_1=1:n_link
            CR_index=find(rand(1,m)<CR);
            temp_para(loop_1,CR_index)=pre_para(loop_1,CR_index);
        end
        %--end--step 2 of DE algorithm: CROSSOVER--DE�㷨 ��2��������---------------
        
        %----begin-step 3 of DE algorithm: SELECTION------DE�㷨 ��3����ѡ��---------------------
        parfor loop_1=1:n_link %  parallel compute. refer to matlab offical guide for more help ���м��㣬��Ҫ���ô����ֲ�ʽ������Դ���������ȡ���ͼ�����Դ���������
            [tr,error_flag]=Nbody_C_7Np_RKF78(temp_para(loop_1,:),mc,Rs,tstop,time_precise,err);
            chi_i_post(loop_1)=sum(((tr(planet_data{1}(:,1)+1,1)-planet_data{1}(:,2))./planet_data{1}(:,3)).^2)+sum(((tr(planet_data{2}(:,1)+1,2)-planet_data{2}(:,2))./planet_data{2}(:,3)).^2)+sum(((tr(planet_data{3}(:,1)+1,3)-planet_data{3}(:,2))./planet_data{3}(:,3)).^2)+sum(((tr(planet_data{4}(:,1)+1,4)-planet_data{4}(:,2))./planet_data{4}(:,3)).^2)+sum(((tr(planet_data{5}(:,1)+1,5)-planet_data{5}(:,2))./planet_data{5}(:,3)).^2)+sum(((tr(planet_data{6}(:,1)+1,6)-planet_data{6}(:,2))./planet_data{6}(:,3)).^2)+sum(((tr(planet_data{7}(:,1)+1,7)-planet_data{7}(:,2))./planet_data{7}(:,3)).^2);
            %error_flag is not used. because if error_flag=1, many zeros value will appear in "tr", and chi_i_post will be too big,
            %and this p_para will be abandon in DE algorithm or DE-MCMC algorithm
            %error_flag�ĳ����в�δ��ʹ�ã�ʵ�����������error_flag����ֵΪ1�����������ʱ�̷���ֵtr����ֺܶ�Ϊ0��ֵ
            %����ʱ��chi_i_post��ܴ���ô��DE,DE-MCMC�㷨�У��������p_para�ᱻ������
        end
        
        
        %         chi_i_post=chi_i_post(find(temp_N2<=N_var));
        %         temp_para=temp_para(find(temp_N2<=N_var),:);
        chi_i_post=chi_i_post(temp_N2<=N_var);
        temp_para=temp_para(temp_N2<=N_var,:);
        
        
        %----begin--delele repeat parameter and reserve the optimal parameter at the same time-----ȥ���ظ��Ĳ�������Ϊ����Ŀ���Ǿ����ı���ȡֵ�Ķ����ԣ����������Ų�����Ҫ��ȥ��----------
        [chi_i_post,index]=sort(chi_i_post);
        temp_para=temp_para(index,:);
        for j=1:m
            [temp_var1,row_index]=unique(temp_para(:,j),'stable');
            temp_para=temp_para(row_index,:);
            chi_i_post=chi_i_post(row_index);
        end
         %----end--delele repeat parameter and reserve the optimal parameter at the same time-----ȥ���ظ��Ĳ�������Ϊ����Ŀ���Ǿ����ı���ȡֵ�Ķ����ԣ����������Ų�����Ҫ��ȥ��----------
        
        %---select the optimal n_link group in this and last fitting parameter ---�ӱ��κ��ϴε����в�����ѡ���ŵ�n_link�飬chi_iԽС��Խ��
        temp_chi=[chi_i,chi_i_post];
        temp_para=[p_para;temp_para];
        [temp_minchi,index]=sort(temp_chi);
        p_para=temp_para(index(1:n_link),:);
        chi_i=temp_chi(index(1:n_link));
                
        %----begin-step 3 of DE algorithm: SELECTION------DE�㷨 ��3����ѡ��---------------------
        
        
        
        %--begin--update gamma----���ݽ��ܸ��ʡ�Ŀ�����  ��������һ��gammaֵ----
        pre_acc_p=acc_p;%ǰһ��ֵ
        acc_p=sum(index(1:n_link)>n_link)/n_link; %accept probabillity ������ܸ���
        %         pre_acc_p=[pre_acc_p(2:gamma_reset_count),acc_p]; %���gamma_reset_count��acc_pֵ
        
        if  (acc_p<(aim_acc_p*lower_acc_ratio))
            gamma=sqrt(lower_acc_ratio)*gamma;
        else
            if (acc_p>(aim_acc_p*upper_acc_ratio))
                gamma=sqrt(upper_acc_ratio)*gamma;
            else
                gamma=gamma*sqrt(acc_p/aim_acc_p);
            end
        end
        %--end--update gamma----���ݽ��ܸ��ʡ�Ŀ�����  ��������һ��gammaֵ----
        
        
        loop_out=loop_out+1;
        
        % reboot machine to adviod running more and more slowly 
        %-----��Ϊmatlab����ʱ��Խ�����ͻ�Խ�������Ե���һ�����в��裨����������ʵ���������������һ�»���������ǰ�����ݽ��б��棬
        %�������ϵͳҪ����Ӧ�����ã�ȷ���������ڲ���¼������£����ܼ����ں�̨�������б��ű���
        if (mod(loop_out,reboot_count)==0)  %%����һ��
            save(startup_file,'-v7.3');%save datas before reboot ����ı��������Ϊ��������������������
            pause(10);%pause a little time because saving data costs times ��Ϊ���������Ҫʱ�䣬������ʱһ��ʱ�䣬����������󣬱���ʱ���������Ҫ�����ʱ��
            delete(gcp('nocreate')); %����Ҫ�ͷŲ������������������������ͻ
            dos('shutdown /r'); %reboot windows OS 
            pause(63);%pause a little time because rebooting OS  does not execute immediately ����Ǳ���ģ�ִ�����������Լ1�����Զ��������������63�룬������ӣ�����������У����ܻ������ͻ
        end
    end %ending DE algorithm
    
    save([file_path,'DE_data20190301.mat'],'-v7.3');
    loop_out
    disp('begin burn-in step')
    
    %pre-precess for next phase(burn-in) ��һ�׶Σ�burn-in)��Ԥ����
    DE_phase=0;
    
    %----sort---����-------------------------------------
    [min_chi,min_index]=sort(chi_i);
    p_para=p_para(min_index,:);
    chi_i=chi_i(min_index);
    %--delele repeat parameter and reserve the optimal parameter at the same time--ȥ�� ʹ����Ԫ�ض���һ����
    temp_para=p_para;
    temp_chi=chi_i;
    source_index=1:n_link;
    for j=1:m
        [temp_var1,row_index]=unique(temp_para(:,j),'stable');
        temp_para=temp_para(row_index,:);
        temp_chi=temp_chi(row_index);
        source_index=source_index(row_index);
    end
    
    
    % select conve_N groups fitting parameter--ѡȡconve_N�������
    %conve_Nֵ���ö����Ҳ����Ǻ�ȷ����������̫���˲�����������һ�ַ�����������DE�㷨�����󣬷ֱ�ȡconve_NΪ10��20��30��40�����ö�̨�����ֱ���㣬����һ������
    temp_length=length(temp_chi);
    if (temp_length>=conve_N)
        p_para=p_para(source_index(1:conve_N),:);
        chi_i=chi_i(source_index(1:conve_N));
    else  % ��ֹ���ȥ�غ󲻹�conve_N��
        diff_index=setdiff([1:n_link],source_index);
        p_para=[temp_para;p_para(diff_index(1:(conve_N-temp_length)),:)];
        chi_i=[temp_chi,chi_i(diff_index(1:(conve_N-temp_length)))];
    end
    %-------------------------
    
    
    n_link=conve_N;%number of MCMC links ѡȡconve_N������ɷ�����
    temp_para=zeros(n_link,m);
    chi_i_post=zeros(1,n_link);
    cooling=1; %������ܱȽϹؼ���ֱ�����ó�1
    gamma=gamma_ini*(1+1e-4*randn);
    loop_out=1;
    N_burn=1*1e5;%step of burn-in phase 
    reboot_count=1.5*1e4; %update reboot_count for next phase, pls based on the facts 
    ncores=n_link;
    
    save(startup_file,'-v7.3');
    
end

%% Phase two: burn-in phase, using DE-MCMC for iterative computing N_burn step. convergence is not computed in this phase
%burn-in �׶Σ�ʹ��DE-MCMC�㷨������N_burn��������׶β������������������ֹ����Ϊ�ﵽ����������(N_burn)
% DE-MCMC�㷨�������DE�㷨����Ҫ���������ϲ�����ѡ�񷽷��ϣ�DE�㷨ÿ�ζ��ǰ�����ֵѡ����DE-MCMC�㷨����һ������ĳ�ָ���ѡ

if (burn_phase==1)
    
    
    delete(gcp('nocreate'));
    pool=parpool('local', ncores);%notice: call local machine cores��pls ensure that local cores>=conve_N��
    %����ǵ��õ��Ǳ�����������ȷ������������С�� ncores(=conve_N)���������м��㣬һ����˵�Ƚ��ȶ�������Ҳ�����Ϊƿ��
    addAttachedFiles(pool,[file_path,'Nbody_C_7Np_RKF78.mexw64']);
    
    
    while(1)
        temp_N1=zeros(1,n_link);%������ܸ��ʵ�ʱ���õ��˱���
        temp_N2=zeros(1,n_link);
        pre_para=p_para;
        %     gamma=gamma_ini*(1+1e-4*randn);
        %         cooling=N_burn/loop_out; %cooling�������������һ��ֵ =1
        
        
        parfor loop_1=1:n_link %���д���
            
            while (1) %ȷ������Чֵ��Χ
                temp_var2=[1:loop_1-1,loop_1+1:n_link];
                temp_var2=temp_var2(randperm(n_link-1));
                %                 R1=temp_var2(1);
                R2=temp_var2(2);
                R3=temp_var2(3);
                %                         b=range_level.*randn(1,m);
                %                         temp_para(loop_1,:)=pre_para(loop_1,:)+gamma*(pre_para(R1,:)-pre_para(R2,:))+(-b+2*b.*rand(1,m));
                temp_para(loop_1,:)=pre_para(loop_1,:)+gamma*(pre_para(R2,:)-pre_para(R3,:));
                
                if ((temp_para(loop_1,:)>=range_lower)&(temp_para(loop_1,:)<=range_upper))
                    break;
                end
                temp_N2(loop_1)=temp_N2(loop_1)+1;
                if temp_N2(loop_1)>N_var
                    %                                 temp_para(loop_1,:)=range_lower+(range_upper-range_lower).*rand(1,m);
                    temp_para(loop_1,:)=pre_para(loop_1,:);
                    break;
                end
            end%ȷ������Чֵ��Χ
        end
        
        parfor loop_1=1:n_link %���м���
            [tr,error_flag]=Nbody_C_7Np_RKF78(temp_para(loop_1,:),mc,Rs,tstop,time_precise,err);
            chi_i_post(loop_1)=sum(((tr(planet_data{1}(:,1)+1,1)-planet_data{1}(:,2))./planet_data{1}(:,3)).^2)+sum(((tr(planet_data{2}(:,1)+1,2)-planet_data{2}(:,2))./planet_data{2}(:,3)).^2)+sum(((tr(planet_data{3}(:,1)+1,3)-planet_data{3}(:,2))./planet_data{3}(:,3)).^2)+sum(((tr(planet_data{4}(:,1)+1,4)-planet_data{4}(:,2))./planet_data{4}(:,3)).^2)+sum(((tr(planet_data{5}(:,1)+1,5)-planet_data{5}(:,2))./planet_data{5}(:,3)).^2)+sum(((tr(planet_data{6}(:,1)+1,6)-planet_data{6}(:,2))./planet_data{6}(:,3)).^2)+sum(((tr(planet_data{7}(:,1)+1,7)-planet_data{7}(:,2))./planet_data{7}(:,3)).^2);
        end
        
        
        %-----selection of DE-MCMC algorithm, which is different from DE algorithm-------
        %-----�˴�ΪDE-MCMC�㷨��ѡ����ϲ����ķ�������DE�㷨��һ��-------
        d=min(exp(0.5*(chi_i-chi_i_post)/cooling),1);
        index=(d>=rand(1,n_link))&(temp_N2<=N_var);
        p_para(index,:)=temp_para(index,:);
        chi_i(index)=chi_i_post(index);
        temp_N1(index)=1;
        
        
        re_chi=chi_i./(n-m); %reduced chiƽ��
        
        
        acc_p=sum(temp_N1)/n_link;
        
        %---update gamma of DE-MCMC---DE-MCMC�㷨��Ŀ����ܸ���ֱ������Ϊ0.8��lower_acc_ratioԼΪ0.8; upper_acc_ratioԼΪ1.2;
        %     if acc_p>0
        if  (acc_p<0.2)
            gamma=0.8944*gamma;
        else
            if (acc_p>0.31)
                gamma=1.1136*gamma;
            else
                gamma=gamma*sqrt(acc_p/0.25);
            end
        end
        %     end
        
        if (mod(loop_out,gamma_reset_count)==0)
            %             gamma=gamma_ini*(1+1e-4*randn);
            gamma=1.0;
            %         gamma_ini=1;
            
            %             loop_out,min(re_chi),max(re_chi),acc_p
            
        end
        
        %--temp---------
        %         para_all_burn(:,:,loop_out)=p_para;
        %----------
        
        loop_out=loop_out+1;
        
        if (loop_out>N_burn) %�ﵽ��������������ֹburn-in �׶�
            break;
        end
        
        % reboot machine to adviod running more and more slowly 
        %-----��Ϊmatlab����ʱ��Խ�����ͻ�Խ�������Ե���һ�����в��裨����������ʵ���������������һ�»���������ǰ�����ݽ��б��棬
        %�������ϵͳҪ����Ӧ�����ã�ȷ���������ڲ���¼������£����ܼ����ں�̨�������б��ű���
        if (mod(loop_out,reboot_count)==0)  %%����һ��
            save(startup_file,'-v7.3'); %����ı��������Ϊ��������������������
            pause(10);%��Ϊ���������Ҫʱ�䣬������ʱһ��ʱ�䣬����������󣬱���ʱ���������Ҫ�����ʱ��
            delete(gcp('nocreate')); %����Ҫ�ͷŲ������������������������ͻ
            dos('shutdown /r');
            pause(63);%����Ǳ���ģ�ִ�����������Լ1�����Զ��������������63�룬������ӣ�����������У����ܻ������ͻ
            
        end
        
    end
    save([file_path,'burn_data20190301.mat'],'-v7.3');
    clear para_all_burn
    loop_out
    disp('begin calculate convergence')
    
    %pre-precess for next phase ��һ�׶�Ԥ����
    burn_phase=0;
    loop_out=1;
    gamma=gamma_ini*(1+1e-4*randn);
    
    para_all=zeros(n_link,m,5*1e5); %��Ҫ����ÿһ������ϲ���ֵ
    para_sum=zeros(n_link,m); %���в�������ϲ�����Ӧ���
    
    R_index=1;
    R=zeros(1e2,m);
    re_chi_index=1;
    re_chi_array=zeros(1e3,n_link);
    
    save(startup_file,'-v7.3');
    
end
%% Phase three: convergence phase, using DE-MCMC for iterative computing untill reaching convergence condition. 
% convergence is judged in this phase, which is defferent from phase two
% �������㣬��burn-in�Ĳ��������Ҫ������������


delete(gcp('nocreate'));
pool=parpool('local', ncores);%notice: call local machine cores��pls ensure that local cores>=conve_N��
addAttachedFiles(pool,[file_path,'Nbody_C_7Np_RKF78.mexw64']);


while(1)
    temp_N1=zeros(1,n_link);
    temp_N2=zeros(1,n_link);
    pre_para=p_para;
    %     gamma=gamma_ini*(1+1e-4*randn);
    %         cooling=N_burn/loop_out; %cooling�������������һ��ֵ =1
    
    
    parfor loop_1=1:n_link
        
        while (1) %ȷ������Чֵ��Χ
            temp_var2=[1:loop_1-1,loop_1+1:n_link];
            temp_var2=temp_var2(randperm(n_link-1));
            %                 R1=temp_var2(1);
            R2=temp_var2(2);
            R3=temp_var2(3);
            %                         b=range_level.*randn(1,m);
            %                         temp_para(loop_1,:)=pre_para(loop_1,:)+gamma*(pre_para(R1,:)-pre_para(R2,:))+(-b+2*b.*rand(1,m));
            temp_para(loop_1,:)=pre_para(loop_1,:)+gamma*(pre_para(R2,:)-pre_para(R3,:));
            
            if ((temp_para(loop_1,:)>=range_lower)&(temp_para(loop_1,:)<=range_upper))
                break;
            end
            temp_N2(loop_1)=temp_N2(loop_1)+1;
            if temp_N2(loop_1)>N_var
                %                                 temp_para(loop_1,:)=range_lower+(range_upper-range_lower).*rand(1,m);
                temp_para(loop_1,:)=pre_para(loop_1,:);
                break;
            end
        end%ȷ������Чֵ��Χ
    end
    
    parfor loop_1=1:n_link
        [tr,error_flag]=Nbody_C_7Np_RKF78(temp_para(loop_1,:),mc,Rs,tstop,time_precise,err);
        chi_i_post(loop_1)=sum(((tr(planet_data{1}(:,1)+1,1)-planet_data{1}(:,2))./planet_data{1}(:,3)).^2)+sum(((tr(planet_data{2}(:,1)+1,2)-planet_data{2}(:,2))./planet_data{2}(:,3)).^2)+sum(((tr(planet_data{3}(:,1)+1,3)-planet_data{3}(:,2))./planet_data{3}(:,3)).^2)+sum(((tr(planet_data{4}(:,1)+1,4)-planet_data{4}(:,2))./planet_data{4}(:,3)).^2)+sum(((tr(planet_data{5}(:,1)+1,5)-planet_data{5}(:,2))./planet_data{5}(:,3)).^2)+sum(((tr(planet_data{6}(:,1)+1,6)-planet_data{6}(:,2))./planet_data{6}(:,3)).^2)+sum(((tr(planet_data{7}(:,1)+1,7)-planet_data{7}(:,2))./planet_data{7}(:,3)).^2);
    end
    
    %         for loop_1=1:n_link
    %             d=min(exp(0.5*(chi_i(loop_1)-chi_i_post(loop_1))/cooling),1);
    %             if (d>=rand&&temp_N2(loop_1)<=N_var)
    %                 p_para(loop_1,:)=temp_para(loop_1,:);
    %                 chi_i(loop_1)=chi_i_post(loop_1);
    %                 temp_N1(loop_1)=1;
    %             end
    %         end
    
    
    d=min(exp(0.5*(chi_i-chi_i_post)/cooling),1);
    index=(d>=rand(1,n_link))&(temp_N2<=N_var);
    p_para(index,:)=temp_para(index,:);
    chi_i(index)=chi_i_post(index);
    temp_N1(index)=1;
    
    
    
    
    re_chi=chi_i./(n-m);
    
    
    acc_p=sum(temp_N1)/n_link;
    %     if acc_p>0
    if  (acc_p<0.2)
        gamma=0.8944*gamma;
    else
        if (acc_p>0.31)
            gamma=1.1136*gamma;
        else
            gamma=gamma*sqrt(acc_p/0.25);
        end
    end
    %     end
    
    if (mod(loop_out,gamma_reset_count)==0)
        %         gamma=gamma_ini*(1+1e-4*randn);
        gamma=1.0;
        %             gamma=gamma_ini;
        %         gamma_ini=1;
        
    end
    
    
    %--begin--convergence compute and judge---��������----
    para_all(:,:,loop_out)=p_para;  %record fitting parameter for all steps ��¼ÿһ������ϲ���ֵ
    para_sum=para_sum+p_para; %sum of parameter ���
    
    if (mod(loop_out,5*1e4)==0) %ÿ5*1e4������һ��Rֵ���ж��Ƿ�����
        
        s_mean=sum(para_sum)/(loop_out*n_link);  %m����;
        s_j_mean=para_sum/loop_out; %n_link*m����
        B=loop_out/(n_link-1)*sum((s_j_mean-s_mean).^2);%m����;
        s_t=sum((para_all(:,:,1:loop_out)-s_j_mean).^2,3);%n_link*m����
        W=1/(n_link*(loop_out-1))*sum(s_t);%m����;
        var=(loop_out-1)/loop_out*W+B/loop_out;%sigmma ƽ��
        
        R(R_index,:)= (n_link+1)/n_link*var./W-(loop_out-1)/(loop_out*n_link);
        if R(R_index,:)<1.1% convergence condition ���в���Rֵ��С��1.1������������������ֹ������
            break;
        end
        R_index=R_index+1;
    end
    %--end--convergence compute and judge---��������----
    
    loop_out=loop_out+1;
    
    % reboot machine to adviod running more and more slowly 
    %-----��Ϊmatlab����ʱ��Խ�����ͻ�Խ�������Ե���һ�����в��裨����������ʵ���������������һ�»���������ǰ�����ݽ��б��棬
    %�������ϵͳҪ����Ӧ�����ã�ȷ���������ڲ���¼������£����ܼ����ں�̨�������б��ű���
    if (mod(loop_out,reboot_count)==0)
        %         time_stamp{re_chi_index}=datestr(now, 'yyyy-mm-dd HH:MM:SS');
        re_chi_array(re_chi_index,:)=re_chi;
        re_chi_index=re_chi_index+1;
        
        save(startup_file,'-v7.3');
        pause(20);%��Ϊ���������Ҫʱ�䣬������ʱһ��ʱ�䣬����������󣬱���ʱ���������Ҫ�����ʱ��
        delete(gcp('nocreate'));
        dos('shutdown /r'); %reboot windows system.  just for windows
        pause(63);%����Ǳ���ģ�ִ�����������Լ1�����Զ��������������63�룬������ӣ�����������У����ܻ������ͻ
    end
    
    
    
end

conve_phase=0;
save([file_path,'FinalResult_RKF78_20190301.mat'],'-v7.3');

