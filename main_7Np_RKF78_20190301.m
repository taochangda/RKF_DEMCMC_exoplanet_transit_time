% This is the matlab script for fitting parameters of planets, 
% using Differential Evolution (DE) algorithm,DE-MCMC(Markov Chain Monte Carlo) algorithm
% the method of integration is RKF78(Runge-Kutta-Fehlberg7(8) )
%Version 1.0
% Tao Changda, June, 2019, All rights reserved.


clear all
file_path='C:\works\Nbody\'; %path of files 文件所在路径
startup_file=[file_path,'for_startup20190301.mat']; %path of source data 根据实际情况
load(startup_file);

%% Phase one: Differential Evolution（DE）algorithm, finding many groups optimal fitting parameter
%DE算法，找出较优值
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
            
            %DE算法的终止条件，目前我并没有找到确切的终止条件，这个需要再进一步探索，目前我取值max_dist为0.05，或者0.001。这个不是固定的标准。
            %1)，如果一直计算到re_chi不能再减小然后再终止计算，比如max_dist很小(比如1e-12），且连续多次接受概率为0(比如连续100次为0，sum(pre_acc_p)==0)，那有可能最后的参数都收敛到同一个值，这些无法计算收敛分布
            %2)，但如果max_dist过大，那也有可能DE算法过早的中止，re_chi的值这时候可能还会比较大。
            break;
        end
        
        
        temp_N1=zeros(1,n_link); %DE算法中，这个参数没用到。
        temp_N2=zeros(1,n_link); %record steps of mutation 记录每次迭代的变异步骤数
        pre_para=p_para; %record the last group fitting parameter 这步是必须的，保存上次参数值
        
        % import half of groups random values after certain  iteration
        % step. I am not sure it is useful. what I wanted is to control the DE algorithm not to converge the local best
        % 引入外部随机数，这一步我不是很确定是否有用，我的本意是经过一定步骤后引入一半(n_link/2)的随机数，以免DE算法容易陷入局部最优解
        if (mod(loop_out,gamma_reset_count)==half_gamma_reset_count)
            pre_para=[p_para(1:n_link/2,:);range_lower+(range_upper-range_lower).*rand(n_link/2,m)];
        end
        
        %reset gamma to aviod a too small gamma value 重置gamma，为避免gamma过小
        if (mod(loop_out,gamma_reset_count)==0)
            save([file_path,'DE_data20190301.mat'],'loop_out','re_chi','pre_acc_p','acc_p','gamma'); %为方便观察结果
            %             gamma=gamma_ini*(1+1e-4*randn);
            %             gamma=rand;
            gamma=1.0;
        end
        
        
        %--begin--step 1 of DE algorithm: MUTATIAON------DE算法 第1步：变异-------------------
        parfor loop_1=1:n_link % parallel compute. refer to matlab offical guide for more help   并行计算，需要调用大量分布式计算资源，计算快慢取决和计算资源和网络情况
            
            while (1) %确保在有效值域范围
                %R2,R3: random value in [1:n_link], and R2~=R3~=loop_1 在1：n_link中随机取值，且R2~=R3~=loop_1
                temp_var2=[1:loop_1-1,loop_1+1:n_link];
                temp_var2=temp_var2(randperm(n_link-1));
                %                 R1=temp_var2(1);
                R2=temp_var2(2);
                R3=temp_var2(3);
                %                 temp_para(loop_1,:)=pre_para(R1,:)+gamma*(pre_para(R2,:)-pre_para(R3,:));
                temp_para(loop_1,:)=pre_para(loop_1,:)+gamma*(pre_para(R2,:)-pre_para(R3,:));
                temp_N2(loop_1)=temp_N2(loop_1)+1;
                
                %all(temp_para(loop_1,:)-pre_para) 避免重复
                %                 if ((temp_para(loop_1,:)>=range_lower)&(temp_para(loop_1,:)<=range_upper)&all(temp_para(loop_1,:)-pre_para))
                if ((temp_para(loop_1,:)>=range_lower)&(temp_para(loop_1,:)<=range_upper)) %ensure valid value of the fitting parameter 要确保所有参数都在有效值范围内，这个是向量的比较
                    break;
                end
                if temp_N2(loop_1)>N_var %equal to last value if exceed  max steps 超过最大步数，直接等于原来的值
                    temp_para(loop_1,:)=pre_para(loop_1,:);
                    break;
                end
            end%确保在有效值域范围
            
        end
        %--end--step 1 of DE algorithm: MUTATIAON------DE算法 第1步：变异-------------------
        
        %--begin--step 2 of DE algorithm: CROSSOVER--DE算法 第2步：交叉---------------
        for loop_1=1:n_link
            CR_index=find(rand(1,m)<CR);
            temp_para(loop_1,CR_index)=pre_para(loop_1,CR_index);
        end
        %--end--step 2 of DE algorithm: CROSSOVER--DE算法 第2步：交叉---------------
        
        %----begin-step 3 of DE algorithm: SELECTION------DE算法 第3步：选择---------------------
        parfor loop_1=1:n_link %  parallel compute. refer to matlab offical guide for more help 并行计算，需要调用大量分布式计算资源，计算快慢取决和计算资源和网络情况
            [tr,error_flag]=Nbody_C_7Np_RKF78(temp_para(loop_1,:),mc,Rs,tstop,time_precise,err);
            chi_i_post(loop_1)=sum(((tr(planet_data{1}(:,1)+1,1)-planet_data{1}(:,2))./planet_data{1}(:,3)).^2)+sum(((tr(planet_data{2}(:,1)+1,2)-planet_data{2}(:,2))./planet_data{2}(:,3)).^2)+sum(((tr(planet_data{3}(:,1)+1,3)-planet_data{3}(:,2))./planet_data{3}(:,3)).^2)+sum(((tr(planet_data{4}(:,1)+1,4)-planet_data{4}(:,2))./planet_data{4}(:,3)).^2)+sum(((tr(planet_data{5}(:,1)+1,5)-planet_data{5}(:,2))./planet_data{5}(:,3)).^2)+sum(((tr(planet_data{6}(:,1)+1,6)-planet_data{6}(:,2))./planet_data{6}(:,3)).^2)+sum(((tr(planet_data{7}(:,1)+1,7)-planet_data{7}(:,2))./planet_data{7}(:,3)).^2);
            %error_flag is not used. because if error_flag=1, many zeros value will appear in "tr", and chi_i_post will be too big,
            %and this p_para will be abandon in DE algorithm or DE-MCMC algorithm
            %error_flag的程序中并未被使用，实际情况，对于error_flag返回值为1的情况，凌星时刻返回值tr会出现很多为0的值
            %，此时的chi_i_post会很大，那么在DE,DE-MCMC算法中，此组参数p_para会被舍弃。
        end
        
        
        %         chi_i_post=chi_i_post(find(temp_N2<=N_var));
        %         temp_para=temp_para(find(temp_N2<=N_var),:);
        chi_i_post=chi_i_post(temp_N2<=N_var);
        temp_para=temp_para(temp_N2<=N_var,:);
        
        
        %----begin--delele repeat parameter and reserve the optimal parameter at the same time-----去掉重复的参数，因为我们目的是尽量的保留取值的多样性，但保留最优参数不要被去掉----------
        [chi_i_post,index]=sort(chi_i_post);
        temp_para=temp_para(index,:);
        for j=1:m
            [temp_var1,row_index]=unique(temp_para(:,j),'stable');
            temp_para=temp_para(row_index,:);
            chi_i_post=chi_i_post(row_index);
        end
         %----end--delele repeat parameter and reserve the optimal parameter at the same time-----去掉重复的参数，因为我们目的是尽量的保留取值的多样性，但保留最优参数不要被去掉----------
        
        %---select the optimal n_link group in this and last fitting parameter ---从本次和上次的所有参数中选最优的n_link组，chi_i越小就越优
        temp_chi=[chi_i,chi_i_post];
        temp_para=[p_para;temp_para];
        [temp_minchi,index]=sort(temp_chi);
        p_para=temp_para(index(1:n_link),:);
        chi_i=temp_chi(index(1:n_link));
                
        %----begin-step 3 of DE algorithm: SELECTION------DE算法 第3步：选择---------------------
        
        
        
        %--begin--update gamma----根据接受概率、目标概率  来更新下一步gamma值----
        pre_acc_p=acc_p;%前一个值
        acc_p=sum(index(1:n_link)>n_link)/n_link; %accept probabillity 计算接受概率
        %         pre_acc_p=[pre_acc_p(2:gamma_reset_count),acc_p]; %最近gamma_reset_count的acc_p值
        
        if  (acc_p<(aim_acc_p*lower_acc_ratio))
            gamma=sqrt(lower_acc_ratio)*gamma;
        else
            if (acc_p>(aim_acc_p*upper_acc_ratio))
                gamma=sqrt(upper_acc_ratio)*gamma;
            else
                gamma=gamma*sqrt(acc_p/aim_acc_p);
            end
        end
        %--end--update gamma----根据接受概率、目标概率  来更新下一步gamma值----
        
        
        loop_out=loop_out+1;
        
        % reboot machine to adviod running more and more slowly 
        %-----因为matlab运行时间越长，就会越慢，所以到达一定运行步骤（步骤数根据实际情况定）会重启一下机器，重启前对数据进行保存，
        %另外操作系统要作相应的配置，确保重启后（在不登录的情况下）还能继续在后台继续运行本脚本。
        if (mod(loop_out,reboot_count)==0)  %%重启一把
            save(startup_file,'-v7.3');%save datas before reboot 这里的保存变量是为了重启后的载入继续运行
            pause(10);%pause a little time because saving data costs times 因为保存变量需要时间，所以暂时一段时间，如果变量过大，保存时间过长，需要调大此时间
            delete(gcp('nocreate')); %必须要释放并行任务，以免重启后并行任务冲突
            dos('shutdown /r'); %reboot windows OS 
            pause(63);%pause a little time because rebooting OS  does not execute immediately 这个是必须的，执行上面命令后约1分钟自动重启，这里加上63秒，如果不加，程序接着运行，可能会引起冲突
        end
    end %ending DE algorithm
    
    save([file_path,'DE_data20190301.mat'],'-v7.3');
    loop_out
    disp('begin burn-in step')
    
    %pre-precess for next phase(burn-in) 下一阶段（burn-in)的预处理
    DE_phase=0;
    
    %----sort---排序-------------------------------------
    [min_chi,min_index]=sort(chi_i);
    p_para=p_para(min_index,:);
    chi_i=chi_i(min_index);
    %--delele repeat parameter and reserve the optimal parameter at the same time--去重 使所有元素都不一样。
    temp_para=p_para;
    temp_chi=chi_i;
    source_index=1:n_link;
    for j=1:m
        [temp_var1,row_index]=unique(temp_para(:,j),'stable');
        temp_para=temp_para(row_index,:);
        temp_chi=temp_chi(row_index);
        source_index=source_index(row_index);
    end
    
    
    % select conve_N groups fitting parameter--选取conve_N组参数，
    %conve_N值设置多少我并不是很确定，但好像太大了不好收敛。有一种方法，可以在DE算法结束后，分别取conve_N为10，20，30，40。。用多台机器分别计算，看如一组收敛
    temp_length=length(temp_chi);
    if (temp_length>=conve_N)
        p_para=p_para(source_index(1:conve_N),:);
        chi_i=chi_i(source_index(1:conve_N));
    else  % 防止最后去重后不够conve_N组
        diff_index=setdiff([1:n_link],source_index);
        p_para=[temp_para;p_para(diff_index(1:(conve_N-temp_length)),:)];
        chi_i=[temp_chi,chi_i(diff_index(1:(conve_N-temp_length)))];
    end
    %-------------------------
    
    
    n_link=conve_N;%number of MCMC links 选取conve_N条马尔可夫链，
    temp_para=zeros(n_link,m);
    chi_i_post=zeros(1,n_link);
    cooling=1; %这个可能比较关键，直接设置成1
    gamma=gamma_ini*(1+1e-4*randn);
    loop_out=1;
    N_burn=1*1e5;%step of burn-in phase 
    reboot_count=1.5*1e4; %update reboot_count for next phase, pls based on the facts 
    ncores=n_link;
    
    save(startup_file,'-v7.3');
    
end

%% Phase two: burn-in phase, using DE-MCMC for iterative computing N_burn step. convergence is not computed in this phase
%burn-in 阶段，使用DE-MCMC算法，先算N_burn步，这个阶段并不计算收敛结果，终止条件为达到最大迭代步数(N_burn)
% DE-MCMC算法跟上面的DE算法的主要差别在于拟合参数的选择方法上，DE算法每次都是按最优值选，而DE-MCMC算法根据一定规则按某种概率选

if (burn_phase==1)
    
    
    delete(gcp('nocreate'));
    pool=parpool('local', ncores);%notice: call local machine cores（pls ensure that local cores>=conve_N）
    %这个是调用的是本机器核数，确本机器核数不小于 ncores(=conve_N)，本机并行计算，一般来说比较稳定，网络也不会成为瓶颈
    addAttachedFiles(pool,[file_path,'Nbody_C_7Np_RKF78.mexw64']);
    
    
    while(1)
        temp_N1=zeros(1,n_link);%计算接受概率的时候用到此变量
        temp_N2=zeros(1,n_link);
        pre_para=p_para;
        %     gamma=gamma_ini*(1+1e-4*randn);
        %         cooling=N_burn/loop_out; %cooling用上面计算的最后一个值 =1
        
        
        parfor loop_1=1:n_link %并行处理
            
            while (1) %确保在有效值域范围
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
            end%确保在有效值域范围
        end
        
        parfor loop_1=1:n_link %并行计算
            [tr,error_flag]=Nbody_C_7Np_RKF78(temp_para(loop_1,:),mc,Rs,tstop,time_precise,err);
            chi_i_post(loop_1)=sum(((tr(planet_data{1}(:,1)+1,1)-planet_data{1}(:,2))./planet_data{1}(:,3)).^2)+sum(((tr(planet_data{2}(:,1)+1,2)-planet_data{2}(:,2))./planet_data{2}(:,3)).^2)+sum(((tr(planet_data{3}(:,1)+1,3)-planet_data{3}(:,2))./planet_data{3}(:,3)).^2)+sum(((tr(planet_data{4}(:,1)+1,4)-planet_data{4}(:,2))./planet_data{4}(:,3)).^2)+sum(((tr(planet_data{5}(:,1)+1,5)-planet_data{5}(:,2))./planet_data{5}(:,3)).^2)+sum(((tr(planet_data{6}(:,1)+1,6)-planet_data{6}(:,2))./planet_data{6}(:,3)).^2)+sum(((tr(planet_data{7}(:,1)+1,7)-planet_data{7}(:,2))./planet_data{7}(:,3)).^2);
        end
        
        
        %-----selection of DE-MCMC algorithm, which is different from DE algorithm-------
        %-----此处为DE-MCMC算法的选择拟合参数的方法，与DE算法不一样-------
        d=min(exp(0.5*(chi_i-chi_i_post)/cooling),1);
        index=(d>=rand(1,n_link))&(temp_N2<=N_var);
        p_para(index,:)=temp_para(index,:);
        chi_i(index)=chi_i_post(index);
        temp_N1(index)=1;
        
        
        re_chi=chi_i./(n-m); %reduced chi平方
        
        
        acc_p=sum(temp_N1)/n_link;
        
        %---update gamma of DE-MCMC---DE-MCMC算法，目标接受概率直接设置为0.8，lower_acc_ratio约为0.8; upper_acc_ratio约为1.2;
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
        
        if (loop_out>N_burn) %达到最大迭代步数后终止burn-in 阶段
            break;
        end
        
        % reboot machine to adviod running more and more slowly 
        %-----因为matlab运行时间越长，就会越慢，所以到达一定运行步骤（步骤数根据实际情况定）会重启一下机器，重启前对数据进行保存，
        %另外操作系统要作相应的配置，确保重启后（在不登录的情况下）还能继续在后台继续运行本脚本。
        if (mod(loop_out,reboot_count)==0)  %%重启一把
            save(startup_file,'-v7.3'); %这里的保存变量是为了重启后的载入继续运行
            pause(10);%因为保存变量需要时间，所以暂时一段时间，如果变量过大，保存时间过长，需要调大此时间
            delete(gcp('nocreate')); %必须要释放并行任务，以免重启后并行任务冲突
            dos('shutdown /r');
            pause(63);%这个是必须的，执行上面命令后约1分钟自动重启，这里加上63秒，如果不加，程序接着运行，可能会引起冲突
            
        end
        
    end
    save([file_path,'burn_data20190301.mat'],'-v7.3');
    clear para_all_burn
    loop_out
    disp('begin calculate convergence')
    
    %pre-precess for next phase 下一阶段预处理
    burn_phase=0;
    loop_out=1;
    gamma=gamma_ini*(1+1e-4*randn);
    
    para_all=zeros(n_link,m,5*1e5); %需要保留每一步的拟合参数值
    para_sum=zeros(n_link,m); %所有步数的拟合参数对应求和
    
    R_index=1;
    R=zeros(1e2,m);
    re_chi_index=1;
    re_chi_array=zeros(1e3,n_link);
    
    save(startup_file,'-v7.3');
    
end
%% Phase three: convergence phase, using DE-MCMC for iterative computing untill reaching convergence condition. 
% convergence is judged in this phase, which is defferent from phase two
% 收敛计算，与burn-in的差别在于需要进行收敛计算


delete(gcp('nocreate'));
pool=parpool('local', ncores);%notice: call local machine cores（pls ensure that local cores>=conve_N）
addAttachedFiles(pool,[file_path,'Nbody_C_7Np_RKF78.mexw64']);


while(1)
    temp_N1=zeros(1,n_link);
    temp_N2=zeros(1,n_link);
    pre_para=p_para;
    %     gamma=gamma_ini*(1+1e-4*randn);
    %         cooling=N_burn/loop_out; %cooling用上面计算的最后一个值 =1
    
    
    parfor loop_1=1:n_link
        
        while (1) %确保在有效值域范围
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
        end%确保在有效值域范围
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
    
    
    %--begin--convergence compute and judge---收敛计算----
    para_all(:,:,loop_out)=p_para;  %record fitting parameter for all steps 记录每一步的拟合参数值
    para_sum=para_sum+p_para; %sum of parameter 求和
    
    if (mod(loop_out,5*1e4)==0) %每5*1e4步计算一次R值并判断是否收敛
        
        s_mean=sum(para_sum)/(loop_out*n_link);  %m个数;
        s_j_mean=para_sum/loop_out; %n_link*m个数
        B=loop_out/(n_link-1)*sum((s_j_mean-s_mean).^2);%m个数;
        s_t=sum((para_all(:,:,1:loop_out)-s_j_mean).^2,3);%n_link*m个数
        W=1/(n_link*(loop_out-1))*sum(s_t);%m个数;
        var=(loop_out-1)/loop_out*W+B/loop_out;%sigmma 平方
        
        R(R_index,:)= (n_link+1)/n_link*var./W-(loop_out-1)/(loop_out*n_link);
        if R(R_index,:)<1.1% convergence condition 所有参数R值均小于1.1则满足收敛条件，终止，计算
            break;
        end
        R_index=R_index+1;
    end
    %--end--convergence compute and judge---收敛计算----
    
    loop_out=loop_out+1;
    
    % reboot machine to adviod running more and more slowly 
    %-----因为matlab运行时间越长，就会越慢，所以到达一定运行步骤（步骤数根据实际情况定）会重启一下机器，重启前对数据进行保存，
    %另外操作系统要作相应的配置，确保重启后（在不登录的情况下）还能继续在后台继续运行本脚本。
    if (mod(loop_out,reboot_count)==0)
        %         time_stamp{re_chi_index}=datestr(now, 'yyyy-mm-dd HH:MM:SS');
        re_chi_array(re_chi_index,:)=re_chi;
        re_chi_index=re_chi_index+1;
        
        save(startup_file,'-v7.3');
        pause(20);%因为保存变量需要时间，所以暂时一段时间，如果变量过大，保存时间过长，需要调大此时间
        delete(gcp('nocreate'));
        dos('shutdown /r'); %reboot windows system.  just for windows
        pause(63);%这个是必须的，执行上面命令后约1分钟自动重启，这里加上63秒，如果不加，程序接着运行，可能会引起冲突
    end
    
    
    
end

conve_phase=0;
save([file_path,'FinalResult_RKF78_20190301.mat'],'-v7.3');

