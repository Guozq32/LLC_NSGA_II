function [new_pop,data_Para_new,data_Para_cell_new]=replacement(population_inter_sorted,data_Para_sorted,data_Para_cell_sorted, front)
global pop_size 
%% Description
% The next generation population is formed by appending each front subsequently until the
% population size exceeds the current population size. If When adding all the individuals
% of any front, the population exceeds the population size, then the required number of 
% remaining individuals alone are selected from that particular front based
% on crowding distance.
%% code starts
index=0;
ii=1;
while index < pop_size
    l_f=length(front(ii).fr);  %%l_f是[2pop_size，1]，记录支配级别
    if index+l_f < pop_size  %% 如果分层仍小于总种群数，则把此分层作为新种群
        new_pop(index+1:index+l_f,:)= population_inter_sorted(index+1:index+l_f,:);
        data_Para_new(index+1:index+l_f,:)= data_Para_sorted(index+1:index+l_f,:);
        data_Para_cell_new(index+1:index+l_f,:)=data_Para_cell_sorted(index+1:index+l_f,:);
        index=index+l_f;
    else
            temp1=population_inter_sorted(index+1:index+l_f,:);  % population_inter_sorted V+1列到V+M列为目标函数值; V+M+1列:zeros(pop_size1,1)表示误差为0  % V+M+1列对应层级；V+M+4列对应总拥挤度
            [temp2,sort]=sortrows(temp1,size(temp1,2)); %% 按照拥挤度排序,取列数即为最后的总拥挤度。获取参数按行排序后的映射关系sort。
            temp3=data_Para_sorted(index+1:index+l_f,:);%%取与temp1相同的对应参数组的行数。
            temp4=temp3(sort,:);%%按照sort映射关系进行映射得到与参数相同排序的数据。
            temp5=data_Para_cell_sorted(index+1:index+l_f,:);%%取与temp1相同的对应参数组的行数。
            temp6=temp5(sort,:);%%按照sort映射关系进行映射得到与参数相同排序的数据。
            new_pop(index+1:pop_size,:)= temp2(l_f-(pop_size-index)+1:l_f,:); %% 提取拥挤度高的参数，即提取相邻元素距离近的元素
            data_Para_new(index+1:pop_size,:)= temp4(l_f-(pop_size-index)+1:l_f,:); %% 提取拥挤度高的参数对应数据，即提取相邻元素距离近的元素
            data_Para_cell_new(index+1:pop_size,:)= temp6(l_f-(pop_size-index)+1:l_f,:); %% 提取拥挤度高的参数对应数据，即提取相邻元素距离近的元素
            index=index+l_f;
    end
    ii=ii+1;
end
