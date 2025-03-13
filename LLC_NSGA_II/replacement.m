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
    l_f=length(front(ii).fr);  %%l_f is [ 2pop _ size, 1 ], which records the dominance level
    if index+l_f < pop_size  %% If the stratification is still less than the total population, the stratification is regarded as a new population.
        new_pop(index+1:index+l_f,:)= population_inter_sorted(index+1:index+l_f,:);
        data_Para_new(index+1:index+l_f,:)= data_Para_sorted(index+1:index+l_f,:);
        data_Para_cell_new(index+1:index+l_f,:)=data_Para_cell_sorted(index+1:index+l_f,:);
        index=index+l_f;
    else
            temp1=population_inter_sorted(index+1:index+l_f,:); 
            [temp2,sort]=sortrows(temp1,size(temp1,2)); %% According to the crowding degree, the number of columns is the final total crowding degree. Get the mapping relationship sort after the parameters are sorted by row.
            temp3=data_Para_sorted(index+1:index+l_f,:);%%Take the same number of rows of the corresponding parameter group as temp1.
            temp4=temp3(sort,:);%%According to the sort mapping relationship, the data with the same order as the parameters are obtained.
            temp5=data_Para_cell_sorted(index+1:index+l_f,:);
            temp6=temp5(sort,:);
            new_pop(index+1:pop_size,:)= temp2(l_f-(pop_size-index)+1:l_f,:); %% Extract the parameters with high crowding degree, that is, extract the elements close to the adjacent elements.
            data_Para_new(index+1:pop_size,:)= temp4(l_f-(pop_size-index)+1:l_f,:); %% Extract the corresponding data of the parameters with high crowding, that is, extract the elements close to the adjacent elements.
            data_Para_cell_new(index+1:pop_size,:)= temp6(l_f-(pop_size-index)+1:l_f,:); %%Extract the corresponding data of the parameters with high crowding, that is, extract the elements close to the adjacent elements.
            index=index+l_f;
    end
    ii=ii+1;
end
