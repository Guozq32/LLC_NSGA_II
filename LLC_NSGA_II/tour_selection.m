function [parent_selected] = tour_selection(pool)
%% Description

% 1. Parents are selected from the population pool for reproduction by using binary tournament selection
%    based on the rank and crowding distance. 
% 2. An individual is selected if the rank is lesser than the other or if
%    crowding distance is greater than the other.
% 3. Input and output are of same size [pop_size, V+M+3].


%% Binary Tournament Selection
[pop_size, distance]=size(pool); %% 计算pool的行数和列数
rank=distance-1; %% 倒数第二列对应层级，得出对应层级的列标号
candidate=[randperm(pop_size);randperm(pop_size)]';  %% canditate 为 popsize×2 的随机矩阵，矩阵中的数为1到pop_size之间的证书随机置换，根据随机数挑选父样本

for i = 1: pop_size
    parent=candidate(i,:);                                  % Two parents indexes are randomly selected
 if pool(parent(1),rank)~=pool(parent(2),rank)              % For parents with different rank %%如果层级不同
    if pool(parent(1),rank)<pool(parent(2),rank)            % Checking the rank of two individuals %%随机筛选出层级高，即层技术低的作为父样本
        mincandidate=pool(parent(1),:);   %% mincandidate是行向量，对应pool中某个样本
    elseif pool(parent(1),rank)>pool(parent(2),rank)   %%随机筛选出层级高，即层技术低的最为父样本
        mincandidate=pool(parent(2),:);
    end
parent_selected(i,:)=mincandidate;                          % Minimum rank individual is selected finally
 else                                                       % for parents with same ranks   %%如果层级相同
    if pool(parent(1),distance)>pool(parent(2),distance)    % Checking the distance of two parents %%随机筛选出层级相同但拥挤度距离小的作为父样本
        maxcandidate=pool(parent(1),:);
    elseif pool(parent(1),distance)< pool(parent(2),distance)  %%随机筛选出层级高，即层技术低的最为父样本
        maxcandidate=pool(parent(2),:);
    else
        temp=randperm(2);   %% 如果两个样本层级和拥挤度都相同，则随机选取其中一个
        maxcandidate=pool(parent(temp(1)),:);
    end 
parent_selected(i,:)=maxcandidate;                          % Maximum distance individual is selected finally
end
end


    
    
    
    
    
    
    
    
    