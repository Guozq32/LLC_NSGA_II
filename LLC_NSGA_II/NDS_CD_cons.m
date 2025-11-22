%% Description
% 1. This function is to perform Deb's fast elitist non-domination sorting and crowding distance assignment. 
% 2. Input is in the variable 'population' with size: [size(popuation), V+M+1]
% 3. This function returns 'chromosome_NDS_CD' with size [size(population),V+M+3]
% 4. A flag 'problem_type' is used to identify whether the population is fully feasible (problem_type=0) or fully infeasible (problem_type=1) 
%    or partly feasible (problem_type=0.5). 

%% Reference:
%Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal, and T. Meyarivan, " A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II", 
%IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL. 6, No. 2, APRIL 2002. 


%% function begins
function [chromosome_NDS_CD,data_Para_sorted,data_Para_cell_sorted,front]=NDS_CD_cons(population,data_Para_inter,data_Para_cell_inter) %% front是一个列向量，front(rank).fr对应样本分层前的索引号； front(rank).fr是行向量；
global V M problem_type 

%% Initialising structures and variables
chromosome_NDS_CD1=[];
infpop=[];
front.fr=[];
struct.sp=[];
data_Parac_infeasible_sorted=[];
data_Parac_cell_infeasible_sorted=[];
rank=1;


%% Segregating feasible and infeasible solutions

if all(population(:,V+M+1)==0)  %% all类似于与运算  %% V+M+1列对应约束条件的误差
    problem_type=0;   %% 判断样本是否满足约束条件，所有误差均等于0，所以所有样均满足条件
    chromosome=population(:,1:V+M);                         % All Feasible chromosomes;  %%如果是无限制条件优化，删去最后一列  
    pop_size1=size(chromosome,1);  %% 因此样本大小只考虑优化条件即可 %% 计算矩阵第1个维度长度，如果是2则计算第二个维度长度
    data_Parac_feasible=data_Para_inter; %%所有参数均满足可行解条件，因此直接获取默认排序对应数据
    data_Parac_cell_feasible=data_Para_cell_inter;  %%所有参数均满足可行解条件，因此直接获取默认排序对应数据
elseif all(population(:,V+M+1)~=0)
    problem_type=1;   %% 有条件优化且所有的个体样本均满足优化条件，即满足不等式的约束条件
    pop_size1=0;      %% 因此样本大小为0
      infchromosome=population(:,1:V+M);         % infeasible chromosomes;  %% 不满足条件的样本
    data_Parac_infeasible=data_Para_inter; %%按照infeas_index顺序抽取不可行解对应的数据。
    data_Parac_cell_infeasible=data_Para_cell_inter;  %%按照infeas_index顺序抽取不可行解对应的数据。  
else
    problem_type=0.5;   %% 有条件优化且部分个体样本满足不等式约束条件部分不满足
    feas_index=find(population(:,V+M+1)==0);  %%部分参数满足条件，feas_index为抽取可行解在总体解中的顺序。
    chromosome=population(feas_index,1:V+M);                % Feasible chromosomes;   %% 确定满足条件的样本
    data_Parac_feasible=data_Para_inter(feas_index,:);   %%按照feas_index顺序抽取可行解对应的数据。
    data_Parac_cell_feasible=data_Para_cell_inter(feas_index,:);   %%按照feas_index顺序抽取可行解对应的数据。
    pop_size1=size(chromosome,1);
    infeas_index=find(population(:,V+M+1)~=0);%%部分参数不满足条件，infeas_index为抽取不可行解在总体解中的顺序。
    infchromosome=population(infeas_index,1:V+M+1);         % infeasible chromosomes;  %% 不满足条件的样本
    data_Parac_infeasible=data_Para_inter(infeas_index,:); %%按照infeas_index顺序抽取不可行解对应的数据。
    data_Parac_cell_infeasible=data_Para_cell_inter(infeas_index,:);  %%按照infeas_index顺序抽取不可行解对应的数据。
end

%% Handling feasible solutions 
if problem_type==0 | problem_type==0.5
    pop_size1 = size(chromosome,1);
    f1 = chromosome(:,V+1);   % objective function values
    f2 = chromosome(:,V+2); 
%Non- Domination Sorting 
% First front
for p=1:pop_size1
    struct(p).sp=find(((f1(p)-f1)<0 &(f2(p)-f2)<0) | ((f2(p)-f2)==0 &(f1(p)-f1)<0) | ((f1(p)-f1)==0 &(f2(p)-f2)<0));  %% 样本p支配的其他样本的索引号
    n(p)=length(find(((f1(p)-f1)>0 &(f2(p)-f2)>0) | ((f2(p)-f2)==0 &(f1(p)-f1)>0) | ((f1(p)-f1)==0 &(f2(p)-f2)>0)));  %% 支配P的样本的个数组成的向量
end

front(1).fr=find(n==0); %% 寻找等于0的索引号，样本的列索引号，生成第一层分类子集 
% Creating subsequent fronts
while (~isempty(front(rank).fr))
    front_indiv=front(rank).fr;      %% 层级rank样本的索引号组成的向量
    n(front_indiv)=inf;     %% 把支配此分过层级的样本数置位无穷
    chromosome(front_indiv,V+M+1)=rank; %% 第V+M+1列记录样本层数
    rank=rank+1;
    front(rank).fr=[];  %%front(rank).fr对应样本分层前的索引号； front(rank).fr是行向量；每个rank对应一个front(rank).fr，front是一个列向量
   for i = 1:length(front_indiv)  %% 层级rank的个数
        temp=struct(front_indiv(i)).sp; %% front_indiv(i)支配支配其他样本的索引号
        n(temp)=n(temp)-1;  %% 已经把减去当前分过层级的样本书，
   end 
        q=find(n==0); %% 寻找n=0对应的索引号
        front(rank).fr=[front(rank).fr q]; %%生成对应rank+1层级的分类子集   %%这里计算的rank比实际的有效样本的分层+1 
   
end
[chromosome_sorted,sort_feasible]=sortrows(chromosome,V+M+1);    % 对所有可行解按照最后一列层级数排序，sort_feasible为排序的对应映射关系。
data_Parac_feasible_sorted=data_Parac_feasible(sort_feasible,:); % 对可行解对应的数据按照sort_feasible进行重新排序，保证与初始参数的一致。
data_Parac_cell_feasible_sorted=data_Parac_cell_feasible(sort_feasible,:);% 对可行解对应的数据按照sort_feasible进行重新排序，保证与初始参数的一致。
%Crowding distance Assignment
rowsindex=1;
for i = 1:length(front)-1 %%确定多少个层级
 l_f=length(front(i).fr); %%第i层的样本数量

 if l_f > 2
     
  sorted_indf1=[];
  sorted_indf2=[];
  sortedf1=[];
  sortedf2=[];
  % sorting based on f1 and f2;
[sortedf1 sorted_indf1]=sortrows(chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+1));	%% 层级i对应的行rowsindex:(rowsindex+l_f-1); V+1列对应f1的结果；sortedf1对应返回排序的向量，sorted_indf1对应排序前的索引号
[sortedf2 sorted_indf2]=sortrows(chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+2));   %% 层级i对应的行rowsindex:(rowsindex+l_f-1); V+2列对应f2的结果；sortedf2对应返回排序的向量，sorted_indf2对应排序前的索引号

f1min=chromosome_sorted(sorted_indf1(1)+rowsindex-1,V+1); %% 因为索引号只针对chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+1)这部分进行标号，所以需要叫上rowindex生成绝对的索引标号；计算出层级i中f1的最大值和最小值
f1max=chromosome_sorted(sorted_indf1(end)+rowsindex-1,V+1);

chromosome_sorted(sorted_indf1(1)+rowsindex-1,V+M+2)=inf;   %% 第V+M+2列记录f1的拥挤度
chromosome_sorted(sorted_indf1(end)+rowsindex-1,V+M+2)=inf; %% 第V+M+3列记录f2的拥挤度

f2min=chromosome_sorted(sorted_indf2(1)+rowsindex-1,V+2); %% 同理计算层级i中f2的最大值和最小值
f2max=chromosome_sorted(sorted_indf2(end)+rowsindex-1,V+2);

chromosome_sorted(sorted_indf2(1)+rowsindex-1,V+M+3)=inf;
chromosome_sorted(sorted_indf2(end)+rowsindex-1,V+M+3)=inf;

 for j = 2:length(front(i).fr)-1
     if  (f1max - f1min == 0) | (f2max - f2min == 0)
         chromosome_sorted(sorted_indf1(j)+rowsindex-1,V+M+2)=inf;
         chromosome_sorted(sorted_indf2(j)+rowsindex-1,V+M+3)=inf;
     else
         chromosome_sorted(sorted_indf1(j)+rowsindex-1,V+M+2)=(chromosome_sorted(sorted_indf1(j+1)+rowsindex-1,V+1)-chromosome_sorted(sorted_indf1(j-1)+rowsindex-1,V+1))/(f1max-f1min); %第j个样本f1的拥挤度定义为j+1和j-1之间f1值的差值并对其标幺化拥挤度
         chromosome_sorted(sorted_indf2(j)+rowsindex-1,V+M+3)=(chromosome_sorted(sorted_indf2(j+1)+rowsindex-1,V+2)-chromosome_sorted(sorted_indf2(j-1)+rowsindex-1,V+2))/(f2max-f2min);
     end
 end

  else
    chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+M+2:V+M+3)=inf; %% V+M+2 和 V+M+3 两个优化目标的拥挤度对应拥挤度
  end
 rowsindex = rowsindex + l_f;
end
chromosome_sorted(:,V+M+4) = sum(chromosome_sorted(:,V+M+2:V+M+3),2); %V+M+4列记录V+M+2列和V+M+3列的和最为样本的拥挤度
chromosome_NDS_CD1 = [chromosome_sorted(:,1:V+M) zeros(pop_size1,1) chromosome_sorted(:,V+M+1) chromosome_sorted(:,V+M+4)]; % Final Output Variable % V+1列到V+M列为目标函数值; V+M+1列:zeros(pop_size1,1)表示误差为0  % V+M+1列对应层级；V+M+4列对应拥挤度
% front(1).fr=feas_index(find(n==0)); %% 寻找等于0的索引号，样本的列索引号，生成第一层分类子集 
end

%% Handling infeasible solutions
if problem_type==1 | problem_type==0.5
[infpop,sort_infeasible]=sortrows(infchromosome,V+M+1); %%对所有不可行解按照误差值排序，sort_infeasible为排序的对应映射关系。
data_Parac_infeasible_sorted=data_Parac_infeasible(sort_infeasible,:); % 对不可行解对应的数据按照sort_infeasible进行重新排序，保证与初始参数的一致。
data_Parac_cell_infeasible_sorted=data_Parac_cell_infeasible(sort_infeasible,:); % 对不可行解对应的数据按照sort_infeasible进行重新排序，保证与初始参数的一致。
infpop=[infpop(:,1:V+M+1) (rank:rank-1+size(infpop,1))' inf*(ones(size(infpop,1),1))]; %% size(infpop,1)计算出矩阵infpop第一位的大小，即行数； %%由于不满足约束条件，(rank:rank-1+size(infpop,1))'每个元素定义一个层级；%表示inf*(ones(size(infpop,1),1))拥挤度为无穷
for kk = (size(front,2)):(size(front,2))+(length(infchromosome))-1;
 front(kk).fr= pop_size1+1; %%继续吧不可行的样本继续分层，并标注分层后的索引号
end
end
chromosome_NDS_CD = [chromosome_NDS_CD1;infpop];  %%把不可行的样本列在可行排序的样本之后生成所有的排序样本
data_Para_sorted=[data_Parac_feasible_sorted;data_Parac_infeasible_sorted]; %%把不可行的数据列在可行排序的数据之后生成所有样本对应的数据。
data_Para_cell_sorted=[data_Parac_cell_feasible_sorted;data_Parac_cell_infeasible_sorted]; %%把不可行的数据列在可行排序的数据之后生成所有样本对应的数据。

