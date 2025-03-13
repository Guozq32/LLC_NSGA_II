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
function [chromosome_NDS_CD,data_Para_sorted,data_Para_cell_sorted,front]=NDS_CD_cons(population,data_Para_inter,data_Para_cell_inter) %% Front is a column vector, front ( rank ).fr corresponds to the index number of the sample before stratification ; front ( rank ).fr is a row vector ;
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

if all(population(:,V+M+1)==0)  %% all is similar to the AND operation  %% The error of V + M + 1 column corresponding to the constraint condition
    problem_type=0;   %% To determine whether the sample meets the constraints, all errors are equal to 0, so all samples meet the conditions.
    chromosome=population(:,1:V+M);                         % All Feasible chromosomes;  %%If it is an unconstrained optimization, delete the last column 
    pop_size1=size(chromosome,1);   
    data_Parac_feasible=data_Para_inter; %%All parameters satisfy the feasible solution condition, so the default sort data is obtained directly.
    data_Parac_cell_feasible=data_Para_cell_inter;  
elseif all(population(:,V+M+1)~=0)
    problem_type=1;   %% Conditional optimization and all individual samples meet the optimization conditions, that is, meet the inequality constraints.
    pop_size1=0;      
      infchromosome=population(:,1:V+M);         % infeasible chromosomes;
    data_Parac_infeasible=data_Para_inter; %%Extract the data corresponding to the infeasible solution in the order of infeas _ index.
    data_Parac_cell_infeasible=data_Para_cell_inter;  
else
    problem_type=0.5;   %% Conditional optimization and some individual samples meet the inequality constraints are partially not satisfied.
    feas_index=find(population(:,V+M+1)==0);  %%Some parameters meet the conditions, feas _ index is the order of extracting feasible solutions in the overall solution.
    chromosome=population(feas_index,1:V+M);               
    data_Parac_feasible=data_Para_inter(feas_index,:);   %%Extract the data corresponding to the feasible solution according to the feas _ index order.
    data_Parac_cell_feasible=data_Para_cell_inter(feas_index,:);  
    pop_size1=size(chromosome,1);
    infeas_index=find(population(:,V+M+1)~=0);%%Some parameters do not meet the conditions, and infeas _ index is the order of extracting infeasible solutions in the overall solution.
    infchromosome=population(infeas_index,1:V+M+1);         % infeasible chromosomes;  %% Samples that do not meet the conditions
    data_Parac_infeasible=data_Para_inter(infeas_index,:); 
    data_Parac_cell_infeasible=data_Para_cell_inter(infeas_index,:);  
end

%% Handling feasible solutions 
if problem_type==0 | problem_type==0.5
    pop_size1 = size(chromosome,1);
    f1 = chromosome(:,V+1);   % objective function values
    f2 = chromosome(:,V+2); 
%Non- Domination Sorting 
% First front
for p=1:pop_size1
    struct(p).sp=find(((f1(p)-f1)<0 &(f2(p)-f2)<0) | ((f2(p)-f2)==0 &(f1(p)-f1)<0) | ((f1(p)-f1)==0 &(f2(p)-f2)<0));  %% The index number of other samples dominated by sample p
    n(p)=length(find(((f1(p)-f1)>0 &(f2(p)-f2)>0) | ((f2(p)-f2)==0 &(f1(p)-f1)>0) | ((f1(p)-f1)==0 &(f2(p)-f2)>0)));  %% The vector composed of the number of samples dominating P
end

front(1).fr=find(n==0); %% Find the index number equal to 0, the column index number of the sample, and generate the first layer classification subset. 
% Creating subsequent fronts
while (~isempty(front(rank).fr))
    front_indiv=front(rank).fr;      %% The vector composed of the index number of the hierarchical rank sample
    n(front_indiv)=inf;     %% Set the number of samples that dominate this hierarchy to infinity.
    chromosome(front_indiv,V+M+1)=rank; %% Column V + M + 1 records the number of sample layers.
    rank=rank+1;
    front(rank).fr=[];  
   for i = 1:length(front_indiv)  %% Number of rank levels
        temp=struct(front_indiv(i)).sp; %% Front _ indiv ( i ) dominates the index number of other samples.
        n(temp)=n(temp)-1;  
   end 
        q=find(n==0); %% Find the index number corresponding to n = 0.
        front(rank).fr=[front(rank).fr q]; %%Generate a classification subset corresponding to the rank + 1 level
   
end
[chromosome_sorted,sort_feasible]=sortrows(chromosome,V+M+1);    % For all feasible solutions, sort _ feasible is the corresponding mapping relation of sorting according to the last column of hierarchical series.
data_Parac_feasible_sorted=data_Parac_feasible(sort_feasible,:); % The data corresponding to the feasible solution is reordered according to sort _ feasible to ensure consistency with the initial parameters.
data_Parac_cell_feasible_sorted=data_Parac_cell_feasible(sort_feasible,:);% The data corresponding to the feasible solution is reordered according to sort _ feasible to ensure consistency with the initial parameters.
%Crowding distance Assignment
rowsindex=1;
for i = 1:length(front)-1 
 l_f=length(front(i).fr); 

 if l_f > 2
     
  sorted_indf1=[];
  sorted_indf2=[];
  sortedf1=[];
  sortedf2=[];
  % sorting based on f1 and f2;
[sortedf1 sorted_indf1]=sortrows(chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+1));	
[sortedf2 sorted_indf2]=sortrows(chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+2));  

f1min=chromosome_sorted(sorted_indf1(1)+rowsindex-1,V+1); 
f1max=chromosome_sorted(sorted_indf1(end)+rowsindex-1,V+1);

chromosome_sorted(sorted_indf1(1)+rowsindex-1,V+M+2)=inf;   
chromosome_sorted(sorted_indf1(end)+rowsindex-1,V+M+2)=inf; 

f2min=chromosome_sorted(sorted_indf2(1)+rowsindex-1,V+2); 
f2max=chromosome_sorted(sorted_indf2(end)+rowsindex-1,V+2);

chromosome_sorted(sorted_indf2(1)+rowsindex-1,V+M+3)=inf;
chromosome_sorted(sorted_indf2(end)+rowsindex-1,V+M+3)=inf;

 for j = 2:length(front(i).fr)-1
     if  (f1max - f1min == 0) | (f2max - f2min == 0)
         chromosome_sorted(sorted_indf1(j)+rowsindex-1,V+M+2)=inf;
         chromosome_sorted(sorted_indf2(j)+rowsindex-1,V+M+3)=inf;
     else
         chromosome_sorted(sorted_indf1(j)+rowsindex-1,V+M+2)=(chromosome_sorted(sorted_indf1(j+1)+rowsindex-1,V+1)-chromosome_sorted(sorted_indf1(j-1)+rowsindex-1,V+1))/(f1max-f1min); 
         chromosome_sorted(sorted_indf2(j)+rowsindex-1,V+M+3)=(chromosome_sorted(sorted_indf2(j+1)+rowsindex-1,V+2)-chromosome_sorted(sorted_indf2(j-1)+rowsindex-1,V+2))/(f2max-f2min);
     end
 end

  else
    chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+M+2:V+M+3)=inf; 
  end
 rowsindex = rowsindex + l_f;
end
chromosome_sorted(:,V+M+4) = sum(chromosome_sorted(:,V+M+2:V+M+3),2); 
chromosome_NDS_CD1 = [chromosome_sorted(:,1:V+M) zeros(pop_size1,1) chromosome_sorted(:,V+M+1) chromosome_sorted(:,V+M+4)]; 

end

%% Handling infeasible solutions
if problem_type==1 | problem_type==0.5
[infpop,sort_infeasible]=sortrows(infchromosome,V+M+1); %For all infeasible solutions, sort _ infeasible is the corresponding mapping relation sorted by error value.
data_Parac_infeasible_sorted=data_Parac_infeasible(sort_infeasible,:); % The data corresponding to the infeasible solution is reordered according to sort _ infeasible to ensure consistency with the initial parameters.
data_Parac_cell_infeasible_sorted=data_Parac_cell_infeasible(sort_infeasible,:); % The data corresponding to the infeasible solution is reordered according to sort _ infeasible to ensure consistency with the initial parameters.
infpop=[infpop(:,1:V+M+1) (rank:rank-1+size(infpop,1))' inf*(ones(size(infpop,1),1))]; 
for kk = (size(front,2)):(size(front,2))+(length(infchromosome))-1;
 front(kk).fr= pop_size1+1; 
end
end
chromosome_NDS_CD = [chromosome_NDS_CD1;infpop];  %%The infeasible sample columns are generated after the feasible sorted samples to generate all sorted samples.
data_Para_sorted=[data_Parac_feasible_sorted;data_Parac_infeasible_sorted]; %%The infeasible data is columned after the feasible sorted data to generate the data corresponding to all samples.¡£
data_Para_cell_sorted=[data_Parac_cell_feasible_sorted;data_Parac_cell_infeasible_sorted]; %%The infeasible data column is generated after the feasible sorted data to generate the data corresponding to all samples.

