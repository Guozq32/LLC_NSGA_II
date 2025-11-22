function err_norm  = normalisation(error_pop)
  
%% Description
% 1. This function normalises the constraint violation of various individuals, since the range of 
%    constraint violation of every chromosome is not uniform.
% 2. Input is in the matrix error_pop with size [pop_size, number of constraints].
% 3. Output is a normalised vector, err_norm of size [pop_size,1]
 
%% Error Nomalisation
[N,nc]=size(error_pop); %% 结果对应矩阵的行和列  返回误差行数为样本个数
 con_max=0.001+max(error_pop); %% con_max 生成行向量为每一列最大值组成的行向量
 con_maxx=repmat(con_max,N,1);
 cc=error_pop./con_maxx; %% error_pop 对所有元素进行标幺化，每一列元素除以每一列的最大值
 err_norm=sum(cc,2);      % finally sum up all violations  %%把cc中每一行的所有元素相加
 

 