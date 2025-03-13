function err_norm  = normalisation(error_pop)
  
%% Description
% 1. This function normalises the constraint violation of various individuals, since the range of 
%    constraint violation of every chromosome is not uniform.
% 2. Input is in the matrix error_pop with size [pop_size, number of constraints].
% 3. Output is a normalised vector, err_norm of size [pop_size,1]
 
%% Error Nomalisation
[N,nc]=size(error_pop); %% Results The number of rows and columns of the corresponding matrix return error is the number of samples.
 con_max=0.001+max(error_pop); %% con_max Generate a row vector consisting of the maximum value of each column
 con_maxx=repmat(con_max,N,1);
 cc=error_pop./con_maxx; %% error_pop All elements are normalized, and each column element is divided by the maximum value of each column
 err_norm=sum(cc,2);      
 

 