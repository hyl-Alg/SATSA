clc;
clear;
disp('begin');
global initial_flag
initial_flag = 0;

SearchAgents_no=100; % Number of search agents
Max_iteration=1000; %Maximum number of iterations
Function_name=3;
Par = Cal_par(Function_name);
D   = Par.n;
lb   = Par.xmin;
ub   = Par.xmax;
[Best_score,Best_pos]=SATSA(SearchAgents_no,Max_iteration,lb,ub,D,Function_name);
display(['The optimal value for the optimization operation design problem of the alkylation device is ',num2str(Best_score)]);




