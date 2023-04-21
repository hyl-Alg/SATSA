clear all 
clc
D=30;
Xmin=-100;
Xmax=100;
pop_size=100;
iter_max=1000;
func_num=1;
fhd=str2func('cec17_func');
[gbest,gbestval,fitcount]=SATSA(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
display(['The best value for SATSA to solve F1 is ',num2str(gbestval)]);




