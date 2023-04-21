function Z = pen_fun(h,g)
Z=0;
lam=10^15;
% % Apply inequality constraints应用不等式约束
for k=1:length(g)
Z=Z+ lam*g(k)^2*getH(g(k));
end
% Apply equality constraints应用等式约束
for k=1:length(h)
Z=Z+lam*h(k)^2*getHeq(h(k));
end
% Test if inequalities hold检验不等式是否成立 Index function H(g) for
% inequalities不等式的指数函数 h (g)
function H=getH(g)
if g<=0
H=0;
else
H=1;
end
end
% Index function for equalities等式的指数函数
function H=getHeq(geq)
if geq==0
H=0;
else
H=1;
end
% ----------------- end ------------------------------
end

end
