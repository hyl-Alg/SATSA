function Z = pen_fun(h,g)
Z=0;
lam=10^15;
% % Apply inequality constraintsӦ�ò���ʽԼ��
for k=1:length(g)
Z=Z+ lam*g(k)^2*getH(g(k));
end
% Apply equality constraintsӦ�õ�ʽԼ��
for k=1:length(h)
Z=Z+lam*h(k)^2*getHeq(h(k));
end
% Test if inequalities hold���鲻��ʽ�Ƿ���� Index function H(g) for
% inequalities����ʽ��ָ������ h (g)
function H=getH(g)
if g<=0
H=0;
else
H=1;
end
end
% Index function for equalities��ʽ��ָ������
function H=getHeq(geq)
if geq==0
H=0;
else
H=1;
end
% ----------------- end ------------------------------
end

end
