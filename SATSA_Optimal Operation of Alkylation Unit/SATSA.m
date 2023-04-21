function [Best_score,Best_pos] = SATSA(nPop,Max_iteration,lb,ub,dim,Function_name)
%Initialization parameters
N=nPop;
D=dim;
total_iter=Max_iteration;
ST=0.2;
dmax=ub;
dmin=lb;
low=ceil(N*0.1);
high=ceil(N*0.25);
g=zeros(1,30);
maxRet=5;
% Parameter setting of Nelder-Mead 
alfa=1;
beta=0.5;
gamma=2;
for run=1:1
    tic;
    trees=zeros(N,D);
    obj_trees=zeros(1,N);
    Ret=zeros(1,N);
    for i=1:N
        for j=1:D
            trees(i,j)=dmin(j)+(dmax(j)-dmin(j))*rand;
        end
        obj_trees(i)=func(trees(i,:),Function_name);
    end
  
    [minimum]=min(obj_trees);
    mins=zeros(1,total_iter);
    iter=1;
    while(iter<Max_iteration)
        iter=iter+1;
        for i=1:N
            ns=fix(low+(high-low)*rand)+1;
            
            if(ns>high)
                ns=high;
            end
            
            [minimum,min_indis]=min(obj_trees);
            bestParams=trees(min_indis,:);
            % Determine the number of ns 
            if Ret(i)>=maxRet
                ns=fix((low+(high-low)*rand)*1.5)+1;
            end
            seeds=zeros(ns,D);
            obj_seeds=zeros(1,ns);
            % Generate ns seeds
            for j=1:ns
                komsu=fix(rand*N)+1;
                while(i==komsu)
                    komsu=fix(rand*N)+1;
                end
                seeds(j,:)=trees(i,:);
                for d=1:D
                    if(rand<ST)
                        sigma=(rand-0.5)+((Max_iteration-iter)/Max_iteration)*(randn+0.5);
                        seeds(j,d)=trees(i,d)+(bestParams(d)-trees(komsu,d))*sigma;
                        if(seeds(j,d)>dmax(d))
                            seeds(j,d)=dmax(d);
                        end
                        if(seeds(j,d)<dmin(d))
                            seeds(j,d)=dmin(d);
                        end
                    else
                        if Ret(i)<=maxRet
                            seeds(j,d)=trees(i,d)+(trees(i,d)-trees(komsu,d))*(rand-0.5)*2;
                        else
                            seeds(j,d)=bestParams(d)+(trees(i,d)-trees(komsu,d))*rand;
                        end
                        if(seeds(j,d)>dmax(d))
                            seeds(j,d)=dmax(d);
                        end
                        if(seeds(j,d)<dmin(d))
                            seeds(j,d)=dmin(d);
                        end
                    end
                end
                obj_seeds(j)=func(seeds(j,:),Function_name);
            end
            % Update the position of ith tree
            [mintohum,mintohum_indis]=min(obj_seeds);
            if(mintohum<obj_trees(i))
                trees(i,:)=seeds(mintohum_indis,:);
                obj_trees(i)=mintohum;
                Ret(i)=0;
            else
                Ret(i)=Ret(i)+1;
            end
        end
        
        % Update the position of tree by Simplex method
        [fVal,SortIndex]=sort(obj_trees);
        Num=fix(1/2*N*((Max_iteration-iter)/Max_iteration))+1;   % Calculate the value of num by Eq.(10)
        z=0;
        for n=(N-Num+1):N
            z=z+1;
            I(z)=SortIndex(n);
            Val(z)=fVal(n);
        end
        CenT=trees(SortIndex(1:(N-Num)),:);
        Trees=trees(I,:); 
        xc=mean(CenT);
        
        for m=1:size(Trees,1)
            ybest=Val(m);
            xworst=Trees(m,:);
            xr=(1+alfa)*xc-alfa*xworst;   % Reflection operation performed by Eq.(11)
            Flagub=xr>ub;
            Flaglb=xr<lb;
            xr=(xr.*(~(Flagub+Flaglb)))+ub.*Flagub+lb.*Flaglb;
            
            freflec=func(xr,Function_name);
            if freflec<ybest
                xe=gamma*xr+(1-gamma)*xc; % Extension operation by Eq. (12)
                Flag1ub=xe>ub;
                Flag1lb=xe<lb;
                xe=(xe.*(~(Flag1ub+Flag1lb)))+ub.*Flag1ub+lb.*Flag1lb;
                fexp=func(xe,Function_name);
                if fexp<ybest
                    trees(I(m),:)=xe;
                    obj_trees(I(m))=fexp;
                else
                    trees(I(m),:)=xr;
                    obj_trees(I(m))=freflec;
                end
            else
                xcont=beta*xworst+(1-beta)*xc; % contraction operation according to Eq. (13)
                Flag2ub=xcont>ub;
                Flag2lb=xcont<lb;
                xcont=(xcont.*(~(Flag2ub+Flag2lb)))+ub.*Flag2ub+lb.*Flag2lb;
                fcont=func(xcont,Function_name);
                if fcont<ybest
                    trees(I(m),:)=xcont;
                    obj_trees(I(m))=fcont;
                else
                    trees(I(m),:)=1/2*(trees(I(m),:)+bestParams);
                    obj_trees(I(m))=func(trees(I(m),:),Function_name);
                end
            end
        end
        [min_tree,min_tree_index]=min(obj_trees);
        if(min_tree<minimum)
            minimum=min_tree;
            bestParams=trees(min_tree_index,:);
        end
        mins(iter)=minimum;
    end

    fprintf('Iter=%d .... min=%g \n',iter,minimum);
end
Best_pos = bestParams;
Best_score = minimum;
end