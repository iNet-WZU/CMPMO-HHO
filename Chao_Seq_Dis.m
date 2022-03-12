function  [X,ch1] = Chao_Seq_Dis(ch1,X,dim,ub,lb)

    if size(ub,2)==1
        ub=ones(1,dim)*ub;
        lb=ones(1,dim)*lb;
    end
    N=size(X,1);
   
    for i=1:N       
        j=randperm(dim,1);
        X(i,j)=lb(j)+ch1(1,j)*(ub(j)-lb(j));
        ch1(1,j)=(4.*ch1(1,j)).*(1-ch1(1,j));
   end
   
    for i=1:size(X,1)  
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;    
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
    end 

end
