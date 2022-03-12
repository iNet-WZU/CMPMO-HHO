function [Rabbit_Location_all,Rabbit_Energy_all,PF] = Update_Archive(Function_name,pre_Rabbit_Location_all,Rabbit_Location_all,ch_x,N,numObj,Z,Zmin)
% function [Rabbit_Location_all,Rabbit_Energy_all,PF] = Update_Archive(Function_name,Rabbit_Location_all,ch_x,N,numObj,Z,Zmin)


%     [ch_x,ch_all] = Chao_Seq_Dis(ch_all,pre_Rabbit_Location_all,dim,ub,lb);  
    Mix_loc=[pre_Rabbit_Location_all;Rabbit_Location_all;ch_x];%;array_num
%     Mix_loc=[Rabbit_Location_all;ch_x];%;array_num
    Mix_loc_unrep=unique(Mix_loc,'rows');
    [sub_fitness,  ~] =getCMPMOHHOFcn(Function_name, Mix_loc_unrep, numObj);
    Zmin       = min([Zmin;sub_fitness],[],1);
    Rabbit_Location_all = Environmental_Selection(Function_name,Mix_loc_unrep,N,numObj,Z,Zmin);
    
    [Rabbit_Energy_all,PF] =getCMPMOHHOFcn(Function_name, Rabbit_Location_all, numObj);
end