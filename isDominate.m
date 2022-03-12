function  isDomed=isDominate(sol1, sol2, flagMinMax)
%     a = all(fcon1<0);
%     b = all(fcon2<0);
%     if a && b
        if flagMinMax % max sol1>=sol2
            isDomed = (all(sol1>=sol2) && any(sol1>sol2));
        else
            isDomed = (all(sol1<=sol2) && any(sol1<sol2));
        end
%     elseif a && ~b
%         isDomed = 1;
%     elseif ~a && b
%         isDomed = 0;
%     else
%         x = sum(sum(fcon1(fcon1>0)));
%         y = sum(sum(fcon2(fcon2>0)));
%         isDomed = (x<y);
%     end

%     if flagMinMax % max
%         isDomed = (sum(sol1>=sol2, 2)==size(sol2, 2) && sum(sol1>sol2, 2)>0);
%     else
%         isDomed = (sum(sol1<=sol2, 2)==size(sol2, 2) && sum(sol1<sol2, 2)>0);
%     end
    