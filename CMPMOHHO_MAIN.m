% test function
clear, clc
close all;
score=[];
warning('off');


for i=1:1
    disp(' ')
    disp('i=')
    disp(i)
   IGD_iter=[];
t_ind=[];

%% ZDT1、2、3
%   [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('ZDT3', 'maxIterations',300, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',0, 'dimVar',30, 'numObj',2,  'numgroup', 2, 'lb',  zeros(1,30), 'ub',ones(1,30)   );

%% ZDT4
%   [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('ZDT4', 'maxIterations',300, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',0, 'dimVar',10, 'numObj',2,  'numgroup', 2, 'lb', [0,zeros(1,10-1)-5], 'ub',[1,zeros(1,10-1)+5]  );

%% ZDT6
%   [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('ZDT6', 'maxIterations',300, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',0, 'dimVar',10, 'numObj',2,  'numgroup', 2, 'lb',  zeros(1,10), 'ub',ones(1,10)   );

%% DTLZ1
%   [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('DTLZ1', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',0, 'dimVar',7, 'numObj',3,  'numgroup',3, 'lb', zeros(1,7), 'ub',ones(1,7) );

%% DTLZ2、3、4、5、6
%       [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('DTLZ6', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',0, 'dimVar',12, 'numObj',3,  'numgroup',3, 'lb', zeros(1,12), 'ub',ones(1,12) );
%% DTLZ7
%   [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('DTLZ7', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',0, 'dimVar',22, 'numObj',3,  'numgroup',3, 'lb', zeros(1,22), 'ub',ones(1,22) );
%% DTLZ8
%   [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('DTLZ8', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',1, 'dimVar',30, 'numObj',3,  'numgroup',3, 'lb', zeros(1,30), 'ub',ones(1,30) );
%% DTLZ9
%   [Archive_fitness,Archive, PF] = ...
%         CMOHHO('DTLZ1', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',1, 'dimVar',20, 'numObj',2,  'numgroup',2, 'lb', zeros(1,20), 'ub',ones(1,20) );

%% UF1、2、5、6、7
%     [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('UF2', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',0, 'dimVar',30, 'numObj', 2,  'numgroup', 2, 'lb', [0,zeros(1,30-1)-1], 'ub',ones(1,30) );

%% UF3
%     [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('UF3', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',0, 'dimVar',30, 'numObj', 2,  'numgroup', 2, 'lb', zeros(1,30), 'ub',ones(1,30) );

%% UF4
%     [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('UF4', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',0, 'dimVar',30, 'numObj', 2,  'numgroup', 2, 'lb', [0,zeros(1,30-1)-2], 'ub',[1,zeros(1,30-1)+2] );

%% UF8、9、10
%     [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('UF10', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',0, 'dimVar',30, 'numObj', 3,  'numgroup', 3, 'lb', [0,0,zeros(1,30-2)-2], 'ub',[1,1,zeros(1,30-2)+2] );

%% WFG1、2、3、4、5、6、7、8、9
%   [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('WFG2', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',0, 'dimVar',12, 'numObj',3,  'numgroup', 3, 'lb', zeros(1,12), 'ub',(2 : 2 : 2*12));

%% VNT1
%  [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('VNT1', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',1, 'dimVar',2, 'numObj',3,  'numgroup',3,'lb',[-2,-2], 'ub', [2,2] );
%% VNT2、4
%  [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('VNT1', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',1, 'dimVar',2, 'numObj',3,  'numgroup',3,'lb',[-4,-4], 'ub', [4,4] );
%% VNT3
%  [Archive_fitness,Archive, PF] = ...
%         CMPMOHHO('VNT1', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                 'plotFlag',1, 'dimVar',2, 'numObj',3,  'numgroup',3,'lb',[-3,-3], 'ub', [3,3] );
%% MaF1-6
% [Archive_fitness,Archive, PF] = ...
%          CMPMOHHO('MaF2', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
%                  'plotFlag',1, 'dimVar',12, 'numObj',3,  'numgroup',3, 'lb', zeros(1,12), 'ub',ones(1,12) );
%% LSMOP1-3
[Archive_fitness,Archive, PF] = ...
    CMPMOHHO('LSMOP3', 'maxIterations',3000, 'popSize', 100,'minmax', 'min', ...
                 'plotFlag',0, 'dimVar',300, 'numObj',3,  'numgroup',3, 'lb', zeros(1,300), 'ub',[ones(1,3-1),10.*ones(1,300-3+1)] );




    score(i) = IGD(Archive_fitness,PF);% [0,zeros(1,30-1)-2][1,zeros(1,30-1)+2]
    score=roundn(score,-5);

end 
    AVG_score=mean(score(:));
    STD2_score=(std(score,0,2)).^2;
    format shortE 
    disp(AVG_score)
    disp(STD2_score)
% save('ZDT1_IGD.mat','score');
%  save('ZDT2_IGD.mat','score');
% save('ZDT3_IGD.mat','score');
% save('ZDT4_IGD.mat','score');
% save('ZDT6_IGD.mat','score');
% save('DTLZ1_IGD.mat','score');
% save('DTLZ2_IGD.mat','score');
% save('DTLZ3_IGD.mat','score');
% save('DTLZ4_IGD.mat','score');
% save('DTLZ5_IGD.mat','score');
% save('DTLZ6_IGD.mat','score');
% save('DTLZ7_IGD.mat','score');
% save('UF1_IGD.mat','score');
% save('UF2_IGD.mat','score');
% save('UF3_IGD.mat','score');
% save('UF4_IGD.mat','score');
% save('UF5_IGD.mat','score');
% save('UF6_IGD.mat','score');
% save('UF7_IGD.mat','score');
% save('UF8_IGD.mat','score');
% save('UF9_IGD.mat','score');
% save('UF10_IGD.mat','score');
% save('WFG1_IGD.mat','score');
% save('WFG2_IGD.mat','score');
% save('WFG3_IGD.mat','score');
% save('WFG4_IGD.mat','score');
% save('WFG5_IGD.mat','score');
% save('WFG6_IGD.mat','score');
% save('WFG7_IGD.mat','score');
% save('WFG8_IGD.mat','score');
% save('WFG9_IGD.mat','score');
% save('MaF2_IGD.mat','score');
% save('MaF2_IGD.mat','score');
% save('MaF3_IGD.mat','score');
% save('MaF4_IGD.mat','score');
% save('MaF5_IGD.mat','score');
% save('MaF6_IGD.mat','score');
%  save('LSMOP1_IGD.mat','score');
%save('LSMOP2_IGD.mat','score');
save('LSMOP3_IGD.mat','score');


