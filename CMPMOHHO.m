
function [Archive_fitness,Archive,PF] =CMPMOHHO(Function_name, varargin)

% disp('CMPMO-HHO is now tackling your problem')

if ~mod(nargin, 2)
    error('MATLAB:narginchk:notEnoughInputs', ...
        'I have no idea about this, you can guess it');
end


%% Parameter processing
for ind = 1:2:nargin-1
    switch lower(varargin{ind})
        case 'lb'
            lb = varargin{ind + 1};
        case 'ub'
            ub = varargin{ind + 1};
        case 'numobj'
            numObj = varargin{ind + 1};
        case 'maxiterations'
            maxIterations = varargin{ind + 1};
        case 'popsize'
            popSize = varargin{ind + 1};
        case 'dimvar'
            dimVar = varargin{ind + 1};
        case 'numgroup'
            numGroup = varargin{ind + 1};
        case 'minmax'
            if strcmp(varargin{ind + 1}, 'min')
                flagMinMax = 0;
            else
                flagMinMax = 1;
            end
        case 'plotflag'
            plotFlag = varargin{ind + 1};
        otherwise
            error('The function don''t support this parameter');
    end
end
if ~exist('flagMinMax', 'var')
    error('You must specify the parameter ''MinMax''')
end

%% Initialization
% disp('CMPMOHHO is now tackling your problem')
% tic

%初始化参数
t=0; % Loop counter
T=maxIterations;
N=popSize;
% D=dimVar;
dim=dimVar;
M=numGroup;
% iter=1;
% IGD_iter=[];
% t_ind=[];
Archive = [];



% if (plotFlag)
%     figure(1); hold on;
%     h = animatedline;
%     h.LineStyle = 'none'; h.Marker = '.'; h.Color = 'r';
% end




%Initialize the locations of Harris' hawks,M populations 
    Harris_Hawks_all=initialization(M*N,dimVar,ub,lb);
    % initialize the location and Energy of the rabbit
    Rabbit_Location_all=zeros(M,dim);
    if flagMinMax
        Rabbit_Energy_all = -inf*ones(1,M);
        %     Archive_fitness = -inf*ones(M*N, numObj);
    else
        Rabbit_Energy_all = inf*ones(1,M);
        %     Archive_fitness = inf*ones(M*N, numObj);
    end
    

% for m=1:M



    %Perturbation initialization
    ch_all=zeros(1,dim);
    for d=1:dim
        x = rand();
        while 1
            if(x~=0.25&&x~=0.5&&x~=0.75&&x~=1)  
                ch_all(1,d) = x;
                break;
            else
                x=rand;
            end
        end
    end
    
%         ch_all(m,:)=ch;
% end



% pre_Harris_Hawks_all=Harris_Hawks_all;

%Generate reference point and initialize ideal points
[Z,N] = UniformPoint(N,numObj);
[obj_init,PF] = getCMPMOHHOFcn(Function_name, Harris_Hawks_all, numObj);
Zmin  = min(obj_init,[],1);
% [Archive,Archive_fitness,PF] = Update_Archive(Function_name,pre_Harris_Hawks_all,Harris_Hawks_all,ch_x,N,numObj,Z,Zmin);
%     [ch_x1,ch_all] = Chao_Seq_Dis(ch_all,Harris_Hawks_all,dim,ub,lb);
% [Archive,Archive_fitness,~] = Update_Archive(Function_name,Archive,Harris_Hawks_all,ch_x1,N,numObj,Z,Zmin);

while t<T
    
    
%     X=[];
    
    for m = 1:M

        
        X = Harris_Hawks_all((m-1)*N+1:m*N,:);
        
        for i=1:size(X,1)
        % Check boundries
        FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        [fitness,~]=getCMPMOHHOFcn(Function_name, X(i,:), numObj);
        % Update the location of Rabbit
        if fitness(:,m)<Rabbit_Energy_all(:,m)
            Rabbit_Energy_all(:,m)=fitness(:,m);
            Rabbit_Location_all(m,:)=X(i,:);
        end
        end
         Rabbit_Location =  Rabbit_Location_all(m,:);


         
        if ~isempty(Archive)
            Loction_mix = [X;Archive];
            X = Environmental_Selection(Function_name,Loction_mix,N,numObj,Z,Zmin);
            
        end
        

        
        
        %
        
        
        %      t
        %     length(find(rankedFitness(1:N,1)==1))
        E1=2*(1-(t/T)); % factor to show the decreaing energy of rabbit
        % Update the location of Harris' hawks
        for i=1:size(X,1)
            E0=2*rand()-1; %-1<E0<1
            Escaping_Energy=E1*(E0);  % escaping energy of rabbit
            
            
            %             %         Rabbit_Location = Rabbit_Location_all(rank1_rand_index,:);
            if abs(Escaping_Energy)>=1
                
                %% Exploration:
                % Harris' hawks perch randomly based on 2 strategy:
                q=rand();
                rand_Hawk_index = floor(N*rand()+1);
                X_rand = X(rand_Hawk_index, :);
                if q<0.5
                    % perch based on other family members
                    X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:));
                elseif q>=0.5
                    %                 perch on a random tall tree (random site inside group's home range)
                    X(i,:)=(Rabbit_Location(1,:)-mean(X))-rand()*((ub-lb)*rand+lb);
                end
                
            elseif abs(Escaping_Energy)<1
                %% Exploitation:
                % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
                
                %% phase 1: surprise pounce (seven kills)
                % surprise pounce (seven kills): multiple, short rapid dives by different hawks
                
                r=rand(); % probablity of each event

                
                if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege
                    X(i,:)=(Rabbit_Location)-Escaping_Energy*abs(Rabbit_Location-X(i,:));
                end
                
                if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege
                    Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                    X(i,:)=(Rabbit_Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                end
                
                %% phase 2: performing team rapid dives (leapfrog movements)
                if r<0.5 && abs(Escaping_Energy)>=0.5 % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                    
                    Jump_strength=2*(1-rand());
                    X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                    
                    [fobj1,~] = getCMPMOHHOFcn(Function_name, X1, numObj);
                    fobj2= getCMPMOHHOFcn(Function_name, X(i,:), numObj);
                    if fobj1(:,m)<fobj2(:,m)
                        X(i,:)=X1;
                    else
                        X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,dim).*Levy(dim);
                        [fobj3,~] = getCMPMOHHOFcn(Function_name, X2, numObj);
                        if fobj3(:,m)<fobj2(:,m)
                            X(i,:)=X2;
                        end
                    end

                    
                end
                
                if r<0.5 && abs(Escaping_Energy)<0.5  % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                    % hawks try to decrease their average location with the rabbit
                    Jump_strength=2*(1-rand());
                    X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X));
                    
                    [fobj1,~] = getCMPMOHHOFcn(Function_name, X1, numObj);
                    fobj2= getCMPMOHHOFcn(Function_name, X(i,:), numObj);
                    if fobj1(:,m)<fobj2(:,m)
                        X(i,:)=X1;
                    else
                        X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,dim).*Levy(dim);
                        [fobj3,~] = getCMPMOHHOFcn(Function_name, X2, numObj);
                        if fobj3(:,m)<fobj2(:,m)
                            X(i,:)=X2;
                        end
                    end
                    
                end
                %%
            end
            
        end
        
        % Check boundries
        for i=1:size(X,1)
            FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        end
        Harris_Hawks_all((m-1)*size(X,1)+1:m*size(X,1),:) =X;

    end

    %
    %     [ch_x,ch_all] = Chao_Seq_Dis(ch_all,pre_Rabbit_Location_all,dim,ub,lb);
    [ch_x,ch_all] = Chao_Seq_Dis(ch_all,Harris_Hawks_all,dim,ub,lb);
%     [Archive,Archive_fitness,PF] = Update_Archive(Function_name,pre_Harris_Hawks_all,Harris_Hawks_all,ch_x,N,numObj,Z,Zmin);
[Archive,Archive_fitness,PF] = Update_Archive(Function_name,Archive,Harris_Hawks_all,ch_x,N,numObj,Z,Zmin);
    %     Rabbit_Energy_all = Archive_fitness;
    
%     pre_Harris_Hawks_all=Harris_Hawks_all;
    

%     if (plotFlag && mod(t,round(T/60)) ==0)%
%         %         t
%         clearpoints(h);
%         if numObj==3
%             addpoints(h, Archive_fitness(:, 1), Archive_fitness(:, 2),Archive_fitness(:, 3));
%             %               plot3(Archive_fitness(:, 1),Archive_fitness(:, 2),Archive_fitness(:, 3),'.','markersize',10);
%         elseif numObj==2
%             addpoints(h, Archive_fitness(:, 1), Archive_fitness(:, 2));
%         end
%         drawnow
%     end
    
    t=t+1;
    
end

% toc
end%function


% ___________________________________
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end


