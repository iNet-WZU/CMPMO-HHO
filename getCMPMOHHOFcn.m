function [PopObj, P] = getCMPMOHHOFcn(F, PopDec, numObj)

[N,D]  = size(PopDec);
switch F
    
    case 'ZDT1'
        PopObj(:, 1) = PopDec(:, 1);
        g = 1 + 9*mean(PopDec(:,2:end),2);
        h = 1 - (PopObj(:,1)./g).^0.5;
        PopObj(:,2) = g.*h;
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5;
        %             PopCon(:, 1:D) = -PopDec;
        %             PopCon(:, D+1: 2*D) = PopDec - 1;
    case 'ZDT2'
        PopObj(:,1) = PopDec(:,1);
        g = 1 + 9*mean(PopDec(:,2:end),2);
        h = 1 - (PopObj(:,1)./g).^2;
        PopObj(:,2) = g.*h;
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^2;
        %             PopCon(:, 1:D) = -PopDec;
        %             PopCon(:, D+1: 2*D) = PopDec - 1;
    case 'ZDT3'
        PopObj(:,1) = PopDec(:,1);
        g = 1 + 9*mean(PopDec(:,2:end),2);
        h = 1 - (PopObj(:,1)./g).^0.5 - PopObj(:,1)./g.*sin(10*pi*PopObj(:,1));
        PopObj(:,2) = g.*h;
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5 - P(:,1).*sin(10*pi*P(:,1));
        P      = P(NDSort(P,1)==1,:);
        %             PopCon(:, 1:D) = -PopDec;
        %             PopCon(:, D+1: 2*D) = PopDec - 1;
    case 'ZDT4'
        PopObj(:,1) = PopDec(:,1);
        g = 1 + 10*(size(PopDec,2)-1) + sum(PopDec(:,2:end).^2-10*cos(4*pi*PopDec(:,2:end)),2);
        h = 1 - (PopObj(:,1)./g).^0.5;
        PopObj(:,2) = g.*h;
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5;
        %             PopCon(:, 1:D) = -PopDec;
        %             PopCon(:, D+1: 2*D) = PopDec - 1;
    case 'ZDT5'%new add
        D        = ceil(max(D-30,1)/5)*5 + 30;
        u      = zeros(size(PopDec,1),1+(size(PopDec,2)-30)/5);
        u(:,1) = sum(PopDec(:,1:30),2);
        for i = 2 : size(u,2)
            u(:,i) = sum(PopDec(:,(i-2)*5+31:(i-2)*5+35),2);
        end
        v           = zeros(size(u));
        v(u<5)      = 2 + u(u<5);
        v(u==5)     = 1;
        PopObj(:,1) = 1 + u(:,1);
        g           = sum(v(:,2:end),2);
        h           = 1./PopObj(:,1);
        PopObj(:,2) = g.*h;
        P(:,1) = 1 : 31;
        %             P(:,2) = (obj.Global.D-30)./5./P(:,1);
        %             PopCon(:, 1:D) = -PopDec;
        %             PopCon(:, D+1: 2*D) = PopDec - 1;
    case 'ZDT6'
        PopObj(:,1) = 1 - exp(-4*PopDec(:,1)).*sin(6*pi*PopDec(:,1)).^6;
        g = 1 + 9*mean(PopDec(:,2:end),2).^0.25;
        h = 1 - (PopObj(:,1)./g).^2;
        PopObj(:,2) = g.*h;
        
        minf1  = 0.280775;
        P(:,1) = (minf1:(1-minf1)/(N-1):1)';
        P(:,2) = 1 - P(:,1).^2;
        %             PopCon(:, 1:D) = -PopDec;
        %             PopCon(:, D+1: 2*D) = PopDec - 1;
        %         case 'UF1'
        %             D  = size(X,2);
        %             J1 = 3 : 2 : D;
        %             J2 = 2 : 2 : D;
        %             Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        %             PopObj(:,1) = X(:,1)         + 2*mean(Y(:,J1).^2,2);
        %             PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean(Y(:,J2).^2,2);
    case 'VNT1'
        PopObj(:,1) = PopDec(:,1).^2 + (PopDec(:,2)-1).^2;
        PopObj(:,2) = PopDec(:,1).^2 + (PopDec(:,2)+1).^2 + 1;
        PopObj(:,3) = (PopDec(:,1)-1).^2 + PopDec(:,2).^2 + 2;
        %             PopCon(:, D+1: 2*D) = PopDec - 2;
        X = ReplicatePoint(N,2)*4-2;
        P = CalObjVNT1(X);
        P = P(NDSort(P,1)==1,:);%排序这里有问题
    case 'VNT2'
        PopObj(:,1) = (PopDec(:,1)-2).^2/2 + (PopDec(:,2)+1).^2/13 + 3;
        PopObj(:,2) = (PopDec(:,1)+PopDec(:,2)-3).^2/36 + (-PopDec(:,1)+PopDec(:,2)+2).^2/8 - 17;
        PopObj(:,3) = (PopDec(:,1)+2*PopDec(:,2)-1).^2/175 + (2*PopDec(:,2)-PopDec(:,1)).^2/17 - 13;
        
        X = ReplicatePoint(N,2)*8-4;%FUNCTION ReplicatePoint&&CalObj
        P = CalObjVNT2(X);
        P = P(NDSort(P,1)==1,:);
    case 'VNT3'
        temp = PopDec(:,1).^2 + PopDec(:,2).^2;
        PopObj(:,1) = 0.5*temp + sin(temp);
        PopObj(:,2) = (3*PopDec(:,1)-2*PopDec(:,2)+4).^2/8 + (PopDec(:,1)-PopDec(:,2)+1).^2/27 + 15;
        PopObj(:,3) = 1./(temp+1) - 1.1*exp(-temp);
        %             PopCon(:, 1:D) = -2-PopDec;
        %             PopCon(:, D+1: 2*D) = PopDec - 2;
        X = ReplicatePoint(N,2)*6-3;
        P = CalObjVNT3(X);
        P = P(NDSort(P,1)==1,:);
    case 'VNT4'
        PopObj(:,1) = (PopDec(:,1)-2).^2/2 + (PopDec(:,2)+1).^2/13 + 3;
        PopObj(:,2) = (PopDec(:,1)+PopDec(:,2)-3).^2/175 + (2*PopDec(:,2)-PopDec(:,1)).^2/17 - 13;
        PopObj(:,3) = (3*PopDec(:,1)-2*PopDec(:,2)+4).^2/8 + (PopDec(:,1)-PopDec(:,2)+1).^2/27 + 15;
        %             PopCon(:,1) = PopDec(:,2) + 4*PopDec(:,1) - 4;
        %             PopCon(:,2) = -1 - PopDec(:,1);
        %             PopCon(:,3) = PopDec(:,1) - 2 - PopDec(:,2);
        N = ceil(sqrt(N));
        x = linspace(-1,1.2,N);
        X = [];
        for i = 1 : N
            X = [X;repmat(x(i),N,1),linspace(x(i)-2,-4*x(i)+4,N)'];
        end
        P = CalObjVNT4(X);
        P = P(NDSort(P,1)==1,:);
    case 'DTLZ1'
        M      = numObj;
        g      = 100*(D-M+1+sum((PopDec(:,M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,M:end)-0.5)),2));
        PopObj = 0.5*repmat(1+g,1,M).*fliplr(cumprod([ones(size(PopDec,1),1),PopDec(:,1:M-1)],2)).*[ones(size(PopDec,1),1),1-PopDec(:,M-1:-1:1)];
        P = UniformPoint(N,M)/2;
        %             PopCon(:, 1:D) = -PopDec;
        %             PopCon(:, D+1: 2*D) = PopDec - 1;
    case 'DTLZ2'
        M = numObj;
        g = sum((PopDec(:,M:end)-0.5).^2,2);
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        %             PopCon(:, 1:D) = -PopDec;
        %             PopCon(:, D+1: 2*D) = PopDec - 1;
    case 'DTLZ3'
        M      = numObj;
        g      = 100*(D-M+1+sum((PopDec(:,M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,M:end)-0.5)),2));
        %             PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(PopDec,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(PopDec,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        %             PopCon(:, 1:D) = -PopDec;
        %             PopCon(:, D+1: 2*D) = PopDec - 1;
    case 'DTLZ4'
        M      = numObj;
        PopDec(:,1:M-1) = PopDec(:,1:M-1).^100;
        g      = sum((PopDec(:,M:end)-0.5).^2,2);
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1, M);
    case 'DTLZ5'
        M      = numObj;
        g      = sum((PopDec(:,M:end)-0.5).^2,2);
        Temp   = repmat(g,1,M-2);
        PopDec(:,2:M-1) = (1+2*Temp.*PopDec(:,2:M-1))./(2+2*Temp);
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        P = [0:1/(N-1):1;1:-1/(N-1):0]';
        P = P./repmat(sqrt(sum(P.^2,2)),1,size(P,2));
        P = [P(:,ones(1,M-2)),P];
        P = P./sqrt(2).^repmat([M-2,M-2:-1:0],size(P,1),1);
    case 'DTLZ6'
        M      = numObj;
        g      = sum((PopDec(:,M:end)-0.5).^2,2);
        Temp   = repmat(g,1,M-2);
        PopDec(:,2:M-1) = (1+2*Temp.*PopDec(:,2:M-1))./(2+2*Temp);
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        P = [0:1/(N-1):1;1:-1/(N-1):0]';
        P = P./repmat(sqrt(sum(P.^2,2)),1,size(P,2));
        P = [P(:,ones(1,M-2)),P];
        P = P./sqrt(2).^repmat([M-2,M-2:-1:0],size(P,1),1);
    case 'DTLZ7'
        M      = numObj;
        PopObj          = zeros(size(PopDec,1),M);
        g               = 1+9*mean(PopDec(:,M:end),2);
        PopObj(:,1:M-1) = PopDec(:,1:M-1);
        PopObj(:,M)     = (1+g).*(M-sum(PopObj(:,1:M-1)./(1+repmat(g,1,M-1)).*(1+sin(3*pi.*PopObj(:,1:M-1))),2));
        interval     = [0,0.251412,0.631627,0.859401];
        median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
        X            = ReplicatePoint(N,M-1);
        X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
        X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
        P            = [X,2*(M-sum(X/2.*(1+sin(3*pi.*X)),2))];
    case 'DTLZ8'
        M      = numObj;
        PopObj = zeros(N,M);
        for m = 1 : M
            PopObj(:,m) = mean(PopDec(:,(m-1)*D/M+1:m*D/M),2);
        end
        %             PopCon = zeros(size(PopObj,1),M);
        %             PopCon(:,1:M-1) = 1 - repmat(PopObj(:,M),1,M-1) - 4*PopObj(:,1:M-1);
        %             if M == 2
        %                 PopCon(:,M) = 0;
        %             else
        %                 minValue    = sort(PopObj(:,1:M-1),2);
        %                 PopCon(:,M) = 1 - 2*PopObj(:,M) - sum(minValue(:,1:2),2);
        %             end
        if M == 2
            temp = (0:1/(N-1):1)';
            P    = [(1-temp)/4,temp];
        else
            temp = UniformPoint(N/(M-1),3);
            temp(:,3) = temp(:,3) / 2;
            temp = temp(temp(:,1)>=(1-temp(:,3))/4 & temp(:,1)<=temp(:,2) & temp(:,3)<=1/3,:);
            P    = [repmat(temp(:,2),M-1,M-1),repmat(temp(:,3),M-1,1)];
            for i = 1 : M-1
                P((i-1)*size(temp,1)+1:i*size(temp,1),i) = temp(:,1);
            end
            gap  = sort(unique(P(:,M)));
            gap  = gap(2) - gap(1);
            temp = (1/3:gap:1)';
            P    = [P;repmat((1-temp)/4,1,M-1),temp];
            P    = unique(P,'rows');
        end
    case 'DTLZ9'
        M      = numObj;
        PopDec = PopDec.^0.1;
        PopObj = zeros(N,M);
        for m = 1 : M
            PopObj(:,m) = sum(PopDec(:,(m-1)*D/M+1:m*D/M),2);
        end
        %             PopCon = 1 - repmat(PopObj(:,M).^2,1,M-1) - PopObj(:,1:M-1).^2;
        Temp = (0:1/(N-1):1)';
        P    = [repmat(cos(0.5.*pi.*Temp),1,M-1),sin(0.5.*pi.*Temp)];
    case 'WFG1'
        M      = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);
        
        %             z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
        z01 = PopDec./repmat(2:2:D*2,N,1);
        t1 = zeros(N,K+L);
        t1(:,1:K)     = z01(:,1:K);
        t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);
        
        t2 = zeros(N,K+L);
        t2(:,1:K)     = t1(:,1:K);
        t2(:,K+1:end) = b_flat(t1(:,K+1:end),0.8,0.75,0.85);
        
        t3 = zeros(N,K+L);
        t3 = b_poly(t2,0.02);
        
        t4 = zeros(N,M);
        for i = 1 : M-1
            t4(:,i) = r_sum(t3(:,(i-1)*K/(M-1)+1:i*K/(M-1)),2*((i-1)*K/(M-1)+1):2:2*i*K/(M-1));
        end
        t4(:,M) = r_sum(t3(:,K+1:K+L),2*(K+1):2:2*(K+L));
        
        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t4(:,M),A(i)).*(t4(:,i)-0.5)+0.5;
        end
        x(:,M) = t4(:,M);
        
        h      = convex(x);
        h(:,M) = mixed(x);
        PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
        
        %PF
        
        P = UniformPoint(N,M);
        c = ones(size(P,1),M);
        for i = 1 : size(P,1)
            for j = 2 : M
                temp = P(i,j)/P(i,1)*prod(1-c(i,M-j+2:M-1));
                c(i,M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
            end
        end
        x = acos(c)*2/pi;
        temp = (1-sin(pi/2*x(:,2))).*P(:,M)./P(:,M-1);
        a = 0 : 0.0001 : 1;
        E = abs(temp*(1-cos(pi/2*a))-1+repmat(a+cos(10*pi*a+pi/2)/10/pi,size(x,1),1));
        [~,rank] = sort(E,2);
        for i = 1 : size(x,1)
            x(i,1) = a(min(rank(i,1:10)));
        end
        P      = convex(x);
        P(:,M) = mixed(x);
        P      = repmat(2:2:2*M,size(P,1),1).*P;
    case 'WFG2'
        M      = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);
        
        z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
        
        t1 = zeros(N,K+L);
        t1(:,1:K)     = z01(:,1:K);
        t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);
        
        t2 = zeros(N,K+L/2);
        t2(:,1:K) = t1(:,1:K);
        % Same as <t2(:,i)=r_nonsep(t1(:,K+2*(i-K)-1:K+2*(i-K)),2)>
        t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;
        % ---------------------------------------------------------
        
        t3 = zeros(N,M);
        for i = 1 : M-1
            t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t3(:,M) = r_sum(t2(:,K+1:K+L/2),ones(1,L/2));
        
        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
        end
        x(:,M) = t3(:,M);
        
        h      = convex(x);
        h(:,M) = disc(x);
        PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
        
        P = UniformPoint(N,M);
        c = ones(size(P,1),M);
        for i = 1 : size(P,1)
            for j = 2 : M
                temp = P(i,j)/P(i,1)*prod(1-c(i,M-j+2:M-1));
                c(i,M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
            end
        end
        x = acos(c)*2/pi;
        temp = (1-sin(pi/2*x(:,2))).*P(:,M)./P(:,M-1);
        a = 0 : 0.0001 : 1;
        E = abs(temp*(1-cos(pi/2*a))-1+repmat(a.*cos(5*pi*a).^2,size(x,1),1));
        [~,rank] = sort(E,2);
        for i = 1 : size(x,1)
            x(i,1) = a(min(rank(i,1:10)));
        end
        P      = convex(x);
        P(:,M) = disc(x);
        P      = P(NDSort(P,1)==1,:);
        P      = repmat(2:2:2*M,size(P,1),1).*P;
    case 'WFG3'
        M      = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = [1,zeros(1,M-2)];
        
        z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
        
        t1 = zeros(N,K+L);
        t1(:,1:K)     = z01(:,1:K);
        t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);
        
        t2 = zeros(N,K+L/2);
        t2(:,1:K) = t1(:,1:K);
        % Same as <t2(:,i)=r_nonsep(t1(:,K+2*(i-K)-1:K+2*(i-K)),2)>
        t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;
        % ---------------------------------------------------------
        
        t3 = zeros(N,M);
        for i = 1 : M-1
            t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t3(:,M) = r_sum(t2(:,K+1:K+L/2),ones(1,L/2));
        
        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
        end
        x(:,M) = t3(:,M);
        
        h      = linear(x);
        PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
        X = (0:1/(N-1):1)';
        X = [X,zeros(N,M-2)+0.5,zeros(N,1)];
        P = linear(X);
        P = repmat(2:2:2*M,size(P,1),1).*P;
    case 'WFG4'
        M      = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);
        
        z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
        
        t1 = zeros(N,K+L);
        t1 = s_multi(z01,30,10,0.35);
        
        t2 = zeros(N,M);
        for i = 1 : M-1
            t2(:,i) = r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t2(:,M) = r_sum(t1(:,K+1:K+L),ones(1,L));
        
        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
        end
        x(:,M) = t2(:,M);
        
        h = concave(x);
        PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1, M);
        P = repmat(2:2:2* M,size(P,1),1).*P;
    case 'WFG5'
        M      = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);
        
        z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
        
        t1 = zeros(N,K+L);
        t1 = s_decept(z01,0.35,0.001,0.05);
        
        t2 = zeros(N,M);
        for i = 1 : M-1
            t2(:,i) = r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t2(:,M) = r_sum(t1(:,K+1:K+L),ones(1,L));
        
        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
        end
        x(:,M) = t2(:,M);
        
        h = concave(x);
        PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1, M);
        P = repmat(2:2:2* M,size(P,1),1).*P;
    case 'WFG6'
        M      = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);
        
        z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
        
        t1 = zeros(N,K+L);
        t1(:,1:K)     = z01(:,1:K);
        t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);
        
        t2 = zeros(N,M);
        for i = 1 : M-1
            t2(:,i) = r_nonsep(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),K/(M-1));
        end
        % Same as <t2(:,M)=r_nonsep(t1(:,K+1:end),L)>
        SUM = zeros(N,1);
        for i = K+1 : K+L-1
            for j = i+1 : K+L
                SUM = SUM + abs(t1(:,i)-t1(:,j));
            end
        end
        t2(:,M) = (sum(t1(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));
        % -------------------------------------------
        
        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
        end
        x(:,M) = t2(:,M);
        
        h = concave(x);
        PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1, M);
        P = repmat(2:2:2* M,size(P,1),1).*P;
    case 'WFG7'
        M      = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);
        
        z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
        
        t1 = zeros(N,K+L);
        % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,i+1:end),ones(1,K+L-i)),0.98/49.98,0.02,50)>
        Y = (fliplr(cumsum(fliplr(z01),2))-z01)./repmat(K+L-1:-1:0,N,1);
        t1(:,1:K) = z01(:,1:K).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,1:K)).*abs(floor(0.5-Y(:,1:K))+0.98/49.98)));
        % ------------------------------------------------------------------------------------------
        t1(:,K+1:end) = z01(:,K+1:end);
        
        t2 = zeros(N,K+L);
        t2(:,1:K)     = t1(:,1:K);
        t2(:,K+1:end) = s_linear(t1(:,K+1:end),0.35);
        
        t3 = zeros(N,M);
        for i = 1 : M-1
            t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t3(:,M) = r_sum(t2(:,K+1:K+L),ones(1,L));
        
        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
        end
        x(:,M) = t3(:,M);
        
        h = concave(x);
        PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1, M);
        P = repmat(2:2:2* M,size(P,1),1).*P;
    case 'WFG8'
        M      = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);
        
        z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
        
        t1 = zeros(N,K+L);
        t1(:,1:K) = z01(:,1:K);
        % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,1:i-1),ones(1,i-1)),0.98/49.98,0.02,50)>
        Y = (cumsum(z01,2)-z01)./repmat(0:K+L-1,N,1);
        t1(:,K+1:K+L) = z01(:,K+1:K+L).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,K+1:K+L)).*abs(floor(0.5-Y(:,K+1:K+L))+0.98/49.98)));
        % --------------------------------------------------------------------------------------
        
        t2 = zeros(N,K+L);
        t2(:,1:K)     = t1(:,1:K);
        t2(:,K+1:end) = s_linear(t1(:,K+1:end),0.35);
        
        t3 = zeros(N,M);
        for i = 1 : M-1
            t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t3(:,M) = r_sum(t2(:,K+1:K+L),ones(1,L));
        
        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
        end
        x(:,M) = t3(:,M);
        
        h = concave(x);
        PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1, M);
        P = repmat(2:2:2* M,size(P,1),1).*P;
    case 'WFG9'
        M      = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);
        
        z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
        
        t1 = zeros(N,K+L);
        % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,i+1:end),ones(1,K+L-i)),0.98/49.98,0.02,50)>
        Y = (fliplr(cumsum(fliplr(z01),2))-z01)./repmat(K+L-1:-1:0,N,1);
        t1(:,1:K+L-1) = z01(:,1:K+L-1).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,1:K+L-1)).*abs(floor(0.5-Y(:,1:K+L-1))+0.98/49.98)));
        % ------------------------------------------------------------------------------------------
        t1(:,end)     = z01(:,end);
        
        t2 = zeros(N,K+L);
        t2(:,1:K)     = s_decept(t1(:,1:K),0.35,0.001,0.05);
        t2(:,K+1:end) = s_multi(t1(:,K+1:end),30,95,0.35);
        
        t3 = zeros(N,M);
        for i = 1 : M-1
            t3(:,i) = r_nonsep(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),K/(M-1));
        end
        % Same as <t3(:,M)=r_nonsep(t2(:,K+1:end),L)>
        SUM = zeros(N,1);
        for i = K+1 : K+L-1
            for j = i+1 : K+L
                SUM = SUM + abs(t2(:,i)-t2(:,j));
            end
        end
        t3(:,M) = (sum(t2(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));
        % -------------------------------------------
        
        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
        end
        x(:,M) = t3(:,M);
        
        h = concave(x);
        PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1, M);
        P = repmat(2:2:2* M,size(P,1),1).*P;
    case 'UF1'
        J1 = 3 : 2 : D;
        J2 = 2 : 2 : D;
        Y  = PopDec - sin(6*pi*repmat(PopDec(:,1),1,D)+repmat(1:D,size(PopDec,1),1)*pi/D);
        PopObj(:,1) = PopDec(:,1)         + 2*mean(Y(:,J1).^2,2);
        PopObj(:,2) = 1-sqrt(PopDec(:,1)) + 2*mean(Y(:,J2).^2,2);
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5;
    case 'UF2'
        J1 = 3 : 2 : D;
        J2 = 2 : 2 : D;
        X=PopDec;
        Y       = zeros(size(X));
        X1      = repmat(X(:,1),1,length(J1));
        Y(:,J1) = X(:,J1)-(0.3*X1.^2.*cos(24*pi*X1+4*repmat(J1,size(X,1),1)*pi/D)+0.6*X1).*cos(6*pi*X1+repmat(J1,size(X,1),1)*pi/D);
        X1      = repmat(X(:,1),1,length(J2));
        Y(:,J2) = X(:,J2)-(0.3*X1.^2.*cos(24*pi*X1+4*repmat(J2,size(X,1),1)*pi/D)+0.6*X1).*sin(6*pi*X1+repmat(J2,size(X,1),1)*pi/D);
        PopObj(:,1) = X(:,1)         + 2*mean(Y(:,J1).^2,2);
        PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean(Y(:,J2).^2,2);
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5;
    case 'UF3'
        X=PopDec;
        J1 = 3 : 2 : D;
        J2 = 2 : 2 : D;
        Y  = X - repmat(X(:,1),1,D).^(0.5*(1+3*(repmat(1:D,size(X,1),1)-2)/(D-2)));
        PopObj(:,1) = X(:,1)         + 2/length(J1)*(4*sum(Y(:,J1).^2,2)-2*prod(cos(20*Y(:,J1)*pi./sqrt(repmat(J1,size(X,1),1))),2)+2);
        PopObj(:,2) = 1-sqrt(X(:,1)) + 2/length(J2)*(4*sum(Y(:,J2).^2,2)-2*prod(cos(20*Y(:,J2)*pi./sqrt(repmat(J2,size(X,1),1))),2)+2);
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5;
    case 'UF4'
        X=PopDec;
        J1 = 3 : 2 : D;
        J2 = 2 : 2 : D;
        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        hY = abs(Y)./(1+exp(2*abs(Y)));
        PopObj(:,1) = X(:,1)      + 2*mean(hY(:,J1),2);
        PopObj(:,2) = 1-X(:,1).^2 + 2*mean(hY(:,J2),2);
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^2;
    case 'UF5'
        X=PopDec;
        J1 = 3 : 2 : D;
        J2 = 2 : 2 : D;
        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        hY = 2*Y.^2 - cos(4*pi*Y) + 1;
        PopObj(:,1) = X(:,1)   + (1/20+0.1)*abs(sin(20*pi*X(:,1)))+2*mean(hY(:,J1),2);
        PopObj(:,2) = 1-X(:,1) + (1/20+0.1)*abs(sin(20*pi*X(:,1)))+2*mean(hY(:,J2),2);
        P(:,1) = (0:1:20)'/20;
        P(:,2) = 1 - P(:,1);
    case 'UF6'
        X=PopDec;
        J1 = 3 : 2 : D;
        J2 = 2 : 2 : D;
        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        PopObj(:,1) = X(:,1)   + max(0,2*(1/4+0.1)*sin(4*pi*X(:,1)))+2/length(J1)*(4*sum(Y(:,J1).^2,2)-2*prod(cos(20*Y(:,J1)*pi./sqrt(repmat(J1,size(X,1),1))),2)+2);
        PopObj(:,2) = 1-X(:,1) + max(0,2*(1/4+0.1)*sin(4*pi*X(:,1)))+2/length(J2)*(4*sum(Y(:,J2).^2,2)-2*prod(cos(20*Y(:,J2)*pi./sqrt(repmat(J2,size(X,1),1))),2)+2);
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1);
        P(P(:,1)>0 & P(:,1)<1/4 | P(:,1)>1/2 & P(:,1)<3/4,:) = [];
    case 'UF7'
        X=PopDec;
        J1 = 3 : 2 : D;
        J2 = 2 : 2 : D;
        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        PopObj(:,1) = X(:,1).^0.2   + 2*mean(Y(:,J1).^2,2);
        PopObj(:,2) = 1-X(:,1).^0.2 + 2*mean(Y(:,J2).^2,2);
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1);
    case 'UF8'
        X=PopDec;
        J1 = 4 : 3 : D;
        J2 = 5 : 3 : D;
        J3 = 3 : 3 : D;
        Y  = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean(Y(:,J1).^2,2);
        PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean(Y(:,J2).^2,2);
        PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean(Y(:,J3).^2,2);
        P = UniformPoint(N,3);
        P = P./repmat(sqrt(sum(P.^2,2)),1,3);
    case 'UF9'
        X=PopDec;
        J1 = 4 : 3 : D;
        J2 = 5 : 3 : D;
        J3 = 3 : 3 : D;
        Y  = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        PopObj(:,1) = 0.5*(max(0,1.1*(1-4*(2*X(:,1)-1).^2))+2*X(:,1)).*X(:,2)   + 2*mean(Y(:,J1).^2,2);
        PopObj(:,2) = 0.5*(max(0,1.1*(1-4*(2*X(:,1)-1).^2))-2*X(:,1)+2).*X(:,2) + 2*mean(Y(:,J2).^2,2);
        PopObj(:,3) = 1-X(:,2)                                                  + 2*mean(Y(:,J3).^2,2);
        P = UniformPoint(N,3);
        P(P(:,1)>(1-P(:,3))/4 & P(:,1)<(1-P(:,3))*3/4,:) = [];
    case 'UF10'
        X=PopDec;
        J1 = 4 : 3 : D;
        J2 = 5 : 3 : D;
        J3 = 3 : 3 : D;
        Y  = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        Y  = 4*Y.^2 - cos(8*pi*Y) + 1;
        PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean(Y(:,J1),2);
        PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean(Y(:,J2),2);
        PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean(Y(:,J3),2);
        P = UniformPoint(N,3);
        P = P./repmat(sqrt(sum(P.^2,2)),1,3);
        
    case 'MaF1'
        M      = numObj;
        g      = sum((PopDec(:,M:end)-0.5).^2,2);
        PopObj = repmat(1+g,1,M) - repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),PopDec(:,1:M-1)],2)).*[ones(size(g,1),1),1-PopDec(:,M-1:-1:1)];
        P= 1 - UniformPoint(N,M);
    case 'MaF2'
        M      = numObj;
        g = zeros(size(PopDec,1),M);
        for m = 1 : M
            if m < M
                g(:,m) = sum(((PopDec(:,M+(m-1)*floor((D-M+1)/M):M+m*floor((D-M+1)/M)-1)/2+1/4)-0.5).^2,2);
            else
                g(:,m) = sum(((PopDec(:,M+(M-1)*floor((D-M+1)/M):D)/2+1/4)-0.5).^2,2);
            end
        end
        PopObj = (1+g).*fliplr(cumprod([ones(size(g,1),1),cos((PopDec(:,1:M-1)/2+1/4)*pi/2)],2)).*[ones(size(g,1),1),sin((PopDec(:,M-1:-1:1)/2+1/4)*pi/2)];
        P=GetOptimum(M,N);
    case 'MaF3'
        M      = numObj;
        g      = 100*(D-M+1+sum((PopDec(:,M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,M:end)-0.5)),2));
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        PopObj = [PopObj(:,1:M-1).^4,PopObj(:,M).^2];
        
        P    = UniformPoint(N,M).^2;
        temp = sum(sqrt(P(:,1:end-1)),2) + P(:,end);
        P    = P./[repmat(temp.^2,1,size(P,2)-1),temp];
    case 'MaF4'
        M      = numObj;
        g      = 100*(D-M+1+sum((PopDec(:,M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,M:end)-0.5)),2));
        PopObj = repmat(1+g,1,M) - repmat(1+g,1,M).*fliplr(cumprod([ones(size(PopDec,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(PopDec,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        PopObj = PopObj.*repmat(2.^(1:M),size(PopDec,1),1);
        
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        P = (1-P).*repmat(2.^(1:M),size(P,1),1);
    case 'MaF5'
        M      = numObj;
        PopDec(:,1:M-1) = PopDec(:,1:M-1).^100;
        g      = sum((PopDec(:,M:end)-0.5).^2,2);
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        PopObj = PopObj.*repmat(2.^(M:-1:1),size(g,1),1);
        
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        P = P.*repmat(2.^(M:-1:1),size(P,1),1);
    case 'MaF6'
        M      = numObj;
        I      = 2;
        g      = sum((PopDec(:,M:end)-0.5).^2,2);
        Temp   = repmat(g,1,M-I);
        PopDec(:,I:M-1) = (1+2*Temp.*PopDec(:,I:M-1))./(2+2*Temp);
        PopObj = repmat(1+100*g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        
        P = UniformPoint(N,I);
        P = P./repmat(sqrt(sum(P.^2,2)),1,size(P,2));
        P = [P(:,ones(1,M-size(P,2))),P];
        P = P./sqrt(2).^repmat(max([M-I,M-I:-1:2-I],0),size(P,1),1);
    case 'MaF8'
        M      = numObj;
        Points  = [];
        [thera,rho] = cart2pol(0,1);
        [Points(:,1),Points(:,2)] = pol2cart(thera-(1:M)*2*pi/M,rho);
        PopObj = pdist2(PopDec,Points);
        
        [X,Y] = ndgrid(linspace(-1,1,ceil(sqrt(N))));
        ND    = inpolygon(X(:),Y(:),Points(:,1),Points(:,2));
        P     = pdist2([X(ND),Y(ND)],Points);
    case 'MaF9'
        obj.M      = numObj;
        % Parameter setting
        obj.D        = 2;
        obj.lower    = [-10000,-10000];
        obj.upper    = [10000,10000];
        obj.encoding = 'real';
        % Generate vertexes
        obj.Points = [];
        [thera,rho] = cart2pol(0,1);
        [obj.Points(:,1),obj.Points(:,2)] = pol2cart(thera-(1:obj.M)*2*pi/obj.M,rho);
        % Generate invalid polygons
        head = repmat((1:obj.M)',ceil(obj.M/2-2),1);
        tail = repmat(1:ceil(obj.M/2-2),obj.M,1);
        tail = head + tail(:);
        obj.Polygons = cell(1,length(head));
        for i = 1 : length(obj.Polygons)
            obj.Polygons{i} = obj.Points(mod((head(i):tail(i))-1,obj.M)+1,:);
            obj.Polygons{i} = [obj.Polygons{i};repmat(2*Intersection(obj.Points(mod([head(i)-1,head(i),tail(i),tail(i)+1]-1,obj.M)+1,:)),size(obj.Polygons{i},1),1)-obj.Polygons{i}];
        end
        %% Repair invalid solutions
        PopDec = CalDec(obj,PopDec);
        PopObj = CalObj(obj,PopDec);
        
        
        [X,Y] = ndgrid(linspace(-1,1,ceil(sqrt(N))));
        ND    = inpolygon(X(:),Y(:),obj.Points(:,1),obj.Points(:,2));
        P     = obj.CalObj([X(ND),Y(ND)]);
        
    case 'MaF13'
        M      = numObj;
        Y = PopDec - 2*repmat(PopDec(:,2),1,D).*sin(2*pi*repmat(PopDec(:,1),1,D)+repmat(1:D,N,1)*pi/D);
        PopObj(:,1) = sin(PopDec(:,1)*pi/2)+ 2*mean(Y(:,4:3:D).^2,2);
        PopObj(:,2) = cos(PopDec(:,1)*pi/2).*sin(PopDec(:,2)*pi/2) + 2*mean(PopDec(:,5:3:D).^2,2);
        PopObj(:,3) = cos(PopDec(:,1)*pi/2).*cos(PopDec(:,2)*pi/2) + 2*mean(PopDec(:,3:3:D).^2,2);
        PopObj(:,4:M) = repmat(PopObj(:,1).^2+PopObj(:,2).^10+PopObj(:,3).^10+2*mean(Y(:,4:D).^2,2),1,M-3);
        
        P = UniformPoint(N,3);
        P = P./repmat(sqrt(sum(P.^2,2)),1,3);
        P = [P,repmat(P(:,1).^2+P(:,2).^10+P(:,3).^10,1,M-3)];
    case 'LSMOP1'
        M      = numObj;
        PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
        G = zeros(N,M);
        nk = 5;
        % Calculate the number of variables in each subcomponent
        c = 3.8*0.1*(1-0.1);
        for i = 1 : M-1
            c = [c,3.8.*c(end).*(1-c(end))];
        end
        sublen = floor(c./sum(c).*(D-M+1)/nk);
        len    = [0,cumsum(sublen*nk)];
        for i = 1 : 2 : M
            for j = 1 : nk
                G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
            end
        end
        for i = 2 : 2 : M
            for j = 1 : nk
                G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
            end
        end
        G      = G./repmat(sublen,N,1)./nk;
        PopObj = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
        P = UniformPoint(N,M);
    case 'LSMOP2'
        M      = numObj;
        PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
        nk = 5;
        % Calculate the number of variables in each subcomponent
        c = 3.8*0.1*(1-0.1);
        for i = 1 : M-1
            c = [c,3.8.*c(end).*(1-c(end))];
        end
        sublen = floor(c./sum(c).*(D-M+1)/nk);
        len    = [0,cumsum(sublen*nk)];
        G = zeros(N,M);
        for i = 1 : 2 : M
            for j = 1 : nk
                G(:,i) = G(:,i) + Griewank(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
            end
        end
        for i = 2 : 2 : M
            for j = 1 : nk
                G(:,i) = G(:,i) + Schwefel(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
            end
        end
        G      = G./repmat(sublen,N,1)./nk;
        PopObj = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
        P = UniformPoint(N,M);
    case 'LSMOP3'
        M      = numObj;
        nk = 5;
        % Calculate the number of variables in each subcomponent
            c = 3.8*0.1*(1-0.1);
            for i = 1 : M-1
                c = [c,3.8.*c(end).*(1-c(end))];
            end
            sublen = floor(c./sum(c).*(D-M+1)/nk);
            len    = [0,cumsum(sublen*nk)];
            PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
            for i = 1 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Rastrigin(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Rosenbrock(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                end
            end
            G      = G./repmat(sublen,N,1)./nk;
            PopObj = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
            P= UniformPoint(N,M);
        
end%switch

end
function Output = s_linear(y,A)
Output = abs(y-A)./abs(floor(A-y)+A);
end
function Output = disc(x)%wfg2
Output = 1-x(:,1).*(cos(5*pi*x(:,1))).^2;
end
function Output = concave(x)%WFG4
Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end
function Output = s_multi(y,A,B,C)%WFG4
Output = (1+cos((4*A+2)*pi*(0.5-abs(y-C)/2./(floor(C-y)+C)))+4*B*(abs(y-C)/2./(floor(C-y)+C)).^2)/(B+2);
end
function Output = s_decept(y,A,B,C)%wfg5
Output = 1+(abs(y-A)-B).*(floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
end
function Output = r_nonsep(y,A)%wfg6
Output = zeros(size(y,1),1);
for j = 1 : size(y,2)
    Temp = zeros(size(y,1),1);
    for k = 0 : A-2
        Temp = Temp+abs(y(:,j)-y(:,1+mod(j+k,size(y,2))));
    end
    Output = Output+y(:,j)+Temp;
end
Output = Output./(size(y,2)/A)/ceil(A/2)/(1+2*A-2*ceil(A/2));
end
function Output = b_param(y,Y,A,B,C)
Output = y.^(B+(C-B)*(A-(1-2*Y).*abs(floor(0.5-Y)+A)));
end
function Output = b_flat(y,A,B,C)
Output = A+min(0,floor(y-B))*A.*(B-y)/B-min(0,floor(C-y))*(1-A).*(y-C)/(1-C);
Output = roundn(Output,-6);
end
function Output = b_poly(y,a)
Output = y.^a;
end

function Output = r_sum(y,w)
Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = convex(x)
Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = mixed(x)
Output = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
end

function Output = linear(x)
Output = fliplr(cumprod([ones(size(x,1),1),x(:,1:end-1)],2)).*[ones(size(x,1),1),1-x(:,end-1:-1:1)];
end


function PopObj = CalObjVNT1(PopDec)
PopObj(:,1) = PopDec(:,1).^2 + (PopDec(:,2)-1).^2;
PopObj(:,2) = PopDec(:,1).^2 + (PopDec(:,2)+1).^2 + 1;
PopObj(:,3) = (PopDec(:,1)-1).^2 + PopDec(:,2).^2 + 2;
end
function PopObj = CalObjVNT2(PopDec)
PopObj(:,1) = (PopDec(:,1)-2).^2/2 + (PopDec(:,2)+1).^2/13 + 3;
PopObj(:,2) = (PopDec(:,1)+PopDec(:,2)-3).^2/36 + (-PopDec(:,1)+PopDec(:,2)+2).^2/8 - 17;
PopObj(:,3) = (PopDec(:,1)+2*PopDec(:,2)-1).^2/175 + (2*PopDec(:,2)-PopDec(:,1)).^2/17 - 13;
end
function PopObj = CalObjVNT3(PopDec)
temp = PopDec(:,1).^2 + PopDec(:,2).^2;
PopObj(:,1) = 0.5*temp + sin(temp);
PopObj(:,2) = (3*PopDec(:,1)-2*PopDec(:,2)+4).^2/8 + (PopDec(:,1)-PopDec(:,2)+1).^2/27 + 15;
PopObj(:,3) = 1./(temp+1) - 1.1*exp(-temp);
end
function PopObj = CalObjVNT4(PopDec)
PopObj(:,1) = (PopDec(:,1)-2).^2/2 + (PopDec(:,2)+1).^2/13 + 3;
PopObj(:,2) = (PopDec(:,1)+PopDec(:,2)-3).^2/175 + (2*PopDec(:,2)-PopDec(:,1)).^2/17 - 13;
PopObj(:,3) = (3*PopDec(:,1)-2*PopDec(:,2)+4).^2/8 + (PopDec(:,1)-PopDec(:,2)+1).^2/27 + 15;
end
function W = ReplicatePoint(SampleNum,M)
if M > 1
    SampleNum = (ceil(SampleNum^(1/M)))^M;
    Gap       = 0:1/(SampleNum^(1/M)-1):1;
    eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
    eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
else
    W = (0:1/(SampleNum-1):1)';
end
end
%% MaF2
function R = GetOptimum(M,N)
R = UniformPoint(N,M);
c = zeros(size(R,1),M-1);
for i = 1 : size(R,1)
    for j = 2 : M
        temp = R(i,j)/R(i,1)*prod(c(i,M-j+2:M-1));
        c(i,M-j+1) = sqrt(1/(1+temp^2));
    end
end
if M > 5
    c = c.*(cos(pi/8)-cos(3*pi/8)) + cos(3*pi/8);
else
    c(any(c<cos(3*pi/8)|c>cos(pi/8),2),:) = [];
end
R = fliplr(cumprod([ones(size(c,1),1),c(:,1:M-1)],2)).*[ones(size(c,1),1),sqrt(1-c(:,M-1:-1:1).^2)];
end
%% MaF9
function PopDec = CalDec(obj,PopDec)
Invalid = getInvalid(PopDec,obj.Polygons,obj.Points);
while any(Invalid)
    PopDec(Invalid,:) = unifrnd(repmat(obj.lower,sum(Invalid),1),repmat(obj.upper,sum(Invalid),1));
    Invalid           = getInvalid(PopDec,obj.Polygons,obj.Points);
end
end
%% Calculate objective values
function PopObj = CalObj(obj,PopDec)
PopObj = zeros(size(PopDec,1),size(obj.Points,1));
for m = 1 : size(obj.Points,1)
    PopObj(:,m) = Point2Line(PopDec,obj.Points(mod(m-1:m,size(obj.Points,1))+1,:));
end
end
%% LSMOP1
function f = Sphere(x)
f = sum(x.^2,2);
end

function f = Griewank(x)
f = sum(x.^2,2)./4000 - prod(cos(x./repmat(sqrt(1:size(x,2)),size(x,1),1)),2) + 1;
end
%%  LSMOP2
function f = Schwefel(x)
f = max(abs(x),[],2);
end
%%  LSMOP3
function f = Rastrigin(x)
    f = sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end

function f = Rosenbrock(x)
    f = sum(100.*(x(:,1:size(x,2)-1).^2-x(:,2:size(x,2))).^2+(x(:,1:size(x,2)-1)-1).^2,2);
end