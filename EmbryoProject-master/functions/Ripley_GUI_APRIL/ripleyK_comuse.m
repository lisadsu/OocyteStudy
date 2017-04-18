function out=ripleyK_comuse(data,t, doplot,dowait)
% out=ripleyK_comuse(data,t,doplot,dowait)
% Calculate the Ripley K function for a distribution with edge correction in 3D
% Input:
% data - a matrix [X Y Z]
% t - the steps at which to evaluate the K-function. eg 2:1:12 or 1:28
% doplot - whether to plot the resulting K-function
% dowait - whether to show a bar visualizing the progress
% Output:
% a structure with the fields K (the K-function), EK (the expected K-function if this distribution was CSR),
% Correction (the correction terms used at the edges), Geom (a vector with the side lengths, total volume and number of points of the distribution).
%
% License: RipleyGUI is distributed free under the conditions that
% (1) it shall not be incorporated in software that is subsequently sold; 
% (2) the authorship of the software shall be acknowledged in any publication that uses results generated by the software; 
% (3) this notice shall remain in place in each source file. 
if (nargin < 3) || isempty(doplot)
    doplot = 1;
    dowait = 1;
elseif (nargin < 4)
    dowait = 1;
end
[n dim]=size(data);
tx=numel(t);
switch dim
    case 2
        s=(max(data)-min(data));
        a=s(1)*s(2);
        K=zeros(n,tx);
        EK=zeros(1,tx);C=zeros(1,tx);
        D=euclid(data,data);
        for i=1:tx
            EK(i)=pi*power(t(i),2);
            C(i)=1-(4/(3*pi))*((t(i)/s(1))+(t(i)/s(2)))+((11/(3*pi))-1)*((t(i)^2)/(s(1)*s(2)));
            for j=1:n
                K(j,i)=(a*(sum(sum(D(j,:)<=t(i))-1)))/(C(i)*n*(n-1));
            end
        end
        if min(C)>0
            ripK=sum(K);
        elseif min(C)<=0
            tx=find(C>0,1,'last');
            ripK=sum(K(:,1:tx));
            EK=EK(1:tx);
            warning('The edge correction term in 2D causes singularity. Try a smaller value for t.') %#ok<WNTAG>
        end
        if doplot %dont plot if doplot is 0
            subplot(2,2,[1 2])
            scatter(data(:,1),data(:,2),'.r'), drawnow; axis equal,grid on
            subplot(2,2,3)
            plot(ripK-EK,'-sb','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',2)
            drawnow; grid on,hold on,xlabel('t'),ylabel('K(t)-E[K(t)]')
            subplot(2,2,4)
            plot(ripK./EK,'o-r','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',2)
            drawnow; grid on,hold on,xlabel('t'),ylabel('K(t)/E[K(t)]')
        end
        out.EK=EK;
        out.K=ripK;
        out.Correction=C;
        out.Geom=[s(1) s(2) a n];
    case 3
        warning_state = warning('off','MATLAB:quad:MinStepSize');
        % Hide warning MinStepSize for a cleaner output.
        % The warning indicates that the integral calculation in quad is
        % approximative
        s=(max(data)-min(data));
        v=s(1)*s(2)*s(3);
        w=correct3D(data,t,0.1,dowait);
        if w == 0
            out=0;
            return;
        end
        K=zeros(n,tx);
        EK=zeros(1,tx);
        D=euclid(data,data);
        for j=1:n
            for i=1:tx
                EK(i)=(4*pi*power(t(i),3))/3;
                K(j,i)=(v*(sum(sum(D(j,:)<=t(i))-1)/w(j,i)))/n^2;
            end
        end
        ripK=sum(K);
        if doplot %dont plot if doplot is 0
            subplot(2,2,[1 2])
            scatter3(data(:,1),data(:,2),data(:,3),'.r'), drawnow; axis equal,grid on
            subplot(2,2,3)
            plot(ripK-EK,'-sb','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',2)
            drawnow; grid on,hold on,xlabel('t'),ylabel('K(t)-E[K(t)]')
            subplot(2,2,4)
            plot(ripK./EK,'o-r','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',2)
            drawnow; grid on,hold on,xlabel('t'),ylabel('K(t)/E[K(t)]')
        end
        out.EK=EK;
        out.K=ripK;
        out.Correction=w;
        out.Geom=[s(1) s(2) s(3) v n];
        warning(warning_state);
        % Reset the warning state to normal
end
end

% function out=correct3D(data,t,tol,bound)
% bound=[minx miny minz;maxx maxy maxz];
function out=correct3D(data,t,tol,dowait,bound)
n=length(data);
if nargin<5
    nx=min(data(:,1));xx=max(data(:,1));        % n=min x=max
    ny=min(data(:,2));xy=max(data(:,2));
    nz=min(data(:,3));xz=max(data(:,3));
    bound=[nx ny nz;xx xy xz];
end
if nargin<4
    dowait = 0;
end
if dowait
    wait_h = waitbar(0, 'Doing Edge Correction...');
end
for j=1:n
    if dowait
        try
            waitbar(j/n,wait_h);
        catch
            helpdlg('Calculations had been aborted','Waitbar message')
            out=0;
            return;
        end
    end
    for k=1:numel(t)
        if t(k)>euclid(data(j,1),bound(1,1))
            hlx=t(k)-euclid(data(j,1),bound(1,1));
            baslx=sqrt(t(k)^2-euclid(data(j,1),bound(1,1))^2);
            clx=(pi*hlx/6)*((3*baslx^2)+(hlx^2));
        else clx=0;
        end
        if t(k)>euclid(data(j,1),bound(2,1))
            hux=t(k)-euclid(data(j,1),bound(2,1));
            basux=sqrt(t(k)^2-euclid(data(j,1),bound(2,1))^2);
            cux=(pi*hux/6)*((3*basux^2)+(hux^2));
        else cux=0;
        end
        cx(j,k)=clx+cux;
        if t(k)>euclid(data(j,2),bound(1,2))
            hly=t(k)-euclid(data(j,2),bound(1,2));
            basly=sqrt(t(k)^2-euclid(data(j,2),bound(1,2))^2);
            cly=(pi*hly/6)*((3*basly^2)+(hly^2));
        else cly=0;
        end
        if t(k)>euclid(data(j,2),bound(2,2))
            huy=t(k)-euclid(data(j,2),bound(2,2));
            basuy=sqrt(t(k)^2-euclid(data(j,2),bound(2,2))^2);
            cuy=(pi*huy/6)*((3*basuy^2)+(huy^2));
        else cuy=0;
        end
        cy(j,k)=cly+cuy;
        if t(k)>euclid(data(j,3),bound(1,3))
            hlz=t(k)-euclid(data(j,3),bound(1,3));
            baslz=sqrt(t(k)^2-euclid(data(j,3),bound(1,3))^2);
            clz=(pi*hlz/6)*((3*baslz^2)+(hlz^2));
        else clz=0;
        end
        if t(k)>euclid(data(j,3),bound(2,3))
            huz=t(k)-euclid(data(j,3),bound(2,3));
            basuz=sqrt(t(k)^2-euclid(data(j,3),bound(2,3))^2);
            cuz=(pi*huz/6)*((3*basuz^2)+(huz^2));
        else cuz=0;
        end
        cz(j,k)=clz+cuz;
        if and(t(k)>euclid(data(j,1),bound(1,1)),t(k)>euclid(data(j,2),bound(1,2)))
            a1=euclid(data(j,1),bound(1,1));
            b1=euclid(data(j,2),bound(1,2));
            intg1=sqrt(t(k)^2-b1^2-a1^2)*(t(k)^2-b1^2-a1^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b1).*(t(k).^2-x.^2-y.^2>=b1^2);
            q1=dblquad(F,-intg1,intg1,a1,sqrt(t(k)^2-b1^2),tol);
        else q1=0;
        end
        if and(t(k)>euclid(data(j,1),bound(2,1)),t(k)>euclid(data(j,2),bound(1,2)))
            a2=euclid(data(j,1),bound(2,1));
            b2=euclid(data(j,2),bound(1,2));
            intg2=sqrt(t(k)^2-b2^2-a2^2)*(t(k)^2-b2^2-a2^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b2).*(t(k).^2-x.^2-y.^2>=b2^2);
            q2=dblquad(F,-intg2,intg2,a2,sqrt(t(k)^2-b2^2),tol);
        else q2=0;
        end
        if and(t(k)>euclid(data(j,1),bound(1,1)),t(k)>euclid(data(j,2),bound(2,2)))
            a3=euclid(data(j,1),bound(1,1));
            b3=euclid(data(j,2),bound(2,2));
            intg3=sqrt(t(k)^2-b3^2-a3^2)*(t(k)^2-b3^2-a3^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b3).*(t(k).^2-x.^2-y.^2>=b3^2);
            q3=dblquad(F,-intg3,intg3,a3,sqrt(t(k)^2-b3^2),tol);
        else q3=0;
        end
        if and(t(k)>euclid(data(j,1),bound(2,1)),t(k)>euclid(data(j,2),bound(2,2)))
            a4=euclid(data(j,1),bound(2,1));
            b4=euclid(data(j,2),bound(2,2));
            intg4=sqrt(t(k)^2-b4^2-a4^2)*(t(k)^2-b4^2-a4^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b4).*(t(k).^2-x.^2-y.^2>=b4^2);
            q4=dblquad(F,-intg4,intg4,a4,sqrt(t(k)^2-b4^2),tol);
        else q4=0;
        end
        if and(t(k)>euclid(data(j,1),bound(1,1)),t(k)>euclid(data(j,3),bound(1,3)))
            a5=euclid(data(j,1),bound(1,1));
            b5=euclid(data(j,3),bound(1,3));
            intg5=sqrt(t(k)^2-b5^2-a5^2)*(t(k)^2-b5^2-a5^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b5).*(t(k).^2-x.^2-y.^2>=b5^2);
            q5=dblquad(F,-intg5,intg5,a5,sqrt(t(k)^2-b5^2),tol);
        else q5=0;
        end
        if and(t(k)>euclid(data(j,1),bound(2,1)),t(k)>euclid(data(j,3),bound(1,3)))
            a6=euclid(data(j,1),bound(2,1));
            b6=euclid(data(j,3),bound(1,3));
            intg6=sqrt(t(k)^2-b6^2-a6^2)*(t(k)^2-b6^2-a6^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b6).*(t(k).^2-x.^2-y.^2>=b6^2);
            q6=dblquad(F,-intg6,intg6,a6,sqrt(t(k)^2-b6^2),tol);
        else q6=0;
        end
        if and(t(k)>euclid(data(j,1),bound(1,1)),t(k)>euclid(data(j,3),bound(2,3)))
            a7=euclid(data(j,1),bound(1,1));
            b7=euclid(data(j,3),bound(2,3));
            intg7=sqrt(t(k)^2-b7^2-a7^2)*(t(k)^2-b7^2-a7^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b7).*(t(k).^2-x.^2-y.^2>=b7^2);
            q7=dblquad(F,-intg7,intg7,a7,sqrt(t(k)^2-b7^2),tol);
        else q7=0;
        end
        if and(t(k)>euclid(data(j,1),bound(2,1)),t(k)>euclid(data(j,3),bound(2,3)))
            a8=euclid(data(j,1),bound(2,1));
            b8=euclid(data(j,3),bound(2,3));
            intg8=sqrt(t(k)^2-b8^2-a8^2)*(t(k)^2-b8^2-a8^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b8).*(t(k).^2-x.^2-y.^2>=b8^2);
            q8=dblquad(F,-intg8,intg8,a8,sqrt(t(k)^2-b8^2),tol);
        else q8=0;
        end
        if and(t(k)>euclid(data(j,2),bound(1,2)),t(k)>euclid(data(j,3),bound(1,3)))
            a9=euclid(data(j,2),bound(1,2));
            b9=euclid(data(j,3),bound(1,3));
            intg9=sqrt(t(k)^2-b9^2-a9^2)*(t(k)^2-b9^2-a9^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b9).*(t(k).^2-x.^2-y.^2>=b9^2);
            q9=dblquad(F,-intg9,intg9,a9,sqrt(t(k)^2-b9^2),tol);
        else q9=0;
        end
        if and(t(k)>euclid(data(j,2),bound(2,2)),t(k)>euclid(data(j,3),bound(1,3)))
            a10=euclid(data(j,2),bound(2,2));
            b10=euclid(data(j,3),bound(1,3));
            intg10=sqrt(t(k)^2-b10^2-a10^2)*(t(k)^2-b10^2-a10^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b10).*(t(k).^2-x.^2-y.^2>=b10^2);
            q10=dblquad(F,-intg10,intg10,a10,sqrt(t(k)^2-b10^2),tol);
        else q10=0;
        end
        if and(t(k)>euclid(data(j,2),bound(1,2)),t(k)>euclid(data(j,3),bound(2,3)))
            a11=euclid(data(j,2),bound(1,2));
            b11=euclid(data(j,3),bound(2,3));
            intg11=sqrt(t(k)^2-b11^2-a11^2)*(t(k)^2-b11^2-a11^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b11).*(t(k).^2-x.^2-y.^2>=b11^2);
            q11=dblquad(F,-intg11,intg11,a11,sqrt(t(k)^2-b11^2),tol);
        else q11=0;
        end
        if and(t(k)>euclid(data(j,2),bound(2,2)),t(k)>euclid(data(j,3),bound(2,3)))
            a12=euclid(data(j,2),bound(2,2));
            b12=euclid(data(j,3),bound(2,3));
            intg12=sqrt(t(k)^2-b12^2-a12^2)*(t(k)^2-b12^2-a12^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-b12).*(t(k).^2-x.^2-y.^2>=b12^2);
            q12=dblquad(F,-intg12,intg12,a12,sqrt(t(k)^2-b12^2),tol);
        else q12=0;
        end
        q(j,k)=q1+q2+q3+q4+q5+q6+q7+q8+q9+q10+q11+q12;
        if (t(k)>euclid(data(j,1),bound(1,1)) && t(k)>euclid(data(j,2),bound(1,2)) && t(k)>euclid(data(j,3),bound(1,3)))
            as1=euclid(data(j,1),bound(1,1));
            bs1=euclid(data(j,2),bound(1,2));
            cs1=euclid(data(j,3),bound(1,3));
            intgs1=sqrt(t(k)^2-bs1^2-as1^2)*(t(k)^2-bs1^2-as1^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-bs1).*(t(k).^2-x.^2-y.^2>=bs1^2);
            sw1=dblquad(F,cs1,intgs1,as1,sqrt(t(k)^2-bs1^2),tol);
        else sw1=0;
        end
        if (t(k)>euclid(data(j,1),bound(1,1)) && t(k)>euclid(data(j,2),bound(1,2)) && t(k)>euclid(data(j,3),bound(2,3)))
            as2=euclid(data(j,1),bound(1,1));
            bs2=euclid(data(j,2),bound(1,2));
            cs2=euclid(data(j,3),bound(2,3));
            intgs2=sqrt(t(k)^2-bs2^2-as2^2)*(t(k)^2-bs2^2-as2^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-bs2).*(t(k).^2-x.^2-y.^2>=bs2^2);
            sw2=dblquad(F,-intgs2,-cs2,as2,sqrt(t(k)^2-bs2^2),tol);
        else sw2=0;
        end
        if (t(k)>euclid(data(j,1),bound(1,1)) && t(k)>euclid(data(j,2),bound(2,2)) && t(k)>euclid(data(j,3),bound(1,3)))
            as3=euclid(data(j,1),bound(1,1));
            bs3=euclid(data(j,2),bound(2,2));
            cs3=euclid(data(j,3),bound(1,3));
            intgs3=sqrt(t(k)^2-bs3^2-as3^2)*(t(k)^2-bs3^2-as3^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-bs3).*(t(k).^2-x.^2-y.^2>=bs3^2);
            sw3=dblquad(F,cs3,intgs3,as3,sqrt(t(k)^2-bs3^2),tol);
        else sw3=0;
        end
        if (t(k)>euclid(data(j,1),bound(1,1)) && t(k)>euclid(data(j,2),bound(2,2)) && t(k)>euclid(data(j,3),bound(2,3)))
            as4=euclid(data(j,1),bound(1,1));
            bs4=euclid(data(j,2),bound(2,2));
            cs4=euclid(data(j,3),bound(2,3));
            intgs4=sqrt(t(k)^2-bs4^2-as4^2)*(t(k)^2-bs4^2-as4^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-bs4).*(t(k).^2-x.^2-y.^2>=bs4^2);
            sw4=dblquad(F,-intgs4,-cs4,as4,sqrt(t(k)^2-bs4^2),tol);
        else sw4=0;
        end
        if (t(k)>euclid(data(j,1),bound(2,1)) && t(k)>euclid(data(j,2),bound(1,2)) && t(k)>euclid(data(j,3),bound(1,3)))
            as5=euclid(data(j,1),bound(2,1));
            bs5=euclid(data(j,2),bound(1,2));
            cs5=euclid(data(j,3),bound(1,3));
            intgs5=sqrt(t(k)^2-bs5^2-as5^2)*(t(k)^2-bs5^2-as5^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-bs5).*(t(k).^2-x.^2-y.^2>=bs5^2);
            sw5=dblquad(F,cs5,intgs5,as5,sqrt(t(k)^2-bs5^2),tol);
        else sw5=0;
        end
        if (t(k)>euclid(data(j,1),bound(2,1)) && t(k)>euclid(data(j,2),bound(1,2)) && t(k)>euclid(data(j,3),bound(2,3)))
            as6=euclid(data(j,1),bound(2,1));
            bs6=euclid(data(j,2),bound(1,2));
            cs6=euclid(data(j,3),bound(2,3));
            intgs6=sqrt(t(k)^2-bs6^2-as6^2)*(t(k)^2-bs6^2-as6^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-bs6).*(t(k).^2-x.^2-y.^2>=bs6^2);
            sw6=dblquad(F,-intgs6,-cs6,as6,sqrt(t(k)^2-bs6^2),tol);
        else sw6=0;
        end
        if (t(k)>euclid(data(j,1),bound(2,1)) && t(k)>euclid(data(j,2),bound(2,2)) && t(k)>euclid(data(j,3),bound(1,3)))
            as7=euclid(data(j,1),bound(2,1));
            bs7=euclid(data(j,2),bound(2,2));
            cs7=euclid(data(j,3),bound(1,3));
            intgs7=sqrt(t(k)^2-bs7^2-as7^2)*(t(k)^2-bs7^2-as7^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-bs7).*(t(k).^2-x.^2-y.^2>=bs7^2);
            sw7=dblquad(F,cs7,intgs7,as7,sqrt(t(k)^2-bs7^2),tol);
        else sw7=0;
        end
        if (t(k)>euclid(data(j,1),bound(2,1)) && t(k)>euclid(data(j,2),bound(2,2)) && t(k)>euclid(data(j,3),bound(2,3)))
            as8=euclid(data(j,1),bound(2,1));
            bs8=euclid(data(j,2),bound(2,2));
            cs8=euclid(data(j,3),bound(2,3));
            intgs8=sqrt(t(k)^2-bs8^2-as8^2)*(t(k)^2-bs8^2-as8^2>0);
            F=@(x,y)(sqrt(t(k).^2-(x.^2+y.^2))-bs8).*(t(k).^2-x.^2-y.^2>=bs8^2);
            sw8=dblquad(F,-intgs8,-cs8,as8,sqrt(t(k)^2-bs8^2),tol);
        else sw8=0;
        end
        sw(j,k)=sw1+sw2+sw3+sw4+sw5+sw6+sw7+sw8;
    end
end
if dowait
    close(wait_h)
end
for j=1:n
    for k=1:numel(t)
        w(j,k)=1-((cx(j,k)+cy(j,k)+cz(j,k)-q(j,k)+sw(j,k))/(4*pi*power(t(k),3)/3));
    end
end
if t(1)==0
    w(:,1)=1;
end
out=w;
end

function d=euclid(x,y)
d=sqrt(sum(x.^2,2)*ones(1,size(y,1))+ones(size(x,1),1)*sum(y.^2,2)'-2*(x*y'));end