%******************************************************************************
% This programe is implemented via MATLAB 2018.                              *
% Author :  Ping Wang                                                        *
% Contact:  pingwangsky@gmail.com; pingwangsky@163.com                       *
% License:  Copyright (c) 2020 Ping Wang, All rights reserved.               *
% My site:  https://sites.google.com/view/ping-wang-homepage                 *
%*****************************************************************************/
function [Rc,tc]=RDLT_GN1(Xw,xxn)
    % Input checks
    if (length(Xw) < 4)
        error('At least  4  points are supplied.');
    end
    
    if rank(Xw)>2
        % ordinar case;
        [Rc,tc]=RDLT_Ordinary1(Xw,xxn);
    else
        % planar case;
        [Rc,tc]=RDLT_Planar1(Xw,xxn);
    end

    n=length(Xw);
    %refine solution by using the damped Gauss-Newton method;
    Pu=ones(3,1)*xxn(1,:).*Xw;    Pv=ones(3,1)*xxn(2,:).*Xw;
    u=xxn(1,:); v=xxn(2,:);
    Pi=sum(Xw,2)/n; Pui=sum(Pu,2)/n; Pvi=sum(Pv,2)/n;
    ui=sum(u)/n; vi=sum(v)/n;
    P_mean=Xw-Pi*ones(1,n);
    Pu_mean=Pu-Pui*ones(1,n);
    Pv_mean=Pv-Pvi*ones(1,n);
    u_mean=u-ui*ones(1,n);
    v_mean=v-vi*ones(1,n);
    Px=P_mean(1,:).'; Py=P_mean(2,:).'; Pz=P_mean(3,:).';
    Pux=Pu_mean(1,:).';   Puy=Pu_mean(2,:).';   Puz=Pu_mean(3,:).';
    Pvx=Pv_mean(1,:).';   Pvy=Pv_mean(2,:).';   Pvz=Pv_mean(3,:).';
    E=[-u_mean.',Px-Puz,-2*Puy,2*Pz+2*Pux,-2*Py,Px+Puz,2*Py,2*Pz-2*Pux,-Px+Puz,-2*Puy,-Px-Puz;
       -v_mean.',Py-Pvz,-2*Pz-2*Pvy,2*Pvx,2*Px,-Py+Pvz,2*Px,-2*Pvx,Py+Pvz,2*Pz-2*Pvy,-Py-Pvz];
    G=E.'*E;
    solution = RefineGaussNewton(matrix2quaternion(Rc),E,G,tc);
    a=solution(1);
    b=solution(2);
    c=solution(3);
    d=solution(4);
    Rw=[a^2+b^2-c^2-d^2,2*b*c-2*a*d,2*b*d+2*a*c;
        2*b*c+2*a*d,a^2-b^2+c^2-d^2,2*c*d-2*a*b;
        2*b*d-2*a*c,2*c*d+2*a*b,a^2-b^2-c^2+d^2];
    factor=1/(a^2+b^2+c^2+d^2);
    Rc=Rw*factor;
    tc=[-Rw(1,:)*Pi+Rw(3,:)*Pui+ui;-Rw(2,:)*Pi+Rw(3,:)*Pvi+vi;1]*factor;
%     D=zeros(2*n,3);
%     E=zeros(2*n,1);
%     D(1:2:end,:)=[ones(n,1),zeros(n,1),-xxn(1,:).'];
%     D(2:2:end,:)=[zeros(n,1),ones(n,1),-xxn(2,:).'];
%     E(1:2:end,:)=(Rc(3,1:2)*Xw(1:2,:).*xxn(1,:)-Rc(1,1:2)*Xw(1:2,:)).';
%     E(2:2:end,:)=(Rc(3,1:2)*Xw(1:2,:).*xxn(2,:)-Rc(2,1:2)*Xw(1:2,:)).';
%     tc=((D.'*D)\D.')*E;
end

function [Rc,tc]=RDLT_Ordinary1(Xw,xxn)
n=length(Xw);

index=nchoosek(1:n,2)';
W1=Xw(:,index(1,:));
W2=Xw(:,index(2,:));
x1=xxn(:,index(1,:));
x2=xxn(:,index(2,:));

x1=[x1;ones(1,length(x1))];
x2=[x2;ones(1,length(x2))];

nc=cross(x1,x2);    Lw=W2-W1;   mw=cross(W1,W2);

ui=xxn(1,:).';   vi=xxn(2,:).';
Xi=Xw(1,:).';   Yi=Xw(2,:).';   Zi=Xw(3,:).';
A=[Xi,zeros(n,1),-Xi.*ui,Yi,zeros(n,1),-Yi.*ui,Zi,zeros(n,1),-Zi.*ui,ones(n,1),zeros(n,1),-ui;
   zeros(n,1),Xi,-Xi.*vi,zeros(n,1),Yi,-Yi.*vi,zeros(n,1),Zi,-Zi.*vi,zeros(n,1),ones(n,1),-vi];

m=length(W1);
nx=[nc(1,:),nc(1,:)].'; ny=[nc(2,:),nc(2,:)].'; nz=[nc(3,:),nc(3,:)].';
Pt=[W1,W2];
Xi=Pt(1,:).';   Yi=Pt(2,:).';   Zi=Pt(3,:).';
B=[Xi.*nx,Xi.*ny,Xi.*nz,Yi.*nx,Yi.*ny,Yi.*nz,Zi.*nx,Zi.*ny,Zi.*nz,nx,ny,nz];

mx=mw(1,:).';   my=mw(2,:).';   mz=mw(3,:).';
lx=Lw(1,:).';   ly=Lw(2,:).';   lz=Lw(3,:).';
nx=nc(1,:).';   ny=nc(2,:).';   nz=nc(3,:).';
C=[mx.*nz,zeros(m,1),-mx.*nx,my.*nz,zeros(m,1),-my.*nx,mz.*nz,zeros(m,1),-mz.*nx,lx.*nz,zeros(m,1),-lx.*nx,ly.*nz,zeros(m,1),-ly.*nx,lz.*nz,zeros(m,1),-lz.*nx;
   zeros(m,1),mx.*nz,-mx.*ny,zeros(m,1),my.*nz,-my.*ny,zeros(m,1),mz.*nz,-mz.*ny,zeros(m,1),lx.*nz,-lx.*ny,zeros(m,1),ly.*nz,-ly.*ny,zeros(m,1),lz.*nz,-lz.*ny];    
M = [ A(:,1:9), zeros(2*n, 9),A(:,10:12);
      B(:,1:9), zeros(2*m, 9),B(:,10:12);
	  C,zeros(2*m, 3)];
w=-M(:,1:end-1)\M(:,end);
H=reshape(w(1:9),3,3);
Rr= SVDdecomposition(H*Xw,Xw);

D=zeros(2*n,3);
E=zeros(2*n,1);
D(1:2:end,:)=[ones(n,1),zeros(n,1),-xxn(1,:).'];
D(2:2:end,:)=[zeros(n,1),ones(n,1),-xxn(2,:).'];
E(1:2:end,:)=(Rr(3,:)*Xw.*xxn(1,:)-Rr(1,:)*Xw).';
E(2:2:end,:)=(Rr(3,:)*Xw.*xxn(2,:)-Rr(2,:)*Xw).';

Rc=Rr;
tc=((D.'*D)\D.')*E;
end

function [Rc,tc]=RDLT_Planar1(Xw,xxn)
n=length(Xw);

index=combntns(1:n,2)';
W1=Xw(:,index(1,:));
W2=Xw(:,index(2,:));
x1=xxn(:,index(1,:));
x2=xxn(:,index(2,:));

x1=[x1;ones(1,length(x1))];
x2=[x2;ones(1,length(x2))];

nc=cross(x1,x2); 

%点的约束;
ui=xxn(1,:).';   vi=xxn(2,:).';
Xi=Xw(1,:).';   Yi=Xw(2,:).';
A=[Xi,zeros(n,1),-Xi.*ui,Yi,zeros(n,1),-Yi.*ui,ones(n,1),zeros(n,1),-ui;
   zeros(n,1),Xi,-Xi.*vi,zeros(n,1),Yi,-Yi.*vi,zeros(n,1),ones(n,1),-vi];

%点和点的约束;
nx=[nc(1,:),nc(1,:)].'; ny=[nc(2,:),nc(2,:)].'; nz=[nc(3,:),nc(3,:)].';
Pt=[W1,W2];
Xi=Pt(1,:).';   Yi=Pt(2,:).';
B=[Xi.*nx,Xi.*ny,Xi.*nz,Yi.*nx,Yi.*ny,Yi.*nz,nx,ny,nz];

M=[A;B];
M1=M(:,1:8);
M2=M(:,end);
s=-M1\M2;
H=reshape(s(1:6),3,2);
Rr= SVDdecomposition(H*Xw(1:2,:),Xw);

D=zeros(2*n,3);
E=zeros(2*n,1);
D(1:2:end,:)=[ones(n,1),zeros(n,1),-xxn(1,:).'];
D(2:2:end,:)=[zeros(n,1),ones(n,1),-xxn(2,:).'];
E(1:2:end,:)=(Rr(3,1:2)*Xw(1:2,:).*xxn(1,:)-Rr(1,1:2)*Xw(1:2,:)).';
E(2:2:end,:)=(Rr(3,1:2)*Xw(1:2,:).*xxn(2,:)-Rr(2,1:2)*Xw(1:2,:)).';
Rc=Rr;
tc=((D.'*D)\D.')*E;
end

function R=SVDdecomposition(Xc,Xw)
% This routine solves the exterior orientation problem for a point cloud
%  given in both camera and world coordinates. It is described in:
%Umeyama S. Least-Squares Estimation of Transformation Parameters Between 2 Point Patterns[J]. 
%IEEE Transactions on Pattern Analysis & Machine Intelligence, 1991, 13(4):376-380.

%% using SVD's method;
n=size(Xc,2);
% derive the centroid of the two point-clouds
p1=sum(Xc,2)/n;
%compute the matrix H = sum(F'*G^{T})
H1=Xc-p1*ones(1,n);
H2=Xw-p1*ones(1,n);
H=H1*H2';
%decompose this matrix (SVD) to obtain rotation
[U,~,V]=svd(H);
%calculate R and T;
R=U*V';  %or R=(V*U')';
if det(R)<0
    V_prime(:,1)=V(:,1);
    V_prime(:,2)=V(:,2);
    V_prime(:,3)=-V(:,3);
    R=U*V_prime';
end

end

function solution= RefineGaussNewton(solution,E,G,tc)
    %refine the solution by using Gauss Newton method
    maxItr=5;  %maximum allowed iterations;
    lambda=1e-8;    %damped factor;
    maxLambda=1e2;  %max lambda;
    minLambda=1e-8; %min lambda;
    tz=sqrt(tc(3));
    solution=solution./tz;
    a=solution(1); b=solution(2); c=solution(3); d=solution(4);
    w=[1,a^2,a*b,a*c,a*d,b^2,b*c,b*d,c^2,c*d,d^2].';
    obj_pre = w.'*G*w;     %objective function;
    itr=1;
    %iteration
    while itr<=maxItr
        a=solution(1); b=solution(2); c=solution(3); d=solution(4);
        w=[1,a^2,a*b,a*c,a*d,b^2,b*c,b*d,c^2,c*d,d^2].';
        %Jacobian matrixa
        Jac=[ 0,     0,     0,    0;
              2*a,   0,     0,    0;
              b,     a,     0,    0;
              c,     0,     a,    0;
              d,     0,     0,    a;
              0,   2*b,     0,    0;
              0,     c,     b,    0;
              0,     d,     0,    b;
              0,     0,   2*c,    0;
              0,     0,     d,    c;
              0,     0,     0,  2*d];
        Fk=E*w;
        Jk=E*Jac;
        solution_temp = solution;
        while lambda<maxLambda
            %increment
            dk=-((Jk.'*Jk)+lambda*eye(4))\(Jk.'*Fk);
            %update parameter
            solution = solution_temp + dk;
            a=solution(1); b=solution(2); c=solution(3); d=solution(4);
            w=[1,a^2,a*b,a*c,a*d,b^2,b*c,b*d,c^2,c*d,d^2].';
            %evaluate thee rror;
            obj_cur = w.'*G*w;
            %check convergence
            if obj_cur >= obj_pre
                lambda = 10*lambda;
                continue;
            else
                obj_pre = obj_cur; 
                lambda = 0.1*lambda;
                break;
            end 
        end
        if lambda >= maxLambda
            solution = solution_temp;
            break;
        end
        
        if lambda <= minLambda
            lambda = minLambda;
        end
        itr = itr + 1;
    end
    
end

function Q = matrix2quaternion(T)

    % This code follows the implementation suggested by Hartley and Zisserman    
    R = T(1:3, 1:3);   % Extract rotation part of T
    
    % Find rotation axis as the eigenvector having unit eigenvalue
    % Solve (R-I)v = 0;
    [v,d] = eig(R-eye(3));
    
    % The following code assumes the eigenvalues returned are not necessarily
    % sorted by size. This may be overcautious on my part.
    d = diag(abs(d));   % Extract eigenvalues
    [s, ind] = sort(d); % Find index of smallest one
    if d(ind(1)) > 0.001   % Hopefully it is close to 0
        warning('Rotation matrix is dubious');
    end
    
    axis = v(:,ind(1)); % Extract appropriate eigenvector
    
    if abs(norm(axis) - 1) > .0001     % Debug
        warning('non unit rotation axis');
    end
    
    % Now determine the rotation angle
    twocostheta = trace(R)-1;
    twosinthetav = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';
    twosintheta = axis'*twosinthetav;
    
    theta = atan2(twosintheta, twocostheta);
    
    Q = [cos(theta/2); axis*sin(theta/2)];
end