%******************************************************************************
% This programe is implemented via MATLAB 2018.                              *
% Author :  Ping Wang                                                        *
% Contact:  pingwangsky@gmail.com; pingwangsky@163.com                       *
% License:  Copyright (c) 2020 Ping Wang, All rights reserved.               *
% My site:  https://sites.google.com/view/ping-wang-homepage                 *
%*****************************************************************************/
function [Rc,tc]=RDLT2(Xw,xxn)
    % Input checks
    if (length(Xw) < 4)
        error('At least  4  points are supplied.');
    end
    
    if rank(Xw)>2
        % ordinar case;
        [Rc,tc]=RDLT_Ordinary(Xw,xxn);
    else
        % planar case;
        [Rc,tc]=RDLT_Planar(Xw,xxn);
    end
end


function [Rc,tc]=RDLT_Ordinary(Xw,xxn)
n=length(Xw);
% if(n<=15)
%     index=nchoosek(1:n,2)';
%     W1=Xw(:,index(1,:));
%     W2=Xw(:,index(2,:));
%     x1=xxn(:,index(1,:));
%     x2=xxn(:,index(2,:));
% end


    index1=randperm(n);
    index2=randperm(n);
    W1=[Xw(:,1:end),Xw(:,index1)];
    W2=[Xw(:,[2:end,1]),Xw(:,index2)];
    x1=[xxn(:,1:end),xxn(:,index1)];
    x2=[xxn(:,[2:end,1]),xxn(:,index2)];


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

function [Rc,tc]=RDLT_Planar(Xw,xxn)
n=length(Xw);
% if(n<=15)
%     index=combntns(1:n,2)';
%     W1=Xw(:,index(1,:));
%     W2=Xw(:,index(2,:));
%     x1=xxn(:,index(1,:));
%     x2=xxn(:,index(2,:));
% end


    index1=randperm(n);
    index2=randperm(n);
    W1=[Xw(:,1:end),Xw(:,index1)];
    W2=[Xw(:,[2:end,1]),Xw(:,index2)];
    x1=[xxn(:,1:end),xxn(:,index1)];
    x2=[xxn(:,[2:end,1]),xxn(:,index2)];


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