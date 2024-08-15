function [R, t] = DPnP_GN_planar(P, Q)
    n = size(Q, 2);
    %% find points with the maximum and minimum depths using homography 
    H = getH(P(1:2,:), Q);
    r3 = cross(H(:,1), H(:,2));
    r3 = r3/norm(r3);
    v = [r3(1); r3(2)];
    l = v'*Q;
    
    [~, minindex] = max(l);
    P(:, [1 minindex]) = P(:, [minindex, 1]);
    Q(:, [1 minindex]) = Q(:, [minindex, 1]);
    l(:, [1 minindex]) = l(:, [minindex, 1]);
    [~, minindex] = min(l);
    P(:, [2 minindex]) = P(:, [minindex, 2]);
    Q(:, [2 minindex]) = Q(:, [minindex, 2]);

    %% DPnP solution
    P1 = P(:,1); P2 = P(:,2); P3e = P(:,3:end);
    a = norm(P2 - P1);
    b = (P3e - P1)'*(P2 - P1)/a;
    c = sqrt(sum((P3e - P1).^2)' - b.^2);
    m = [Q; ones(1, n)]./sqrt(sum(Q.^2) + 1);
    m1 = m(:,1); m2 = m(:,2); m3e = m(:,3:end);
    f1 = b/a; f2 = -m3e'*m2; f4 = (1-2*f1)*(m1'*m2); f5 = m3e'*m1; f6 = f1-1;
    g1 = (b.^2+c.^2)/a^2; g4 = -2*g1*(m1'*m2); g5 = 2*f5; g6 = g1-1;
    
    f22 = f2.^2; f52 = f5.^2; f25x2 = 2*f2.*f5;
    h1 = f22.*g1 - f1.^2;
    h2 = f22.*g4 + f25x2.*(g1 -f1) - 2*f1.*f4;
    h3 = f52.*(g1 - 2*f1) + f25x2.*(g4 - f4) - 2*f1.*f6 + f22.*g6 - f4.^2;
    h4 = f52.*(g4 - 2*f4) + f25x2.*(g6 - f6) - 2*f4.*f6;
    h5 = f52.*(g6 - 2*f6) - f6.^2;
    
    k8 = h1'*h1;
    k7 = 2*h1'*h2;
    k6 = h2'*h2 + 2*h1'*h3;
    k5 = 2*(h1'*h4 + h2'*h3);
    k4 = h3'*h3 + 2*(h1'*h5 + h2'*h4);
    k3 = 2*(h2'*h5 + h3'*h4);
    k2 = h4'*h4 + 2*h3'*h5;
    k1 = 2*h4'*h5;
    % k0 = h5'*h5;
    
    xs = roots([8*k8 7*k7 6*k6 5*k5 4*k4 3*k3 2*k2 k1]);
    maxreal = max(abs(real(xs)));
    xs(abs(imag(xs))/maxreal > 0.001) = [];
    xs = real(xs);
    xs(xs<0.001) = [];
    % roots of local minimums
    F6 = polyval([56*k8 42*k7 30*k6 20*k5 12*k4 6*k3 2*k2], xs);
    xs(F6 <= -1e-5) = [];
    % solve absolute orientation problem and choose R t
    n_solution = length(xs);
    Rs = zeros(3,3,n_solution); ts = zeros(3, n_solution); reproj_err = zeros(n_solution, 1);
    Pb = mean(P, 2);
    sPt2 = sum((P - Pb).^2, "all");
    for i = 1:n_solution
        dd1 = [1; xs(i); -(xs(i)^2*f1 + xs(i)*f4 + f6)./(xs(i)*f2 + f5)];
        x = xs(i);
        % refine di/d1
        sigma1 = (f5 + f2*x); sigma2 = (g1*x^2 + g4*x + g6);
        K1 = 2; K2 = -3*g5; K3 = sigma1.^2 - sigma2*2 + g5.^2;
        K4 = sigma1.*(f1*x^2 + f4*x + f6) + g5.*sigma2;
        dd3e = dd1(3:end);
        for iter = 1:10
            y0 = K1*dd3e.^3 + K2.*dd3e.^2 + K3.*dd3e + K4;
            k = 3*K1*dd3e.^2 + 2*K2.*dd3e + K3;
            dd3e = dd3e - y0./k;
        end
        dd1 = [dd1(1:2);dd3e];
        
        Qs = m.*repmat(dd1', 3, 1);
        Qsb = mean(Qs, 2);
        Qst = Qs - Qsb;
        d1 = sqrt(sPt2/sum(Qst.^2, "all"));
        d = d1*dd1;
        Qs = m.*repmat(d',3,1);
        [R, t] = absolute_orientation(P,Qs);
        Qc_ = R*P + t; Qc_ = Qc_(1:2, :)./Qc_(3,:);
        reproj_err(i) = norm(Qc_ - Q, 'fro')/n;
        Rs(:,:,i) = R; ts(:,i) = t;
    end

    [~, minindex] = min(reproj_err);
    R = Rs(:,:,minindex); t = ts(:,minindex);
    %% Gauss Newton optimization
    x = P(1,:)'; y = P(2,:)'; z = P(3,:)'; u = Q(1,:)'; v = Q(2,:)';
    Pub = mean(u.*P')';  Pvb = mean(v.*P')';
    ub  = mean(u);    ut  = u - ub; 
    vb  = mean(v);    vt  = v - vb;
    xt   = x - Pb(1);           yt   = y - Pb(2);           zt2 = 2*(z - Pb(3));
    uxt2 = 2*(x.*u - Pub(1));   uyt2 = 2*(y.*u - Pub(2));   uzt = z.*u - Pub(3);
    vxt2 = 2*(x.*v - Pvb(1));   vyt2 = 2*(y.*v - Pvb(2));   vzt = z.*v - Pvb(3);
    E = [-ut, xt-uzt,     -uyt2, zt2+uxt2, -2*yt,  xt+uzt, 2*yt, zt2-uxt2, -xt+uzt,    -uyt2, -xt-uzt;
         -vt, yt-vzt, -zt2-vyt2,     vxt2,  2*xt, -yt+vzt, 2*xt,    -vxt2,  yt+vzt, zt2-vyt2, -yt-vzt];

    q = matrix2quaternion(R)/sqrt(t(3));
    a = q(1); b = q(2); c = q(3); d = q(4);
    F = E*[1, a^2, a*b, a*c, a*d, b^2, b*c, b*d, c^2, c*d, d^2]';
    J = [E(:,2:5)*[2*a, b, c, d]', E(:,[3, 6, 7, 8])*[a, 2*b, c, d]', E(:,[4, 7, 9, 10])*[a, b, 2*c, d]', E(:,[5, 8, 10, 11])*[a, b, c, 2*d]'];
    q = q - (J'*J)\(J'*F);

    a = q(1); b = q(2); c = q(3); d = q(4);
    Rw=[a^2+b^2-c^2-d^2,     2*b*c-2*a*d,     2*b*d+2*a*c
            2*b*c+2*a*d, a^2-b^2+c^2-d^2,     2*c*d-2*a*b
            2*b*d-2*a*c,     2*c*d+2*a*b, a^2-b^2-c^2+d^2];
    factor = 1/(q'*q);
    R2 = Rw*factor;
    t2 = [-Rw(1,:)*Pb+Rw(3,:)*Pub+ub; -Rw(2,:)*Pb+Rw(3,:)*Pvb+vb; 1]*factor;

    Qc_ = R2*P + t2; Qc_ = Qc_(1:2, :)./Qc_(3,:);
    reproj_err2 = norm(Qc_ - Q, 'fro')/size(P,2);
    if reproj_err2 < reproj_err(minindex)
        R = R2; t = t2;
    end

end


function [R,t] = absolute_orientation(A,B)
% Registers two sets of 3DoF data
% Assumes A and B are d,n sets of data
% where d is the dimension of the system 
% typically d = 2,3
% and n is the number of points
% typically n>3
%
% Mili Shah
% July 2014
 
[~,n]=size(A);
 
%Mean Center Data
Ac = mean(A,2);
Bc = mean(B,2);
A = A-repmat(Ac,1,n);
B = B-repmat(Bc,1,n);
 
%Calculate Optimal Rotation
M = A*B';
N = [M(1,1) + M(2,2) + M(3,3) M(2,3)-M(3,2) M(3,1)-M(1,3) M(1,2)-M(2,1);...
    M(2,3)-M(3,2) M(1,1)-M(2,2)-M(3,3) M(1,2)+M(2,1) M(3,1)+M(1,3);...
    M(3,1)-M(1,3) M(1,2)+M(2,1) -M(1,1)+M(2,2)-M(3,3) M(2,3)+M(3,2);...
    M(1,2)-M(2,1) M(3,1)+M(1,3) M(2,3)+M(3,2) -M(1,1)-M(2,2)+M(3,3)];
[u,v]=eig(N);
[~,ind]=max(diag(v));
r = u(:,ind);
r = [r(2:4);r(1)];
R = (r(4)^2-r(1:3)'*r(1:3))*eye(3) + 2*r(1:3)*r(1:3)'+2*r(4)*[0 -r(3) r(2);r(3) 0 -r(1);-r(2) r(1) 0];
%Calculate Optimal Translation
t = Bc - R*Ac;
end

function Q = matrix2quaternion(R)
    % This code follows the implementation suggested by Hartley and Zisserman    
    % R = T(1:3, 1:3);   % Extract rotation part of T
    
    % Find rotation axis as the eigenvector having unit eigenvalue
    % Solve (R-I)v = 0;
    [v,d] = eig(R-eye(3));
    
    % The following code assumes the eigenvalues returned are not necessarily
    % sorted by size. This may be overcautious on my part.
    d = diag(abs(d));   % Extract eigenvalues
    [~, ind] = sort(d); % Find index of smallest one
    
    axis = v(:,ind(1)); % Extract appropriate eigenvector
    
    % Now determine the rotation angle
    twocostheta = trace(R)-1;
    twosinthetav = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';
    twosintheta = axis'*twosinthetav;
    
    theta = atan2(twosintheta, twocostheta);
    
    Q = [cos(theta/2); axis*sin(theta/2)];
end

function H = getH(U,Q)
%n is the number of points:
n = size(U,2);
modelDims = size(U,1);

if modelDims==2
   %we zero-center the model points:
    Pbar = [mean(U,2);0];
    U(1,:) = U(1,:)-Pbar(1);
    U(2,:) = U(2,:)-Pbar(2);
else
    %we rotate the model points onto the plane z=0 and zero center them:
    Pbar = mean(U,2);   
    MCenter = eye(4);
    MCenter(1:3,end) = -Pbar;    
    U_ = MCenter(1:3,:)*[U;ones(1,size(U,2))];  
    %compute the rotation that rotates the model points onto the plane z=0: 
    [modelRotation,~,~] = svd(U_*U_'); modelRotation=modelRotation';
    if det(modelRotation)<0
        modelRotation(3,:)  = -modelRotation(3,:); %this ensures modelRotation is a member of SO3
    end
    
    modelRotation(4,4) = 1;
    Mcorrective = modelRotation*MCenter;
    U = Mcorrective(1:2,:)*[U;ones(1,size(U,2))];
    
end

%compute the model-to-image homography:
H = homographyHarker([U;ones(1,n)],[Q;ones(1,n)]);

%Compute the Jacobian J of the homography at (0,0):
H = H./H(3,3);

end

function H = homographyHarker( DataA, DataB)

    % Normalize the input data:
    [DataA,TA,TAi] = normalizeData( DataA ) ;
    [DataB,TB,TBi] = normalizeData( DataB ) ;
    % Construct the orthogonalized design matrix :
    C1 = -DataB(1,:) .* DataA(1,:) ;
    C2 = -DataB(1,:) .* DataA(2,:) ;
    %
    C3 = -DataB(2,:) .* DataA(1,:) ;
    C4 = -DataB(2,:) .* DataA(2,:) ;
    %
    mC1 = mean( C1 ) ;
    mC2 = mean( C2 ) ;
    mC3 = mean( C3 ) ;
    mC4 = mean( C4 ) ;
    %
    Mx = [ C1' - mC1, C2' - mC2, -DataB(1,:)' ] ;
    My = [ C3' - mC3, C4' - mC4, -DataB(2,:)' ] ;
    %
    Pp = pinv(DataA(1:2,:)') ;
    %
    Bx = Pp * Mx ;
    By = Pp * My ;
    %
    D = [ Mx - DataA(1:2,:)'*Bx ;...
          My - DataA(1:2,:)'*By ] ;
    %
    % Find v_min and backsubstitute :
    %
    [~,~,V] = svd( D, 0 ) ;
    %
    h789 = V(:,end) ;
    h12 = -Bx * h789 ;
    h45 = -By * h789 ;
    h3 = -[mC1, mC2] * h789(1:2) ;
    h6 = -[mC3, mC4] * h789(1:2) ;
    %
    % Reshape vector h to matrix H, and transform :
    %
    H = reshape( [h12; h3; h45; h6; h789], 3, 3)' ;
    %
    H = TB * H * TAi ;
end

function [DataN, T, Ti] = normalizeData( Data )

s = Data(1,:) ;
t = Data(2,:) ;
u = Data(3,:) ;
%
x = s./u ;
y = t./u ;
%
xm = mean( x ) ;
ym = mean( y ) ;
%
xh = x - xm ;
yh = y - ym ;
%
n = length( xh ) ;
%
kappa = sum( xh.^2 + yh.^2 ) ;
%
beta = sqrt( 2*n / kappa ) ;
%
xn = beta * xh ;
yn = beta * yh ;
%
DataN = [ xn; yn; ones(size(xn)) ] ;
%
T = [ 1/beta,   0   , xm ;...
        0   , 1/beta, ym ;...
        0   ,   0   ,  1 ] ;
%
Ti = [ beta ,  0   , -beta * xm ; ...
        0   , beta , -beta * ym ; ...
        0   ,  0   ,       1    ] ;
end