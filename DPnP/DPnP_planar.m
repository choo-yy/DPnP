function [R, t] = DPnP_planar(P, Q)
    n = size(Q, 2);
    %% find points with the maximum and minimum depths using homography 
    H = getH(P(1:2,:), Q(1:2,:));
    r3 = cross(H(:,1), H(:,2));
    r3 = r3/norm(r3);
    v = [r3(1); r3(2)];
    l = v'*Q(1:2,:);
    
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
        Rs(:,:,i) = R; ts(:,i) = t;
    end
    if n_solution == 1
        minindex = 1;
    else
        for i = 1:n_solution
            R = Rs(:,:,i); t = ts(:,i);
            Qc_ = R*P + t; Qc_ = Qc_(1:2, :)./Qc_(3,:);
            reproj_err(i) = norm(Qc_ - Q, 'fro')/n;
        end
        [~, minindex] = min(reproj_err);
    end
    R = Rs(:,:,minindex); t = ts(:,minindex);

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