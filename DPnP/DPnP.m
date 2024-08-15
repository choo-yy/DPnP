function [R, t] = DPnP(P, Q)
    %% estimate R t
    n = size(P,2);
    Pw = [P;ones(1,n)];  PPT = Pw*Pw';
    Pwuv = Pw*Q';        IJp = PPT\Pwuv;
    Ip = IJp(1:3,1);   x0 = IJp(4,1);
    Jp = IJp(1:3,2);   y0 = IJp(4,2);
    
    I_anti = [0, -Ip(3), Ip(2)
              Ip(3), 0, -Ip(1)
             -Ip(2), Ip(1), 0];
    J_anti = [0, -Jp(3), Jp(2)
              Jp(3), 0, -Jp(1)
             -Jp(2), Jp(1), 0];
    tz = 1/2*(sqrt(1+x0^2)/norm(Ip) + sqrt(1+y0^2)/norm(Jp));
    k = (eye(3) - tz*y0*I_anti + tz*x0*J_anti)\cross(Ip,Jp)*tz^2;  k = k/norm(k);
    i = tz*Ip + x0*k; i = i/norm(i);
    j = tz*Jp + y0*k; j = j/norm(j);
    R = [i,j,k]';
    t = [x0*tz; y0*tz; tz];
    Pc = R*P + t;
    dc = vecnorm(Pc);

    [~, mindex] = max(dc);
    P(:, [1, mindex]) = P(:, [mindex, 1]);
    Q(:, [1, mindex]) = Q(:, [mindex, 1]);
    dc(:, [1, mindex]) = dc(:, [mindex, 1]);
    [~, mindex] = min(dc);
    P(:, [2, mindex]) = P(:, [mindex, 2]);
    Q(:, [2, mindex]) = Q(:, [mindex, 2]);


    %% DPnP solution
    P1 = P(:,1); P2 = P(:,2); P3e = P(:,3:end);
    a = norm(P2 - P1);
    b = (P3e - P1)'*(P2 - P1)/a;
    c = sqrt(sum((P3e - P1).^2)' - b.^2);
    m = [Q; ones(1, n)]./sqrt(sum(Q.^2) + 1);
    m1 = m(:,1); m2 = m(:,2); m3e = m(:,3:end);
    f1 = b/a; f2 = -m3e'*m2; f4 = (1-2*f1)*(m1'*m2); f5 = m3e'*m1; f6 = f1-1;
    g1 = (b.^2+c.^2)/a^2; g4 = -2*g1*(m1'*m2); g6 = g1-1;
    
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
    xs(F6 <= 0) = [];
    % solve absolute orientation problem and choose R t
    n_solution = length(xs);
    Rs = zeros(3,3,n_solution); ts = zeros(3, n_solution); reproj_err = zeros(n_solution, 1);
    Pb = PPT(1:3,4)/n;
    sPt2 = sum((P - Pb).^2, "all");
    for i = 1:n_solution
        dd1 = [1; xs(i); -(xs(i)^2*f1 + xs(i)*f4 + f6)./(xs(i)*f2 + f5)];
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