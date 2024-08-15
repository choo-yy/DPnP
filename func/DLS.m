function [R,t] = DLS(XX,xx)
R = eye(3); t = zeros(3,1);
[R, t, cost] = robust_dls_pnp(XX,xx);
if size(R,3)>1
    R = R(:,:,1);
    t = t(:, 1);
elseif size(R,3)==0
    R = eye(3); t = zeros(3,1);
end

return