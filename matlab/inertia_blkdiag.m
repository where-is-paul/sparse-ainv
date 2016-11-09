function [np,nn,nz] = inertia_blkdiag(D)
n = length(D);
normD = norm(D,'fro');
threshold = eps/2^8*normD;
np = 0;
nn = 0;
nz = 0;
kk = 1;
while kk <= n
    if kk < n && D(kk+1,kk) ~= 0
        ev = eig(full(D(kk:kk+1,kk:kk+1)));
        kk = kk + 2;
    else
        ev = D(kk,kk);
        kk = kk + 1;
    end
    nz = nz + sum(abs(ev) < threshold);
    np = np + sum(ev >= threshold);
    nn = nn + sum(ev <= -threshold);
end