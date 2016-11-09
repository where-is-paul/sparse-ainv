function ildl_testscript(problem,pvtsize,droptol,droptype,shift,shiftA_flag,linsolve_tol)
%%% pvtsize = 2 for this study; pivot of larger size to be explored in future

% Test examples (zero shifts) 
% case9: flag_wmn = true, threshold = 0.75
% ildl_testscript('case9.mtx','2','1e-7','relative','0','false','1e-7');
% OPF_3754: flag_grook = true, threshold = 0.1
% ildl_testscript('OPF_3754.mtx','2','5e-4','absolute','0','false','1e-8');

flag_gbk = (0 > 0);
flag_grook = (1 > 0);
flag_wmn = (0 > 0);

thresholds = [0.1];

linsolve_tol = str2num(linsolve_tol);
droptol = str2num(droptol);
pvtsize = str2num(pvtsize);
fprintf(strcat('\nProblem\t',problem,'\n'));
A = mmread(problem);
shift = str2num(shift);
B = A - shift*speye(size(A));
if strcmpi(shiftA_flag,'true')
    A = B;
    fprintf('Coefficient matrix is A-s*I with s = %d. ILDL also constructed for A-s*I.\n\n',shift);
else
    fprintf('Coefficient matrix is A. ILDL constructed for A-s*I with s = %d.\n\n',shift);
end
    
fprintf('Performing LDL factorization of A.\n');
[L,D,~,~] = ldl(A,0.5);
fprintf('LDLT factorization done. Computing the inertia of A.\n');
[np,nn,nz] = inertia_blkdiag(D);
fprintf('The inertia of A is [%d %d %d]\n',np,nn,nz);
nnzl = nnz(L);
%%% nnzinvl = nnzl;  %nnz(invL);

fprintf('nnz(A) = %d, nnz(L) = %d\n\n',nnz(A),nnzl);
fprintf('Computing a favorable diagonally strong permutation.\n');
fprintf('Using hsl_mc64... \n');
n = length(A);
LB = tril(B);
%%% [perm_row,perm_col,info,scale_row,scale_col] = hsl_mc64(LB,5,1);
[~,perm_col,~,scale_row,~] = hsl_mc64(LB,5,1);
S_prpc = spdiags(exp(scale_row),0,n,n);
P_prpc = speye(size(A));
P_prpc = P_prpc(:,perm_col);

fprintf('Diagonally strong permutation computed.\n\n');
SB = P_prpc'*(S_prpc*B*S_prpc)*P_prpc;
SA = P_prpc'*(S_prpc*A*S_prpc)*P_prpc;
qq = symamd(SB);
P_amd = speye(size(A));    P_amd = P_amd(:,qq);

n = length(A);
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);
b = randn(n,1);     b = b/norm(b);      %b = P_prpc'*(S_prpc*b);

fprintf('Solving the original linear system by unpreconditioned SQMR...\n');
[~,flag,relres,steps] = sqmr(A,b,linsolve_tol,1);
if flag == 0
    fprintf('Unpreconditioned SQMR converged in %d steps.\n',steps);
else
    fprintf('Unpreconditioned SQMR did not converge within %d steps ...\n',n);
    fprintf('Minimal residual %d reached at step %d.\n',relres,steps);
end
fprintf('Solving the MC64 preprocessed linear system by unpreconditioned SQMR...\n');
[~,flag,relres,steps] = sqmr(SA,P_prpc'*(S_prpc*b),linsolve_tol,1);
if flag == 0
    fprintf('Unpreconditioned SQMR converged in %d steps.\n\n',steps);
else
    fprintf('Unpreconditioned SQMR did not converge within %d steps ...\n',n);
    fprintf('Minimal residual %d reached at step %d.\n\n',relres,steps);
end

lindata.A = A;      lindata.b = b;      lindata.SB = SB;        lindata.qq = qq;
lindata.P_amd = P_amd;    lindata.P_prpc = P_prpc;	lindata.S_prpc = S_prpc;    lindata.nnzl = nnzl;

if flag_gbk
    mincost = realmax;
    for outer = 1 : length(thresholds)
        threshold = thresholds(outer);
        fprintf('Outer iteration %d: testing GBK pivoting with treshold %.3f.\n',outer,threshold);
        optdata = matrixtest_update_ildl(lindata,pvtsize,'gbk',threshold,droptol,droptype,linsolve_tol);
        if optdata.cost < mincost
            mincost = optdata.cost;
            minoptdata = optdata;
        end
        fprintf('\n');
    end
    fprintf('Optimal preconditioner (GBK) droptol = %.2d, threshold = %.2d\n\tdensity factors = [%.3d %.3d], iter = %d, cost = %.2d\n\n\n',...
        minoptdata.droptol,minoptdata.threshold,minoptdata.density(1),minoptdata.density(2),minoptdata.iter,minoptdata.cost);
end

if flag_grook
    mincost = realmax;
    for outer = 1 : length(thresholds)
        threshold = thresholds(outer);
        fprintf('Outer iteration %d: testing GROOK pivoting with treshold %.3f.\n',outer,threshold);
        optdata = matrixtest_update_ildl(lindata,pvtsize,'grook',threshold,droptol,droptype,linsolve_tol);
        if optdata.cost < mincost
            mincost = optdata.cost;
            minoptdata = optdata;
        end
        fprintf('\n');
    end
    fprintf('Optimal preconditioner (GROOK) droptol = %.2d, threshold = %.2d\n\tdensity factors = [%.3d %.3d], iter = %d, cost = %.2d\n\n\n',...
        minoptdata.droptol,minoptdata.threshold,minoptdata.density(1),minoptdata.density(2),minoptdata.iter,minoptdata.cost);
end

if flag_wmn
    mincost = realmax;
    for outer = 1 : length(thresholds)
        threshold = thresholds(outer);
        fprintf('Outer iteration %d: testing WMN pivoting with treshold %.3f.\n',outer,threshold);
        optdata = matrixtest_update_ildl(lindata,pvtsize,'wmn',threshold,droptol,droptype,linsolve_tol);
        if optdata.cost < mincost
            mincost = optdata.cost;
            minoptdata = optdata;
        end
        fprintf('\n');
    end
    fprintf('Optimal preconditioner (WMN) droptol = %.2d, threshold = %.2d\n\tdensity factors = [%.3d %.3d], iter = %d, cost = %.2d\n\n\n',...
        minoptdata.droptol,minoptdata.threshold,minoptdata.density(1),minoptdata.density(2),minoptdata.iter,minoptdata.cost);
end
end