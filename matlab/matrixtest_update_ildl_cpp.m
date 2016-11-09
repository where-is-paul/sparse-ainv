function matrixtest_update_ildl_cpp(name, pvttype, pptol, droptol, sigma, shiftA_flag,linsolve_tol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
maxiter = 1;

%   linsolve_tol = 5e-8;
%   using linesolve_tol = 5e-8 generally, but sometimes may need other values
%   so used as a parameter in the funciton call now



%   Printing out name of command before actual program output for record-keeping

fprintf('\n--------------------------------------------\n');
fprintf('Using command (maxiter = %d):\n', maxiter);
fprintf('matrixtest_update_ildl_cpp(''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'') \n\n', name, pvttype, pptol, droptol, sigma, shiftA_flag, linsolve_tol)

%



linsolve_tol = str2num(linsolve_tol);

warning off;
tic;
droptol = str2num(droptol);
pptol = str2num(pptol);
fprintf(strcat('\n\n\nProblem\t',name,'\n'));
A = mmread(name);
if nargin > 3
    sigma = str2num(sigma);
    B = A-sigma*speye(size(A));
    fprintf('Using a shifted matrix B to construct the ILDL preconditioner.\n');
    if strcmpi(shiftA_flag,'true')
        A = B;
        fprintf('Matrix A itself is also shifted.\n');
    else
        fprintf('Matrix A itself is NOT shifted.\n');
    end
else
    B = A;
end

%fprintf('Unpreprocessed condest(B) = %d.\n',condest(B));

%     fprintf('Performing LDLT factorization of A.\n');
%     [L,D,~,~] = ldl(A,0.5);
%     fprintf('LDLT factorization done. Computing the inertia of A.\n');
%     [np,nn,nz] = inertia_blkdiag(D);
%     fprintf('The inertia of A is [%d %d %d]\n',np,nn,nz);
% %     fprintf('Computing inv(L).\t');
% %     invL = inv(L);
% %     fprintf('inv(L) computed.\n');
%     fprintf('nnz(A) = %d, nnz(L) = %d\n',nnz(A),nnz(L));
L = A;
fprintf('Computing a favorable diagonally strong permutation.\n');

%     [P_prpc,S_prpc] = symscaling(B);

fprintf('Using hsl_mc64... \n');
n = length(A);
LB = tril(B);
[perm_row,perm_col,info,scale_row,scale_col] = hsl_mc64(LB,5,1);
S_prpc = spdiags(exp(scale_row),0,n,n); %(where n = length(A) is the size of the matrix)
P_prpc = speye(size(A));
P_prpc = P_prpc(:,perm_col);

% -------------------------

fprintf('Diagonally strong permutation computed.\n');
SB = P_prpc'*(S_prpc*B*S_prpc)*P_prpc;	SB = (SB+SB')/2;
%fprintf('Mc64 preprocessed condest(B) = %d.\n',condest(SB));

qq = symamd(SB);
P_amd = speye(size(A));    P_amd = P_amd(:,qq);

n = length(A);
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);
b = randn(n,1);     b = b/norm(b);      %b = P_prpc'*(S_prpc*b);

%     fprintf('Solving the original linear system by unpreconditioned SQMR...\n');
%     [~,flag,relres,steps] = sqmr(A,b,linsolve_tol,1);
%     if flag == 0
%         fprintf('Unpreconditioned SQMR converged in %d steps.\n',steps);
%     else
%         fprintf('Unpreconditioned SQMR did not converge within %d steps ...\n',1);
%         fprintf('Minimal residual %d reached at step %d.\n',relres,steps);
%     end
%     fprintf('Solving the scaled/permuted linear system by unpreconditioned SQMR...\n');
%     [~,flag,relres,steps] = sqmr(SA,P_prpc'*(S_prpc*b),linsolve_tol,1);
%     if flag == 0
%         fprintf('Unpreconditioned SQMR converged in %d steps.\n',steps);
%     else
%         fprintf('Unpreconditioned SQMR did not converge within %d steps ...\n',1);
%         fprintf('Minimal residual %d reached at step %d.\n',relres,steps);
%     end
%     percent1 = blk_offon_diagratio(A,(1:n)');
%     fprintf('Diagonal dominance statistics of the original A:\n\n');
%     disp(100*percent1');
%     percent2 = blk_offon_diagratio(SA,(1:n)');
%     fprintf('Diagonal dominance statistics of the scaled/permuted A:\n\n');
%     disp(100*percent2');
%     nnzL = abs(nonzeros(L));
%     percent3 = zeros(5,1);
%     percent3(1) = nnz(nnzL <= 1e-4);
%     percent3(2) = nnz((1e-4 < nnzL).*(nnzL <= 1e-2));
%     percent3(3) = nnz((1e-2 < nnzL).*(nnzL <= 1e-1));
%     percent3(4) = nnz((1e-1 < nnzL).*(nnzL <= 1));
%     percent3(5) = nnz(1 < nnzL);
%     percent3 = percent3/length(nnzL);
%     fprintf('Magnitude statistics of elements in inv(L):\n\n');
%     disp(100*percent3');

dropworks = false;
failure = 0;
max_sqmrsteps = 1000;
minnnzL = realmax;
fprintf('Constructing iLDLt with a sequence of drop tolerances ...\n');
%     if pvtsize <= 2
%         fprintf('Classical BK-AINV method.\n');
%     else
%         fprintf('New large pivoting AINV method. Max pivot size = %d. Pivot selection opt %d.\n',pvtsize,pvtopt);
%     end
%pvttype = 'rook';
iter = 1;
while iter <= maxiter && 1e-10 <= droptol && droptol < 1 % needed for certain problems such as symmindft_easy.m
    fprintf('\nIteration %d, ILDL uses %s pivoting, droptol = %d, pp_tol = %d\n',iter,pvttype,droptol,pptol);
    %[LL,DD,pp] = ildl_largepvt_063015B(SB(qq,qq),pvtsize,pvtopt,droptol);
    %[LL,DD,pp] = sam_ILDL_cleaned_071715(SB(qq,qq),droptol,3);
    %[LL, DD, pp, S_ildl, ~] = ildl(SB(qq,qq), 1e4, droptol, pptol, 'none');
    t_start1 = tic;
    [LL,DD,pp,S_ildl,~] = ildl(SB(qq,qq), 1e10, droptol, 'amd','n',pvttype,pptol);
    t_elapsed1 = toc(t_start1);
    %         fprintf('The density of L over that of A and MATLAB''s L is [%d %d]. Condest(D) = %d.\n',...
    %             nnz(LL)/nnz(A),nnz(LL)/nnz(L),condest(D));
    fprintf('The density of L over that of A and MATLAB''s L is [%d %d].\n',...
        nnz(LL)/nnz(A),nnz(LL)/nnz(L));
    P_ildl = speye(size(A));   P_ildl = P_ildl(:,pp);
    mfun = @(v) S_prpc*(P_prpc*(P_amd*(S_ildl*(P_ildl*(LL'\(DD\(LL\(P_ildl'*(S_ildl*(P_amd'*(P_prpc'*(S_prpc*v))))))))))));
    fprintf('Solving the original linear system by PSQMR ...\n');
    t_start2 = tic;
    [~,flag,relres,steps] = sqmr(A,b,linsolve_tol,max_sqmrsteps,mfun);
    t_elapsed2 = toc(t_start2);
    if flag ~= 0
        fprintf('PSQMR did not converge within %d steps ...\n',max_sqmrsteps);
        fprintf('Minimal residual %d reached at step %d.\n',relres,steps);
        if ~dropworks
            droptol = droptol/5;
        else
            if mod(iter,2) == 1,	droptol = droptol*3/2;
            else			droptol = droptol*4/3;
            end
            failure = failure+1;
            if failure >= 2
                break;
            end
        end
    else
        fprintf('PSQMR converged at step %d.\n',steps);
        if nnz(LL) < minnnzL
            minnnzL = nnz(LL);
            fprintf('This is the sparsest preconditioner found so far ...\n');
            %save(strcat(name,'.precond.mat'),'P_prpc','S_prpc','LL','DD','pp','qq','droptol');
        end
        dropworks = true;
        if mod(iter,2) == 1,	droptol = droptol*3/2;
        else			droptol = droptol*4/3;
        end
        iter = iter + 1;
    end
    fprintf('Time for ildl factorization %.2f secs.\n',t_elapsed1);
    fprintf('Time for PSQMR linear solve %.2f secs.\n',t_elapsed2);
end
toc;
warning on;
end
