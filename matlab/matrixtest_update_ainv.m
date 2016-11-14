function optdata = matrixtest_update_ainv(lindata,pvtsize,pvtopt,threshold,...
    droptol,droptol_type,linsolve_tol)

warning off;
A = lindata.A;  b = lindata.b;  SB = lindata.SB;    qq = lindata.qq;
P_amd = lindata.P_amd;    P_prpc = lindata.P_prpc;  S_prpc = lindata.S_prpc; 
nnzl = lindata.nnzl;

n = length(A);
tic;
maxiter = 1;

dropworks = false;
failure = 0;
max_sqmrsteps = 1000;
minnnzM = realmax;
mincost = realmax;
fprintf('Constructing AINV with a sequence of drop tolerances ...\n');
if pvtsize <= 2
    fprintf('Basic sparsity-driven partial pivoting method.\n');
else
    fprintf('New large pivoting AINV method. Max pivot size = %d. Pivot selection optdata %d.\n',pvtsize,pvtopt);
end
iter = 1;
while iter <= maxiter && droptol <= 0.999 && droptol >= 1e-10
    fprintf('\nIteration %d, AINV droptol = %d (%s)\n',iter,droptol,droptol_type);
    t_start1 = tic;
    if pvtsize <= 2
        if strcmpi(pvtopt,'gbk')
            [MM,DD,pp] = spainv_sym_all(SB(qq,qq),'dummy',threshold,droptol,droptol_type,'gbk');
            %[MM,DD,pp] = spainv_sym_gbk(SB(qq,qq),threshold,droptol,droptol_type);
        elseif strcmpi(pvtopt,'grook')
            [MM,DD,pp] = spainv_sym_all(SB(qq,qq),'dummy',threshold,droptol,droptol_type,'grook');
            %[MM,DD,pp] = spainv_sym_grook(SB(qq,qq),threshold,droptol,droptol_type);
        elseif strcmpi(pvtopt,'wmn')
            [MM,DD,pp] = spainv_sym_all(SB(qq,qq),2,threshold,droptol,droptol_type,'wmn');
            %[MM,DD,pp] = spainv_sym_mwn(SB(qq,qq),2,threshold,droptol,droptol_type);
        else
            error('Supported pivoting strategies include gbk, grook and wmn only.\n');
        end
    else
        %[MM,DD,pp,~] = spapxinv_sym_uptodate_063015C(SB(qq,qq),droptol,pvtsize,pvtopt);
        error('pivot size must be 2 in this study');
    end
    t_elapsed1 = toc(t_start1);
    fprintf('nnz(M) over nnz of A, L, and dense trig is [%d %d %d]\n',...
        nnz(MM)/nnz(A),nnz(MM)/nnzl,nnz(MM)/(n*(n+1)/2));
    PP_ainv = speye(size(A));   PP_ainv = PP_ainv(:,pp);
    mfun = @(v) S_prpc*(P_prpc*(P_amd*(PP_ainv*(MM*(DD\(MM'*(PP_ainv'*(P_amd'*(P_prpc'*(S_prpc*v))))))))));
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
        if nnz(MM) < minnnzM
            minnnzM = nnz(MM);
            fprintf('This is the sparsest preconditioner found so far ...\n');
        end
        cost = steps*((2*nnz(MM)+nnz(DD))/nnz(A)+1);
        fprintf('Normalized arithmetic cost is %d.\n',cost);
        if cost < mincost
            mincost = cost;
            optdata.cost = cost;
            optdata.iter = steps;
            optdata.droptol = droptol;
            optdata.threshold = threshold;
            optdata.density = [nnz(MM)/nnz(A) nnz(MM)/nnzl];
        end
        dropworks = true;
        if mod(iter,2) == 1,	droptol = droptol*3/2;
        else			droptol = droptol*4/3;
        end
        iter = iter + 1;
    end
    fprintf('Timing: AINV factorization %.2f secs, SQMR solve %.2f secs.\n',t_elapsed1,t_elapsed2);
end
toc;
warning on;
end
