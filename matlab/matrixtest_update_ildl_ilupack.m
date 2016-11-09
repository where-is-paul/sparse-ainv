function matrixtest_update_ildl_ilupack(name, droptol, sigma, shiftA_flag,linsolve_tol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    maxiter = 1;

%   linsolve_tol = 5e-8;
%   using linesolve_tol = 5e-8 generally, but sometimes may need other values
%   so used as a parameter in the funciton call now



%   Printing out name of command before actual program output for record-keeping

fprintf('\n--------------------------------------------\n');
fprintf('Using command (maxiter = %d):\n', maxiter);
fprintf('matrixtest_update_ildl_ilupack(''%s'', ''%s'', ''%s'', ''%s'', ''%s'') \n\n',name, droptol, sigma, shiftA_flag, linsolve_tol)

%



    linsolve_tol = str2num(linsolve_tol);

    warning off;
    tic;
    droptol = str2num(droptol);
    %pptol = str2num(pptol);
    fprintf(strcat('\n\n\nProblem\t',name,'\n'));
    A = mmread(name);
    nmrA = normest(A,2);    %sqrt(norm(A,1)*norm(A,inf));
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
%     fprintf('Performing LDLT factorization of A.\n');
%     [L,D,~,~] = ldl(A,0.5);
%     fprintf('LDLT factorization done. Computing the inertia of A.\n');
%     [np,nn,nz] = inertia_blkdiag(D);
%     fprintf('The inertia of A is [%d %d %d]\n',np,nn,nz);
%     fprintf('Computing inv(L).\t');
%     invL = inv(L);
%     fprintf('inv(L) computed.\n');
    L = A;
    fprintf('nnz(A) = %d, nnz(L) = %d\n',nnz(A),nnz(L));
    
    
%     fprintf('Computing a favorable diagonally strong permutation.\n');
% %     [P_prpc,S_prpc] = symscaling(B);\
%         fprintf('Using hsl_mc64... \n');
%         n = length(A);
%         LB = tril(B);
%         [perm_row,perm_col,info,scale_row,scale_col] = hsl_mc64(LB,5,1);
%         S_prpc = spdiags(exp(scale_row),0,n,n); %(where n = length(A) is the size of the matrix)
%         P_prpc = speye(size(A));
%         P_prpc = P_prpc(:,perm_col);
%         % -------------------------
%     fprintf('Diagonally strong permutation computed.\n');
%     SB = P_prpc'*(S_prpc*B*S_prpc)*P_prpc;
%     %SA = P_prpc'*(S_prpc*A*S_prpc)*P_prpc;
%     qq = symamd(SB);
%     P_amd = speye(size(A));    P_amd = P_amd(:,qq);
    
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
    minnnzL = realmax;
    fprintf('Constructing iLDLt with a sequence of drop tolerances ...\n');
%     if pvtsize <= 2
%         fprintf('Classical BK-AINV method.\n');
%     else
%         fprintf('New large pivoting AINV method. Max pivot size = %d. Pivot selection opt %d.\n',pvtsize,pvtopt);
%     end
    iter = 1;
    while iter <= maxiter && 1e-10 <= droptol && droptol < 1 % needed for certain problems such as symmindft_easy.m
        fprintf('\nIteration %d, ILUPACK iLDLt droptol = %d\n',iter,droptol);
        clear options;
        options.droptol = droptol;  options.droptolS = droptol/10;
        options.nrestart = 1;       options.maxit = 1000;
        options.restol = linsolve_tol*20;
        options.isdefinite = 0;
	t_start1 = tic;
        [PREC,options] = AMGfactor(B,options);
	t_elapsed1 = toc(t_start1);
	fprintf('Time for ildl factorization %.2f secs.\n',t_elapsed1);
        %precnnz = AMGnnz_sym(PREC);
%         fprintf('The density of L over that of A and MATLAB''s L is [%d %d].\n',...
%             precnnz/nnz(A),precnnz/nnz(L));
%         fprintf('Solving the original linear system by PSQMR ...\n');
        [x,~] = AMGsolver(A,PREC,options,b);
        relres = norm(b-A*x)/norm(b);   %norm(b-A*x)/(nmrA * norm(x,inf) + norm(b,inf));
        options.restol = options.restol*(linsolve_tol/relres)*4;
        success = false;
        maxdd = 16;
        for dd = 1 : maxdd
	    t_start2 = tic;
            [x,options_out] = AMGsolver(A,PREC,options,b);
	    t_elapsed2 = toc(t_start2);
            relres = norm(b-A*x)/norm(b);
            if relres <= linsolve_tol && options_out.niter <= options.maxit
		fprintf('Success achieved in iteration %d.\n',dd);
                success = true;
                break;
            end
            if dd < maxdd
                options.restol = options.restol/sqrt(2);
            end
        end
        precnnz = AMGnnz_sym(PREC);
        fprintf('The density of L over that of A and MATLAB''s L is [%d %d].\n',...
        precnnz/nnz(A),precnnz/nnz(L));
        if ~success
            fprintf('PSQMR did not converge within %d steps ...\n',options.maxit);
            fprintf('Minimal residual %d reached at step %d.\n',relres,options_out.niter);
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
            fprintf('PSQMR converged at step %d. Relative res = %d.\n',options_out.niter,relres);
            %fprintf('Relative res = %d, backward res = %d\n',norm(b-A*x)/norm(b),relres);
            if precnnz < minnnzL
                minnnzL = precnnz;
                fprintf('This is the sparsest preconditioner found so far ...\n');
                %save(strcat(name,'.precond.mat'),'P_prpc','S_prpc','LL','DD','pp','qq','droptol');
            end
            dropworks = true;
            if mod(iter,2) == 1,	droptol = droptol*3/2;
	    else			droptol = droptol*4/3;
	    end
            iter = iter + 1;
        end
	fprintf('Time for PSQMR linear solve %.2f secs.\n',t_elapsed2);
	PREC = AMGdelete(PREC);
    end
    toc;
    warning on;
end
