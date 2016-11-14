function [M,D,p] = spainv_sym_mwn(A,nmtype,beta,droptol,droptol_type)
%
% Function spainv_sym_mwn generates an inverse triangular factorization 
% of a real symmetric matrix, using a minimal weighted norm pivoting, a
% modified variant of generalized Bunch-Kaufman pivoting
% 
% Input: 
% A: a real symmetric matrix to be factorized
% beta: the threshold in (0,1] to determine the second candidate pivot 
%       row/column. beta = 1 corresponds to the standard rook pivoting. 
%       The smaller beta is, the larger bound on element growth but more 
%       room to choose a sparser pivot row/column
% droptol:  tolerance used to drop small elements in the returned M factor
% droptol_type:  whether the drop tolerance is absolute or relative  
%
% Output:
% M: unit upper triangular, conceptually equivalent to inv(L'), where L is
% the psychological unit lower triangular factor of A, i.e., A = L*D*L'.
% D: block diagonal matrix consisting of 1x1 and 2x2 diagonal blocks
% p: permutation vector incorporating all pivoting and dynamic reordering 
% such that M'*A(p,p)*M = D numerically if droptol = 0. 
%
%
% Main references:
% 
% J. R. Bunch and L. Kaufman, Some stable methods for calculating inertia
% and solving symmetric linear systems, Math. Comp., 31(137):163-179, 1977.
%
% M. Benzi, C. D. Meyer and M. Tuma, A sparse approximate inverse precondi-
% tioner for the conjugate gradient method, SIAM J. Sci. Comput., Vol. 17 
% (1996), pp 1135--1149.
%
%
% Main idea:
%
% Starting with M = I, spainv_sym_mwn takes n = length(A) steps to find 
% the conjugate unit column vectors for M. Before step dd, assume that dd-1
% such columns have been constructed and stored in M(:,p(1:dd-1)), and we 
% have A-orthogonalized the unprocessed vectors M(:,p(dd:n)) against the
% previously constructed conjugate vectors M(:,p(1:dd-1)).
%
% At step dd, we need to explore M_A = M(:,p(dd:n))'*A*M(:,p(dd:n)), i.e., 
% the restriction of A onto the space conjugate to the previously generated
% (processed) M(:,p(1:dd-1)). The matrix M_A correspond to A^{(k)} in the 
% Bunch-Kaufman paper, where k = n-dd+1.
%
% We follow Alg. A in the BK paper, reorder the columns of M(:,p(dd:n)),
% and choose s = 1 or 2 columns M_p = M(:,p(dd:dd-1+s)) from M(:,p(dd:n)), 
% so that D_s = M_p'*A*M_p is pivot of M_A (columns and rows reordered) in 
% the top-left corner. We then combine M_p with the processed columns
% M(:,p(1:dd-1)), D_s with the processed block diagonal D(1:dd-1,1:dd-1), 
% and A-orthogonalize the unprocessed columns M(:,p(dd+s:n)) against the 
% pivot columns M_p, update dd = dd-1+s, and proceed to the next step.

%%% Most expensive steps: computing the two candidate pivot row/columns Aii
%%% and Arr, A-orthogonalization, dropping small entries in the unprocessed
%%% columns, and count the nnz in each columns for the dynamic reordering 

fprintf('%s - AINV based on minimal weighted triangular norm.\n',mfilename);
if ~(beta >= 0.01 && beta <= 1)
    error('The threshold beta must fall in [0.01,1]\n');
end
normA = norm(A,1);
n = length(A);
normAFro = norm(A,'fro');
% initial unit upper triangular M and block diagonal D (stored as a tri-
% diagonal and hence can be constructed explicitly from the main diagonal 
% and subdiagonal, due to the real symmetry of A) 
M = speye(n,n);
DMain = zeros(n,1);
DSubd = zeros(n,1);
p = 1:n;
% the 'alpha' parameter used in gbk pivoting, minimizing the bound on
% element growth; also the unique positive root of the cubic shown below
alpha = max(roots([beta^4+3*beta^3 2*beta^2+beta -2*beta^2 -1]));
gamma = 1e-1;
eta = 1e-2;
droptol_max = 5e-1;
% flag indicating whether to drop small elements in the active columns of M
% during the conjugation procedure
opt_drop_small_new_fillin = true;
% frequency for dropping the small elements of the active columns of M
drop_new_fillin_step_size = 1;
drop_new_fillin_step_count = 0;
% frequency for performing a complete reordering of the active columns of M 
complete_reorder_step_size = floor(sqrt(n)/2);
complete_reorder_step_count = 0;
% pivot type: type 1 corresponds to pivot row/column 1 (size 1), type 2 
% corresponds to pivot row/column r (size 1), and type 3 corresponds to
% pivot rows/columns 1 and r (size 2)
pvttypesum = zeros(3,1);
progress = zeros(10,1);
fprintf('alpha = %d, beta = %d, increase factor not known.\n',alpha,beta);
dd = 1;
tic;
while dd <= n
    r = 0;
    % form the first pivot row/column of M_A = M(:,p(dd:n))'*A*M(:,pp(dd:n))
    A1 = ((A*M(:,p(dd)))'*M(:,p(dd:end)))';
    A1_skip1 = A1;  A1_skip1(1) = 0;
    % use the threshold beta to find the second candidate pivot row/column
    if dd < n
        omega_1 = norm(A1_skip1,'inf');
        r = find(abs(A1_skip1) >= beta*omega_1,1);
    else
        omega_1 = 0;
    end
    % `half-machine-precision' regularization of a zero coordinate pivot
    % row/column; artificially set the diagonal element half epsilon and 
    % accept this regularized coordinate pivot row/column
    if max([abs(A1(1)) omega_1]) <= 64*sqrt(eps)*normA
        DMain(dd) = sign(0.5+sign(real(A1(1))))*64*sqrt(eps)*normA;
        A1(1) = DMain(dd);
        s = 1;      pvttype = 1;    pvttypesum(1) = pvttypesum(1)+1;
    % the leading 1x1 pivot is greater than or equal to alpha times the
    % selected off-diagonal element corresponding to the second candidate
    % coordinate pivot row/column (which itself is greather than or equal 
    % to beta times the largest off-diagonal element)
    elseif abs(A1(1)) >= alpha*beta*omega_1
        DMain(dd) = A1(1);   s = 1;
        pvttype = 1;    pvttypesum(1) = pvttypesum(1)+1;
    % otherwise, the second candidate coordinate pivot row/column must be 
    % formed for further comparison and pivot selection;
    else
        %  mwn pivoting compares the norm of the coordinate vector of each
        %  candidate pivot (normC1, normCr and normC1r), and the product of
        %  the norm of the coordiante vector and the norm of the pivot row
        %  (ptb_nm1, ptb_nmr and ptb_nm1r). It assigns each quantity with a
        %  different weight, and chooses the pivot associated with minimal 
        %  weighted quantity
        Ar = ((A*M(:,p(dd-1+r)))'*M(:,p(dd:end)))';
        Ar_skipr = Ar;      Ar_skipr(r) = 0;
        normM1 = norm(M(:,p(dd)),nmtype);
        A1_skip1_nm = norm(A1_skip1,nmtype);
        normC1 = A1_skip1_nm/abs(A1(1));
        ptb_nm1 = normM1*normC1;
        normMr = norm(M(:,p(dd-1+r)),nmtype);
        Ar_skipr_nm = norm(Ar_skipr,nmtype);
        normCr = Ar_skipr_nm/abs(Ar(r));
        ptb_nmr = normMr*normCr;
        % an inexpensive approach to compute normC1r and ptb_nm1r if 2-norm
        % is used
        if nmtype == 2
            C11 = (A1_skip1_nm)^2-abs(A1(r))^2;
            Crr = (Ar_skipr_nm)^2-abs(Ar(1))^2;
            A1_skip1(r) = 0;    Ar_skipr(1) = 0;
            C1r = sum(conj(A1_skip1).*Ar_skipr,1);
            Pvt = [A1(1) Ar(1); A1(r) Ar(r)];
            Csqrd = Pvt'\full([C11 C1r; C1r' Crr])/Pvt;
            if isnan(norm(Csqrd,1)) || isinf(norm(Csqrd,1))
                normC1r = Inf;
                ptb_nm1r = Inf;
            else
                M11 = norm(M(:,p(dd)))^2;   Mrr = norm(M(:,p(dd-1+r)))^2;
                M1r = M(:,p(dd))'*M(:,p(dd-1+r));
                normM1r = sqrt(max(eig([M11 M1r; M1r' Mrr])));
                normC1r = sqrt(max(eig(Csqrd)));
                ptb_nm1r = normM1r*normC1r;
            end
        else
            C1r = [A1([2:r-1 r+1:end]) Ar([2:r-1 r+1:end])]/[A1(1) Ar(1); A1(r) Ar(r)];
            normC1r = norm(C1r,nmtype);
            if isnan(normC1r) || isinf(normC1r)
                normC1r = Inf;
                ptb_nm1r = Inf;
            else
                normM1r = norm(M(:,p([dd dd-r+1])),nmtype);
                ptb_nm1r = normM1r*normC1r;
            end
        end
        % Since the single pivot column r is least favorable in terms of
        % the number of new fill-ins it introduces in A-orthogonalization,
        % this pivot column is used only if coordinate pivot column r is 
        % much closer to diagonally dominant than coordiante pivot column 1
        % or columns 1 and r (i.e., normCr is much smaller than normC1 and 
        % normC1r, or ptb_nmr is much smaller than ptb_nm1 and ptb_nm1r). 
        % In this case, using this pivot column introduces a much smaller 
        % orthogonalization corretion to the unprocessed columns of M than 
        % the other two choices of pivot columns, and it is likely to drop 
        % more small elements of the unprocessed columns to offset the
        % relatively high number of new fill-ins introduced
        if ptb_nmr < eta*min([ptb_nm1*gamma ptb_nm1r]) || normCr < eta*min([gamma*normC1 normC1r])
            tmp = p(dd);     p(dd) = p(dd+r-1);    p(dd+r-1) = tmp;
            DMain(dd) = Ar(r);   s = 1;
            pvttype = 2;    pvttypesum(2) = pvttypesum(2)+1;
        % the 2x2 pivot corresponding to coordiante pivot columns 1 and r
        % is of medium favorability. It is chosen if the single pivot
        % column r is not much closer diagonally dominant than the other
        % two pivot options, and if the coordinate pivot columns 1 and r
        % are obviously closer to diagonally dominant than pivot column 1
        elseif ptb_nm1r < gamma*ptb_nm1 || normC1r < gamma*normC1
            tmp = p(dd+1);   p(dd+1) = p(dd+r-1);  p(dd+r-1) = tmp;
            DMain(dd) = A1(1);   DMain(dd+1) = Ar(r);
            DSubd(dd) = A1(r);   s = 2;
            pvttype = 3;    pvttypesum(3) = pvttypesum(3)+1;
        % if coordinate pivot column r and pivot columns 1 and r both fail
        % to be sufficiently closer to diagonally dominant than other
        % options, then we surely choose the most favorable pivot column 1
        % which typically introduces the smallest number of fill-ins
        else
            DMain(dd) = A1(1);  s = 1;
            pvttype = 1;    pvttypesum(1) = pvttypesum(1)+1;
        end
    end
    M_pvtcols = M(:,p(dd:dd-1+s));
    
    % dropping small elements in the pivot columns, one column at a time
    for jj = 1 : s
        % assemble all nonzero elements in the current pivot column
        [nzidx,~,tmpv] = find(M_pvtcols(p(1:dd),jj));
        if ~isempty(nzidx)
            if strcmpi(droptol_type,'relative')
                % sort the nonzero elements in the current pivot column,
                % and drop the smallest ones whose sum is no greater than
                % droptol times the total sum of all nonzero elements; 
                % (note that we could sum up the elements themselves or
                % of their squares)
                [ordtmpv,ordix] = sort(conj(tmpv).*tmpv);   %sort(abs(tmpv));  % 
                cumsum_ordtmpv = cumsum(ordtmpv);
                dropidx = min([find(cumsum_ordtmpv >= cumsum_ordtmpv(end)*droptol,1) ...
                               find(ordtmpv >= droptol_max,1)])-1;
                if dropidx > 0
                    tmpv(ordix(1:dropidx)) = 0;
                    M(p(nzidx),p(dd-1+jj)) = tmpv;
                end
            else
                % absolute drop tolerance means to simply drop all nonzero 
                % elements smaller than droptol
                filter = min([droptol droptol_max]);
                M(p(nzidx),p(dd-1+jj)) = tmpv.*(abs(tmpv) >= filter);
            end
        end
    end
    
    if dd+s > n
        break;
    end
    
    % Orthgonalize the unprocessed vectors against the newly constructed 
    % pivot column vectors M(:,p(dd:dd-1+s)). Note here that we use A1 and
    % Ar formed using the corresponding pivot columns of M before dropping
    % the small elements of these columns.
    
    % Specifically, let M_p = M(:,p(dd:dd-1+s)) be the pivot column(s), and
    % M_u = M(:,p(dd+s:n)) be the unprocessed columns (pivot columns have
    % been reordered to the beginning of the active block of M). Then 
    % A1 = M(:,p(dd))'*A*M(:,p(dd:n)), Ar = M(:,p(dd+1))'*A*M(:,p(dd:n))
    % and PvtCols(1:s,1:s) = M_p'*A*M_p, PvtCols(1:s,s+1:end) = M_p'*A*M_u.
    % Then M_u = M_u - M_p*inv(M_p'*A*M_p)*(M_p'*A*M_u) is the Gram-Schmidt
    % A-orthogonalization step. 
    
    % To make the A-orthogonalization step most efficient, we find the 
    % nonzero columns in the coordinate pivot columns M_p'*A*M_u (or 
    % equivalently, inv(M_p'*A*M_p)*(M_p'*A*M_u)), and simply update the
    % columns of M_u corresponding to these nonzero columns; other columns
    % of M_u are already A-orthogonal to M_p and need not be updated 
    
    if pvttype == 1
        PvtCols = A1';
    elseif pvttype == 2
        PvtCols = Ar';  PvtCols([1 r]) = PvtCols([r 1]);
    else
        PvtCols = [A1 Ar]'; PvtCols(:,[2 r]) = PvtCols(:,[r 2]);
    end
    codnt_mag = sum(abs(PvtCols(:,s+1:end)),1);
    idx = find(codnt_mag >= eps/2*norm(codnt_mag,'inf'));
    GS_coeffs = PvtCols(:,1:s)\PvtCols(:,s+idx);
    updt_col_idx = (dd+s):n;
    updt_col_idx = updt_col_idx(idx);
    M(:,p(updt_col_idx)) = M(:,p(updt_col_idx)) - M_pvtcols*GS_coeffs;
    dd = dd + s;
    
    % After A-orthogonalization, drop the small elements in M_u to keep it
    % relatively sparse. If relative drop tolerance is used, we drop all
    % elements smaller than 0.1*droptol times the square of the 2-norm of
    % the newly-updated columns of M_u from these columns; otherwise we
    % drop those smaller than a quarter of droptol in modulus (a heuristic)
    
    if opt_drop_small_new_fillin
        if strcmpi(droptol_type,'relative')
            M_colsum = sum(conj(M(:,p(updt_col_idx))).*M(:,p(updt_col_idx)),1);
            % sum(abs(M(:,p(updt_col_idx))),1);
            gmean_M_colsum = min(M_colsum); %exp(sum(log(M_colsum))/length(updt_col_idx));  % 
            threshold = min([droptol*gmean_M_colsum/10 droptol_max]);
        else
            threshold = min([droptol/4 droptol_max]);
        end
        %%% The dropping of small elements of the unprocessed columns is 
        %%% rather time consuming; however, less frequent dropping seems to 
        %%% increase the density and deterioriate quality of preconditioner
        if floor(dd/drop_new_fillin_step_size) > drop_new_fillin_step_count
            drop_new_fillin_step_count = floor(dd/drop_new_fillin_step_size);
            M(:,p(updt_col_idx)) = M(:,p(updt_col_idx)).*(abs(M(:,p(updt_col_idx))) >= threshold);
        end
    end
    
    % After dropping small elements in M_u, we perform a dynamic reordering
    % of the columns of M_u. A complete reordering of all columns of M_u is
    % performed every complete_reorder_step_count steps of the inverse 
    % factorization; otherwise the leading r+1 columns of M_u are reordered
    
    if floor(dd/complete_reorder_step_size) > complete_reorder_step_count
        complete_reorder_step_count = floor(dd/complete_reorder_step_size);
        num_cols_reorder = n-dd+1;
    else
        num_cols_reorder = min([r+1 n-dd+1]);
    end
    
    % Reorder the columns of M_u according to their number of nonzeros. 
    % The pool of candidate columns contains all columns of M_u reordered.
    M_spones = M(:,p(dd:dd-1+num_cols_reorder)) ~= 0;
    nnz_column = sum(M_spones,1);
    nnz_row = sum(M_spones,2);
    [nnz_column,idx] = sort(nnz_column);
    p(dd:dd-1+num_cols_reorder) = p(dd-1+idx);
    
    % Then, choose the sparse columns from the pool that are marginally
    % denser than the sparsest column. Specifically, if the densest column
    % in the pool has 2-10 nonzeros, then the columns that are denser than 
    % the sparsest column by 1 element are considered; if the denest column
    % has 11-100 nonzeros, then those denser than the sparsest by 2 entries
    % are considered, and so on (a heuristic). 
    % Among these sparse columns of M_u, the one with the largest inf-norm 
    % is placed at the leading position. Such a sparse column of relatively 
    % large norm is likely to become a single pivot column at the next step
    % of factorization, producing essentially a minimal number of fill-ins
    % during the A-orthogonalization step
    
    nscr = find(nnz_column >= nnz_column(1)+ceil(log10(nnz_column(end)))+1,1)-1;
    if isempty(nscr),   nscr = num_cols_reorder;    end
    nscr = min([nscr 100]);
    Mdg = max(abs(M(:,p(dd:dd-1+nscr))),[],1);
    [~,idx] = sort(Mdg,'descend');
    p(dd:dd-1+nscr) = p(dd-1+idx);
    nscr = min([ceil(nscr/10) n-dd+1]);
    fillin = zeros(nscr,1);
    for zz = 1 : nscr
        fillin_row_idx = find(M(:,p(dd-1+zz)));
        fillin(zz) = num_cols_reorder*length(fillin_row_idx) - sum(nnz_row(fillin_row_idx));
    end
    [~,fillin_idx] = sort(fillin);
    p(dd:dd-1+nscr) = p(dd-1+fillin_idx);
    
    % Every 0.1*n factorization steps, report the progress, including the 
    % current density and the 1-norm of M, the number of each type of pivot 
    % used so far, and relative residual of the completed factorization
    tens_percent_done = floor(dd/n*10);
    if tens_percent_done > 0 && tens_percent_done < 10 && progress(tens_percent_done) == 0
        progress(tens_percent_done) = 1;
        nnzM = nnz(M);
        D = spdiags([DSubd DMain [0; conj(DSubd(1:end-1))]],-1:1,n,n);
        R = M(:,p(1:dd-1))'*A*M(:,p(1:dd-1))-D(1:(dd-1),1:(dd-1));
        fprintf(' %d0%% done. nnz(M) = %d (%.3f x nnz(A)). norm(M,1) = %.3d. Pvt type = [%d %d %d]. Rel res = %.3d.\n',...
            tens_percent_done,nnzM,nnzM/nnz(A),norm(M,1),[pvttypesum(1) pvttypesum(2) pvttypesum(3)],norm(R,'fro')/normAFro);
%         dgM = diag(M(p,p));
%         if norm(dgM-ones(n,1)) ~= 0
%             fprintf('Warning: some diagonal elements of M are missing!\n');
%         end
    end
end
nnzM = nnz(M);
% if nnz(sort(p)-(1:n)) > 0
%     error('Wrong permutation vector p\n');
% end
D = spdiags([DSubd DMain [0; conj(DSubd(1:end-1))]],-1:1,n,n);
fprintf('100%% done. nnz(M) = %d (%.3f x nnz(A)). norm(M,1) = %.3d. Pvt type = [%d %d %d]. Rel res = %.3d\n',...
    nnzM,nnzM/nnz(A),norm(M,1),[pvttypesum(1) pvttypesum(2) pvttypesum(3)],norm(R,'fro')/normAFro);
fprintf('Drop tolerance = %d\n',droptol);
M = M(p,p); % Unit upper triangular

end
