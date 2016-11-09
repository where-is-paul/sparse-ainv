function [L,D,p] = spildl_sym_all(A,normtype,beta,droptol,droptol_type,pvtopt)

% Another variant based on actual application of candidate pivot to the
% corresponding off-diagonal entries
if ~(beta >= 0.01 && beta <= 1)
    error('The threshold beta must fall in [0.01,1]\n');
end

n = length(A);
p = 1:n;
L = speye(size(A));
DMain = zeros(n,1);
DSubd = zeros(n,1);

B = A;
normA = norm(A,1);
if isreal(normtype),    strnmtype = num2str(normtype);
else                    strnmtype = normtype;
end

if droptol < 1024*eps,  droptol = 1024*eps;     end

if strcmpi(pvtopt,'grook')
    alpha = max(roots([4*beta 3 -2*beta -1]));
else
    alpha = max(roots([beta^4+3*beta^3 2*beta^2+beta -2*beta^2 -1]));
end

% used with weighted minimal norm pivoting to assign different level of
% priority for the two 1x1 pivots and the 2x2 pivot
gamma = 1e-1;       eta = 1e-2;
% maximum magnitude of elements that could be dropped
droptol_max = 5e-1;
% relative density threshold of the Schur complement; switch to MA57 once
% SC becomes denser than this threshold
denthreshold = 0.5;
% drop small elements of SC periodically
drop_SCsmall_everynumsteps = 100;
drop_SCsmall_times = 0;
% reorder _all_ columns of SC periodically; partial reordering is done
% between two complete reorderings
reorder_stepsize = floor(sqrt(n)/2);
num_reorder_step = 0;
% clear the nonactive part of B periodically
print_stepsize = 250;
num_print_step = 0;
% sparsest column first enabled by default
sparfirst = (1 > 0);
progress = zeros(10,1);
% maximum pivot size and total number of pivots
act_maxpvt_size = 0;
pvt_count = 0;
% smallest and largest singular values of the block diagonal matrix D
small_sv_D = realmax;
large_sv_D = realmin;
condD = max(1,large_sv_D/small_sv_D);
pvttypesum = zeros(3,1);
% an additional permutation aiming at approximate minimal fill-ins
minfill = (0 > 0);

if strcmpi(pvtopt,'gbk')
    fprintf('%s - ILDL based on generalized Bunch-Kaufman pivoting.\n', mfilename);
    fprintf('alpha = %d, beta = %d, increase factor = %d.\n',alpha,beta,1+1/(alpha*beta^2));
elseif strcmpi(pvtopt,'grook')
    fprintf('%s - ILDL based on generalized rook pivoting.\n', mfilename);
    fprintf('alpha = %d, beta = %d, increase factor = %d.\n',alpha,beta,1+1/(alpha*beta));
else
    fprintf('%s - ILDL based on minimal weighted norm pivoting (%d-norm).\n', mfilename,normtype);
    fprintf('alpha = %d, beta = %d, increase factor not known.\n',alpha,beta);
    fprintf('Thresholds [Gamma Eta] for pivoting = [%.3d %.3d].\n',gamma,eta);
end

if sparfirst
    fprintf('Sparsest column first + partial diagonal pivoting enabled.\n');
else
    fprintf('Sparsest column first + partial diagonal pivoting disabled.\n');
end
if minfill
    fprintf('Approximate minimal fill-in reordering enabled.\n');
else
    fprintf('Approximate minimal fill-in reordering disabled.\n');
end

dd = 1;
while dd < n
    
    ii = 1;
    rr = 0;
    searchdepth = n-dd+1;
    
    Bii = B(p(dd:end),p(dd));
    Bii_skipi = Bii;  Bii_skipi(ii) = 0;
    if dd < n
        omega_ii = norm(Bii_skipi,'inf');
        rr = find(abs(Bii_skipi) >= beta*omega_ii,1);
    else
        omega_ii = 0;
    end
    % 'half-machine-precision' regularization of the leading diagonal pivot
    % if the whole leading pivot column has very small elements
    if max([abs(Bii(ii)) omega_ii]) <= 64*sqrt(eps)*normA
        p([dd dd-1+ii]) = p([dd-1+ii dd]);
        DMain(dd) = sign(0.5+sign(real(Bii(ii))))*64*sqrt(eps)*normA;
        B(p(dd),p(dd)) = DMain(dd);
        s = 1;      pvttypesum(1) = pvttypesum(1)+1;
        % the leading 1x1 pivot is greater than alpha times the selected
        % off-diagonal element corresponding to the second candidate pivot
        % row/column (which itself is greather than beta times the largest
        % off-diagonal element)
    elseif abs(Bii(ii)) >= alpha*beta*omega_ii
        p([dd dd-1+ii]) = p([dd-1+ii dd]);
        DMain(dd) = Bii(ii);   s = 1;
        pvttypesum(1) = pvttypesum(1)+1;
        % otherwise, the second candidate pivot row/column must be formed
        % for further comparison and pivot selection
    else
        if strcmpi(pvtopt,'grook')
            searchstep = 2;
            while searchstep <= searchdepth
                rr = find(abs(Bii_skipi) >= beta*omega_ii,1);
                Brr = B(p(dd:end),p(dd-1+rr));
                Brr_skipr = Brr;    Brr_skipr(rr) = 0;
                omega_rr = norm(Brr_skipr,'inf');
                % 'half-machine-precision' regularization of the second
                % pivot column (column r)
                if max([abs(Brr(rr)) omega_rr]) <= 64*sqrt(eps)*normA
                    p([dd dd-1+rr]) = p([dd-1+rr dd]);
                    DMain(dd) = sign(0.5+sign(real(Brr(rr))))*64*sqrt(eps)*normA;
                    B(p(dd),p(dd)) = DMain(dd);
                    s = 1;      pvttypesum(2) = pvttypesum(2)+1;
                    break;
                    % the diagonal element of the second candidate pivot
                    % column (i.e., the (r,r) element) is relatively large
                elseif abs(Brr(rr)) >= alpha*beta*omega_rr
                    p([dd dd-1+rr]) = p([dd-1+rr dd]);
                    DMain(dd) = Brr(rr);   s = 1;
                    pvttypesum(2) = pvttypesum(2)+1;
                    break;
                    % the (i,r) and (r,i) elements are relatively large
                elseif abs(Brr(ii)) >= beta*omega_rr
                    if ii < rr
                        pvt3_tmpidx = [dd:dd-1+ii-1 dd-1+ii+1:dd-1+rr-1 dd-1+rr+1:n];
                        pvt3_tmpval = p(pvt3_tmpidx);
                    else
                        pvt3_tmpidx = [dd:dd-1+rr-1 dd-1+rr+1:dd-1+ii-1 dd-1+ii+1:n];
                        pvt3_tmpval = p(pvt3_tmpidx);
                    end
                    p(dd:dd+1) = p([dd-1+ii dd-1+rr]);
                    p(dd+2:end) = pvt3_tmpval;
                    DMain(dd) = Bii(ii);   DMain(dd+1) = Brr(rr);
                    DSubd(dd) = Brr(ii);   s = 2;
                    pvttypesum(3) = pvttypesum(3)+1;
                    break;
                    % none of the (i,i), or (r,r), or (i,r) + (r,i) entries
                    % are relatively large to become the pivot; continue 
                    % searching from the current candidate pivot column r
                else
                    ii = rr;    omega_ii = omega_rr;
                    Bii = Brr;  Bii_skipi = Brr_skipr;
                    searchstep = searchstep + 1;
                end
            end
        elseif strcmpi(pvtopt,'gbk')
            Brr = B(p(dd:end),p(dd-1+rr));
            Brr_skipr = Brr;      Brr_skipr(rr) = 0;
            omega_rr = norm(Brr_skipr,'inf');
            % the leading 1x1 pivot ((1,1) element) is relatively large
            if abs(Bii(ii))*omega_rr >= alpha*beta^2*omega_ii^2
                DMain(dd) = Bii(ii);   s = 1;
                pvttypesum(1) = pvttypesum(1)+1;
                % the second 1x1 pivot ((r,r) element) is relatively large
            elseif abs(Brr(rr)) >= alpha*beta*omega_rr
                p([dd dd-1+rr]) = p([dd-1+rr dd]);
                DMain(dd) = Brr(rr);   s = 1;
                pvttypesum(2) = pvttypesum(2)+1;
                % both the leading and the second 1x1 pivots are relatively 
                % small compared to the (1,r) (r,1) elements; use 2x2 pivot
            else
                pvt3_tmpidx = [dd+1:dd-1+rr-1 dd-1+rr+1:n];
                pvt3_tmpval = p(pvt3_tmpidx);
                p(dd+1) = p(dd-1+rr);
                p(dd+2:end) = pvt3_tmpval;
                DMain(dd) = Bii(ii);   DMain(dd+1) = Brr(rr);
                DSubd(dd) = Bii(rr);   s = 2;
                pvttypesum(3) = pvttypesum(3)+1;
            end
        else
            %  mwn pivoting compares the norm of the off-diagonal part of
            %  the candidate columns of L associated with each pivot choice
            %  (normC1, normCr and normC1r). It assigns each quantity with
            %  a different weight, and chooses the pivot associated with 
            %  minimal weighted quantity
            
            Brr = B(p(dd:end),p(dd-1+rr));
            Brr_skipr = Brr;      Brr_skipr(rr) = 0;
            Bii_skipi_nm = norm(Bii_skipi,normtype);
            normCii = Bii_skipi_nm/abs(Bii(ii));
            Brr_skipr_nm = norm(Brr_skipr,normtype);
            normCrr = Brr_skipr_nm/abs(Brr(rr));
            % an inexpensive approach to compute normC1r if 2-norm is used
            if normtype == 2
                Cii = (Bii_skipi_nm)^2-abs(Bii(rr))^2;
                Crr = (Brr_skipr_nm)^2-abs(Brr(ii))^2;
                Bii_skipi(rr) = 0;    Brr_skipr(ii) = 0;
                Cir = sum(conj(Bii_skipi).*Brr_skipr,1);
                Pvt = [Bii(ii) Brr(ii); Bii(rr) Brr(rr)];
                Csqrd = Pvt'\full([Cii Cir; Cir' Crr])/Pvt;
                Csqrd = (Csqrd+Csqrd')/2;
                if isnan(norm(Csqrd,1)) || isinf(norm(Csqrd,1))
                    normCir = Inf;
                else
                    normCir = sqrt(max(eig(Csqrd)));
                end
            else
                %Cir = [Bii([1:ii-1 ii+1:rr-1 rr+1:end]) Brr([1:ii-1 ii+1:rr-1 rr+1:end])]/[Bii(ii) Brr(ii); Bii(rr) Brr(rr)];
                Cir = [Bii([1:ii-1 ii+1:rr-1 rr+1:end]) Brr([1:ii-1 ii+1:rr-1 rr+1:end])]*...
                      ([Brr(rr) -Brr(ii); -Bii(rr) Bii(ii)]/(Bii(ii)*Brr(rr)-Brr(ii)*Bii(rr)));
                normCir = norm(Cir,normtype);
                if isnan(normCir) || isinf(normCir)
                    normCir = Inf;
                end
            end
            % The single pivot column r is least favorable in terms of the
            % number of new fill-ins it introduces in conjugation. Hence,
            % this pivot column is used only if it is much closer to 
            % diagonally dominant than coordiante pivot column 1 or columns
            % 1 and r (normCr is much smaller than normC1 and normC1r).
            % In this case, it is likely to drop more small elements of the
            % corresponding column of L to offset the relatively large 
            % number of newly introduced fill-ins
            if normCrr < eta*min([gamma*normCii normCir])
                p([dd dd-1+rr]) = p([dd-1+rr dd]);
                DMain(dd) = Brr(rr);   s = 1;
                pvttypesum(2) = pvttypesum(2)+1;
                % The 2x2 pivot corresponding to pivot columns 1 and r is 
                % of medium favorability. It is chosen if the single pivot
                % column r is not much closer diagonally dominant than the 
                % other two pivot options, and if the coordinate pivot 
                % columns 1 and r are considerably closer to diagonally 
                % dominant than pivot column 1
            elseif normCir < gamma*normCii
                pvt3_tmpidx = [dd+1:dd-1+rr-1 dd-1+rr+1:n];
                pvt3_tmpval = p(pvt3_tmpidx);
                p(dd+1) = p(dd-1+rr);
                p(dd+2:end) = pvt3_tmpval;
                DMain(dd) = Bii(ii);   DMain(dd+1) = Brr(rr);
                DSubd(dd) = Bii(rr);   s = 2;
                pvttypesum(3) = pvttypesum(3)+1;
                % If pivot column r and pivot columns 1 and r both fail to
                % be sufficiently closer to diagonally dominant than other
                % options, then we choose the most favorable pivot column 1
                % which typically introduces fewest fill-ins
            else
                DMain(dd) = Bii(ii);  s = 1;
                pvttypesum(1) = pvttypesum(1)+1;
            end
        end
    end
    
    pvt_count = pvt_count+1;
    pvtidx = dd:dd-1+s;
    nonpvtidx = dd+s:n;
    svdPvt = svd(full(B(p(pvtidx),p(pvtidx))));
    if min(svdPvt) < small_sv_D,    small_sv_D = min(svdPvt);   end
    if max(svdPvt) > large_sv_D,    large_sv_D = max(svdPvt);   end
    condD = large_sv_D/small_sv_D;
    % Record the largest pivot size actually used
    if s > act_maxpvt_size
        act_maxpvt_size = s;
    end
    
    % Retrieve the optimal pivot obtained from searching for such a pivot
    % from pivot size = 1 up to pvtsize.
    p_nonpvtidx = p(nonpvtidx);
    BlkDiagPvt = B(p(pvtidx),p(pvtidx));
    BlkOffDiag = B(p_nonpvtidx,p(pvtidx));
    % Apply the selected block diagonal pivot (LDL factorization)
    if s == 1
        updated_L_pvtcol = BlkOffDiag*(1/BlkDiagPvt);
    else
        updated_L_pvtcol = BlkOffDiag*([BlkDiagPvt(2,2) -BlkDiagPvt(1,2); -BlkDiagPvt(2,1) BlkDiagPvt(1,1)]/...
            (BlkDiagPvt(1,1)*BlkDiagPvt(2,2)-BlkDiagPvt(2,1)*BlkDiagPvt(1,2)));
    end
    % TL corresponds to the newly formed off-diagonal columns of L obtained
    % by applying a drop tolerance 10 times smaller; later we will use TL
    % to update the Schur complement
    TL = sparse(length(nonpvtidx),s);
    
    for jj = 1 : s
        [nzidx,~,tmpv] = find(updated_L_pvtcol(:,jj));
        if ~isempty(nzidx)
            if strcmpi(droptol_type,'relative')
                [ordtmpv,ordix] = sort(conj(tmpv).*tmpv);   %sort(abs(tmpv));  %
                cumsum_ordtmpv = cumsum(ordtmpv);
                dropidx = min([find(cumsum_ordtmpv >= cumsum_ordtmpv(end)*droptol,1) ...
                    find(ordtmpv >= droptol_max,1)])-1;
                tmpvL = tmpv;
                if dropidx > 0
                    tmpv(ordix(1:dropidx)) = 0;
                end
                L(p_nonpvtidx(nzidx),p(dd-1+jj)) = tmpv;
                dropidx = min([find(cumsum_ordtmpv >= cumsum_ordtmpv(end)*droptol/10,1) ...
                    find(ordtmpv >= droptol_max,1)])-1;
                if dropidx > 0
                    tmpvL(ordix(1:dropidx)) = 0;
                end
                TL(nzidx,jj) = tmpvL;
            else
                filter = min([droptol droptol_max]);
                L(p_nonpvtidx(nzidx),p(dd-1+jj)) = tmpv.*(abs(tmpv) >= filter);
                filter = min([droptol/10 droptol_max]);
                TL(nzidx,jj) = tmpv.*(abs(tmpv) >= filter);
            end
        end
    end
    % Restrict the update only to row/columns of the Schur complement that
    % need to be updated
    nzrowidx = find(sum(TL ~= 0,2));
    TL = TL(nzrowidx,:);
    p_nz_nonpvtidx = p(nonpvtidx(nzrowidx));
    B(p_nz_nonpvtidx,p_nz_nonpvtidx) = B(p_nz_nonpvtidx,p_nz_nonpvtidx)-TL*BlkDiagPvt*TL';
    
    dd = dd+s;
    if dd > length(A)
        break;
    end
    % Drop small elements from the Schur complement periodically
    if ceil(dd/drop_SCsmall_everynumsteps) > drop_SCsmall_times
        drop_SCsmall_times = ceil(dd/drop_SCsmall_everynumsteps);
        colsum = exp(mean(log(max(abs(B(p(dd:n),p(dd:n)))))));
        SCdroptol = min([colsum*droptol min(max(abs(B(p(dd:n),p(dd:n)))))/10])/100;
        B(p(dd:n),p(dd:n)) = B(p(dd:n),p(dd:n)).*(abs(B(p(dd:n),p(dd:n))) >= SCdroptol);
        % If, after the current dropping, the Schur complement is still
        % denser than the prescribed threshold, we do a more aggressive 
        % dropping and switch to exact LDL factorization of the SC, and 
        % then combine the completed incomplete LDL and the exact LDL of
        % the approximate SC
        if nnz(B(p(dd:n),p(dd:n))) >= (n-dd+1)^2*denthreshold
            D = spdiags([DSubd DMain [0; conj(DSubd(1:end-1))]],-1:1,n,n);
            R = A(p(1:dd-1),p(1:dd-1))-L(p(1:dd-1),p(1:dd-1))*D(1:dd-1,1:dd-1)*L(p(1:dd-1),p(1:dd-1))';
            fprintf('Step %d, r-density of SC >= %.3f, nnz(L) = %.3f x nnz(A), r-res = %.3d..\n',...
                dd,denthreshold,nnz(L)/nnz(A),norm(R,'fro')/norm(A,'fro'));
            fprintf('Dropping small entries of SC and switching to exact LDL (MA57).\n');
            colsum = exp(mean(log(max(abs(B(p(dd:n),p(dd:n)))))));
            SCdroptol = min([colsum*droptol min(max(abs(B(p(dd:n),p(dd:n)))))/10]);
            oldnnzSC = nnz(B(p(dd:n),p(dd:n)));
            B(p(dd:n),p(dd:n)) = B(p(dd:n),p(dd:n)).*(abs(B(p(dd:n),p(dd:n))) >= SCdroptol);
            newnnzSC = nnz(B(p(dd:n),p(dd:n)));
            fprintf('%d out of %d (%.2f%%) entries of SC dropped.\n',oldnnzSC-newnnzSC,oldnnzSC,(oldnnzSC-newnnzSC)/oldnnzSC*100);
            [EL,ED,Ep,ES] = ldl(B(p(dd:n),p(dd:n)),'vector');
            ES = ES(Ep,Ep);
            restidx = dd:n;      restidx = restidx(Ep);
            p(dd:n) = p(restidx);
            L(p(dd:n),p(dd:n)) = ES\EL;
            DMain(dd:n) = spdiags(ED,0);
            DSubd(dd:n) = spdiags(ED,-1);
            break;
        end
    end
    
    % Sparsest columns first ordering of the Schur complement (unprocessed 
    % part of B) is crucial for keeping L sparse.
    if sparfirst
        if floor(dd/reorder_stepsize) > num_reorder_step
            num_reorder_step = ceil(dd/reorder_stepsize);
            num_cols_reorder = n-dd+1;
        else
            num_cols_reorder = min([max([ii+1 rr+1]) n-dd+1]);
        end
        % Counting the number of nonzeros in each column of active part of
        % B is rather time-consuming in MATLAB. Any more efficient method?
        Z = B(p(dd:n),p(dd:dd-1+num_cols_reorder)) ~= 0;
        nnz_column = sum(Z,1);
        
        [nnz_column,idx] = sort(nnz_column);
        p(dd:dd-1+num_cols_reorder) = p(dd-1+idx);
        % Perform a local diagonal reordering of the unprocessed columns of B.
        % Here local means reodering the columns inside each block of the same
        % density (nnzs); this won't affect the sparsest-column-first strategy.
        %num_sparsest_col = find(([nnz_column realmax]-nnz_column(1))>0,1)-1;
        nscr = find(nnz_column >= nnz_column(1)+ceil(log10(nnz_column(end)))+1,1)-1;
        if isempty(nscr),   nscr = num_cols_reorder;    end
        nscr = min([nscr 100]);
        
        cols_pointer = dd:dd-1+nscr;
        diagB = abs(spdiags(B(p(cols_pointer),p(cols_pointer)),0));
        maxcolB = max(abs(B(p(cols_pointer),p(cols_pointer))))';
        [~,inneridx] = sort(diagB./maxcolB,'descend');
        p(cols_pointer) = p(cols_pointer(inneridx));
        
        if minfill
            Z = Z(:,inneridx);
            nnz_row = sum(Z,2);
            nscr = min([ceil(nscr/5) n-dd+1]);
            fillin = zeros(nscr,1);
            for zz = 1 : nscr
                fillin_row_idx = find(Z(:,zz));
                fillin(zz) = num_cols_reorder*length(fillin_row_idx) - sum(nnz_row(fillin_row_idx));
            end
            [~,fillin_idx] = sort(fillin);
            p(dd:dd-1+nscr) = p(dd-1+fillin_idx);
        end
    end
    
    if floor(dd/print_stepsize) > num_print_step
        num_print_step = floor(dd/print_stepsize);
        B(:,p(1:dd-1)) = 0;     B(p(1:dd-1),:) = 0;
    end
    % Report the progress
    tens_percent_done = floor(dd/n*10);
    if tens_percent_done > 0 && tens_percent_done < 10 && progress(tens_percent_done) == 0
        progress(tens_percent_done) = 1;
        nnzL = nnz(L);
        D = spdiags([DSubd DMain [0; conj(DSubd(1:end-1))]],-1:1,n,n);
        R = A(p(1:dd-1),p(1:dd-1))-L(p(1:dd-1),p(1:dd-1))*D(1:dd-1,1:dd-1)*L(p(1:dd-1),p(1:dd-1))';
        nnzB = nnz(B(p(dd:end),p(dd:end)));
        fprintf(' %d0%%: nnz(L) = %d (%.3f x nnz(A)), r-res = %.2d, cond(Dk) = %.2d, pvt type = [%d %d %d], nnz(SC) = %d (r-density = %.2d).\n',...
            tens_percent_done,nnzL,nnzL/nnz(A),norm(R,'fro')/norm(A,'fro'),condD,pvttypesum(1),pvttypesum(2),pvttypesum(3),nnzB,nnzB/(n-dd+1)^2);
    end
end

if dd == n
    DMain(n) = B(p(n),p(n));
    if min(abs(B(p(n),p(n)))) < small_sv_D
        small_sv_D = min(svdPvt);
    end
    if max(abs(B(p(n),p(n)))) > large_sv_D
        large_sv_D = max(svdPvt);
    end
    condD = large_sv_D/small_sv_D;
end

nnzL = nnz(L);
L = L(p,p);
D = spdiags([DSubd DMain [0; conj(DSubd(1:end-1))]],-1:1,n,n);

fprintf('100%%: nnz(L) = %d (%.3f x nnz(A)), r-res = %.2d, cond(Dk) = %.2d, pvt type = [%d %d %d], nnz(SC) = -- (r-density = --).\n',...
    nnzL,nnzL/nnz(A),norm(A(p,p)-L*D*L','fro')/norm(A,'fro'),condD,pvttypesum(1),pvttypesum(2),pvttypesum(3));

end