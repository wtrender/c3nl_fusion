function teg_RMA_report(prestr, X0, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, y_red_M, verbose0, labelFull, Btmp, Btypes, contvec)

if verbose0 < 0,
    return;
end;

digprec1 = 3;
digprec2 = 2;
digprec3 = 2;

if ~isnan(F),
    resstr = [prestr labelFull ': F(' num2str(df1, digprec1) ', ' num2str(df2, digprec1) ') = ' num2str(F, digprec2) ', p = ' num2str(p, digprec3) ', eta_p^2 = ' num2str(eta, digprec2)];
    if verbose0 > 1,
        resstr = [resstr '\n' prestr  '\t\tSSM = ' num2str(SSM, digprec) '\tSSE = ' num2str(SSE, digprec) '\tMSM = ' num2str(MSM, digprec) '\tMSE = ' num2str(MSE, digprec)];
        resstr = [resstr '\t epsilon = ' num2str(eps0, digprec)];
    end;
    fprintf([prestr resstr]);
    if p <= 0.05,
        fprintf(' * ');
    end;
    fprintf('\n');
else,
    fprintf([prestr labelFull '\n']);
end;

if verbose0 <= 1,
    return;
end;

fprintf([prestr '\tDescriptives:\n']);

uB = unique(Btmp);
% Pairwise group differences per condition
if length(uB) > 1 && isempty(contvec),
    for iB1 = 1:length(uB),
        f1 = find(Btmp == uB(iB1));
        mcv1 = 1 * nanmean(y_red_M(f1, :), 2);
        for iB2 = 1:(iB1 - 1),
            f2 = find(Btmp == uB(iB2));
            mcv2 = 1 * nanmean(y_red_M(f2, :), 2);
            for iC = 1:size(y_red_M, 2),
                v1 = y_red_M(f1, iC);
                v2 = y_red_M(f2, iC);
                if size(y_red_M, 2) > 1,
                    v1 = v1 - mcv1;
                    v2 = v2 - mcv2;
                end;
                m1 = nanmean(v1);
                m2 = nanmean(v2);
                [H,P,CI,STATS] = ttest2(v1, v2);
                fprintf(['\tGroup ' num2str(iB1, digprec) ' - ' num2str(iB2, digprec) ' W = ' num2str(iC) ': ' num2str(m1) ' - ' num2str(m2) ', t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P) ' (W0)\n']);
                for iC2 = 1:(iC - 1),
                    v1b = y_red_M(f1, iC2);
                    v2b = y_red_M(f2, iC2);
                    if size(y_red_M, 2) > 1,
                        v1b = v1b - mcv1;
                        v2b = v2b - mcv2;
                    end;
                    m1 = nanmean(v1) - nanmean(v1b);
                    m2 = nanmean(v2) - nanmean(v2b);
                    [H,P,CI,STATS] = ttest2(v1 - v1b, v2 - v2b);
                    fprintf(['\tGroup ' num2str(iB1, digprec) ' - ' num2str(iB2, digprec) ' DW = ' num2str(iC) ' - ' num2str(iC2) ': ' num2str(m1) ' - ' num2str(m2) ', t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P) '\n']);
                end;
            end;
        end;
    end;
end;

for iGroup = 0:size(uB, 1),
    if iGroup > 0,
        fprintf(['\tGroup ' num2str(iGroup, digprec) ':\t group-code = ' num2str(Btypes(iGroup, :), digprec) '\n']);
        fg = find(Btmp == uB(iGroup));
    end;
    if iGroup == 0,
        if ~isempty(Btmp),
            continue;
        else
            fprintf('\t');
            fg = 1:size(y_red_M, 1);
        end;
    end;
    y_red_Mg = y_red_M(fg, :);
    y_red_Mw = y_red_Mg - mean(y_red_Mg, 2) * ones(1, size(y_red_Mg, 2));
    if isempty(contvec),
        f = ~isnan(mean(y_red_Mg, 2));
        myr = mean(y_red_Mg(f, :));
        if ~isempty(X0),
            [pvec00, t] = teg_ttest(y_red_Mw(f, :));
        else
            [pvec00, t] = teg_ttest(y_red_Mg(f, :));
        end;
        %         myr = t;
        df00 = length(f) - 1;
        fprintf(['\tMeans (df = ' num2str(df00, digprec) '):\t']);
    else
        f = ~isnan(contvec(fg)) & ~isnan(mean(y_red_Mg, 2));
        [myr, pvec00] = corr(y_red_Mg(f, :), contvec(fg(f)));
        df00 = length(f) - 1;
        fprintf(['\tCorrelations (df = ' num2str(df00, digprec) '):\t']);
    end;
    for ib = 1:length(myr),
        fprintf([num2str(myr(ib), digprec)]);
        fprintf([' (p = ' num2str(pvec00(ib), digprec) ')']); % Note the p is for the subject-mean corrected values
        if ib < length(myr),
            fprintf('; ');
        else
            fprintf('\n');
        end;
    end;
    if ~isempty(contvec),
        for iC = 1:size(y_red_Mg, 2),
            for iC2 = 1:(iC - 1),
                d0 = y_red_Mg(:, iC) - y_red_Mg(:, iC2);
                [myr, pvec00] = corr(d0(f), contvec(fg(f)));
                fprintf(['\t\t' num2str(iC) '-' num2str(iC2) ': r = ' num2str(myr(1, 1)) ', p = ' num2str(pvec00(1, 1), digprec) '\n']);
            end;
        end;
    else,
        for iC = 1:size(y_red_Mg, 2),
            for iC2 = 1:(iC - 1),
                d0 = y_red_Mg(:, iC) - y_red_Mg(:, iC2);
                [p, t, df0, effs] = teg_ttest(d0);
                fprintf(['\t\t' num2str(iC) '-' num2str(iC2) ': t(' num2str(df0) ') = ' num2str(t) ', p = ' num2str(p, digprec) ', Cohen d = ' num2str(effs) '\n']);
            end;
        end;
    end;
end;

y_red_Mw = y_red_M - mean(y_red_M, 2) * ones(1, size(y_red_M, 2));
myr = mean(y_red_M);
seyr = var(y_red_Mw) .^ 0.5 / sqrt(size(y_red_M, 1));
if verbose0 > 0,
    if length(myr) > 2,
        % Paired t-tests
        fprintf([prestr '\tPost-hoc differences:\n']);
        for col1 = 1:size(y_red_M, 2),
            for col2 = 1:(col1 - 1),
                vec = y_red_M(:, col1) - y_red_M(:, col2);
                [p, t, df, effs] = teg_ttest(vec);
                if p < 0.05,
                    tresstr = [num2str(col1) ' - ' num2str(col2) ': t(' num2str(df) ') = ' num2str(t) ', p = ' num2str(p) ', Cohen d = ' num2str(effs)];
                    fprintf([prestr '\t\t' tresstr '\n']);
                end;
            end;
        end;
        % vs 0
        fprintf('\tPost-hoc test vs 0:\n');
        for col1 = 1:size(y_red_M, 2),
            vec = y_red_Mw(:, col1);
            [p, t, df] = teg_ttest(vec);
            if p < 0.05,
                tresstr = [num2str(col1) ': t(' num2str(df) ') = ' num2str(t) ', p = ' num2str(p)];
                fprintf([prestr '\t\t' tresstr '\n']);
            end;
        end;
    else,
        fprintf('\n');
    end;
    fprintf('\n');
else,
    fprintf('\n');
end;
