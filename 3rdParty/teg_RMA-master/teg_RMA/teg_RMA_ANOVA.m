function [F, dfM, dfE, p, SSM, SSE, MSM, MSE, eta, eps0, bModel, X_ind, dfM_adj, dfE_adj] = teg_RMA_ANOVA(M, X0, B, Bcoder, contvec, perm_test, nIts_perm)

% Handle non-within-subject tests not involving a between-group factor.
if isempty(X0) && isempty(Bcoder),
    F = NaN; dfM = NaN; dfE = NaN; p = NaN; SSM = NaN; SSE = NaN; MSM = NaN; MSE = NaN; eta = NaN; eps0 = NaN; 
    bModel = NaN; X_ind = NaN; dfM_adj = NaN; dfE_adj = NaN;
    y = mean(M')';
    if ~isempty(contvec),
        [C, P] = corr(y, contvec);
        p = P(1, 1);
    else,
        [p, t, df, effs] = teg_ttest(y, 0);
    end;
    return;
end;

nSubj = size(M, 1);
subjM = mean(M, 2);

% Here, the matrix is transformed to M_red in which variation not due to the factors
% involved in the current effecty is removed.
%
% The reduced matrix has the same size as the original, but columns within the same combination of levels will be
% identical.
number_of_dims = 1;
if ~isempty(X0),
    M_red = NaN * zeros(size(M));
    for iSubj = 1:size(M, 1),
        y = M(iSubj, :);
        y = y(:) - mean(y(:));
        
        X = X0;
        
        if length(find(X ~= 0)) > 0,
            b = inv(X'*X)*X'*y;
            model = X*b;
        else
            model = ones(size(y));
        end;
        M_red(iSubj, :) = model;
        
        number_of_dims = size(X, 2);
    end;
else,
    M_red = mean(M')';
end;

% The X0 matrix is expanded for the F-test by repeating for all subjects.
% Note that M_red will be vectorized to create the dependent variable y.
% X0 is nCols x nDummy and must be expanded to match this y.
if ~isempty(X0),
    newX0 = [];
    for iDummy = 1:size(X0, 2),
        tmp = X0(:, iDummy);
        tmp = tmp';
        tmp = ones(nSubj, 1) * tmp;
        newX0 = [newX0 tmp(:)];
    end
    X0 = newX0;
else,
    X0 = ones(size(M_red, 1), 1);
end;

% Encode interactions with the continuous variable.
if ~isempty(contvec),
    contvec = contvec - mean(contvec);
    origSize = size(X0);
    tmpX = reshape(X0, nSubj, length(X0(:)) / nSubj);
    tmpX = tmpX .* (contvec * ones(1, size(tmpX, 2)));
    tmpX = reshape(tmpX, origSize);
    X0 = tmpX;
end;

% Encode interactions with the between-group factor.
if ~isempty(Bcoder),
    Bdummy = teg_RMA_B_to_BX(B, Bcoder);
    
    % Adjust predictors
    nAdd = size(X0, 1) / size(Bdummy, 1);
    BdummyOrig = Bdummy;
    for iW = 2:nAdd,
        Bdummy = [Bdummy; BdummyOrig];
    end;
    Bdummy = teg_demean(Bdummy); % This prevents group size differences causing BxW interaction due to a W effect
    
    new_X0 = [];
    for iB = 1:size(Bdummy, 2),
        tmp = X0 .* (Bdummy(:, iB) * ones(1, size(X0, 2)));
        new_X0 = [new_X0 tmp];
    end;
    X0 = new_X0;
    
    X0 = X0 - ones(size(X0, 1), 1) * mean(X0);
end;

% The reduced matrix of dependent variables is now columnized and centered.
y = M_red(:) - mean(M_red(:));

% Acquire statistical parameters.
%
% If permutation tests are not used, an F-test is applied using
% Greenhouse-Geisser correction.
% If permutation tests are used, the nominal F-test values under the null hypothesis are empirically sampled.
if perm_test == 0,
    [F, dfM, dfE, p, SSM, SSE, MSM, MSE, eta, eps0, bModel, dfM_adj, dfE_adj] = inner(X0, y, M_red, number_of_dims, Bcoder, nSubj);
else,
    Fv = [];
    for iIt = 1:nIts_perm,
        if ~isempty(contvec),
            contvec = contvec - mean(contvec);
            origSize = size(X0);
            tmpX = reshape(X0, nSubj, length(X0(:)) / nSubj);
            tmpX = tmpX .* (contvec * ones(1, size(tmpX, 2)));
            tmpX = reshape(tmpX, origSize);
            X = tmpX;
        else,
            X = X0;
        end;
        % Permute y per independent observation
        y_perm_sgn = -1 + 2 * floor(2 * rand(nSubj, 1));
        y_perm_sgn = y_perm_sgn * ones(1, length(y)/length(y_perm_sgn));
        if size(y_perm_sgn, 2) > 1,
            dfdfd=0;
        end;
        y_perm_sgn = y_perm_sgn(:);
        y_perm = y .* y_perm_sgn;
        F = inner(X, y_perm, M_red, number_of_dims, Bcoder, nSubj, correct_for_cov, Cont);
        Fv = [Fv; F];
    end;

    if ~isempty(contvec),
        contvec = contvec - mean(contvec);
        origSize = size(X0);
        tmpX = reshape(X0, nSubj, length(X0(:)) / nSubj);
        tmpX = tmpX .* (contvec * ones(1, size(tmpX, 2)));
        tmpX = reshape(tmpX, origSize);
        X = tmpX;
    else,
        X = X0;
    end;
    
    [F, dfM, dfE, p, SSM, SSE, MSM, MSE, eta, eps0, bModel, dfM_adj, dfE_adj] = inner(X, y, M_red, number_of_dims, Bcoder, nSubj);
    p = length(find(Fv >= F)) / length(Fv);
end;

function [F, dfM, dfE, p, SSM, SSE, MSM, MSE, eta, eps0, bModel, dfM_adj, dfE_adj] = inner(X, y, M_red, number_of_dims, Bcoder, nSubj)

y = y - mean(y);
f = find(min(var(X)) > 0);
if length(f) > 0,
    X(:, f) = teg_demean(X(:, f));
end;
if isnan(rcond(X' * X)),
    b = zeros(size(X, 2), 1);
else,
    b = inv(X'*X)*X'*y;
end;
model = X*b;
err = y - model;
SSM = sum(model.^2);
SSE = sum(err.^2);

bModel = b;

% Calculate epsilon
if size(M_red, 2) > 1,
    C = cov(M_red);
    [O2L, L] = eig(C);
    L = diag(L);
    
    % Imprecision: "Zero" values in terms of dimensionality can still well above the machine-minimum eps
    buffer0_zero = max(L) * (10 * eps);
    
    fnonzero = find(L > buffer0_zero);
    d = length(fnonzero);
    
    eps0 = (sum(L) .^ 2) / (d * sum(L .^ 2));
else
    eps0 = 1;
end;

nGroups = size(Bcoder, 2) + 1;
dfM = size(X, 2);
dfE = (nSubj - 1 - max(0, (nGroups - 1))) * (number_of_dims); % Take account of matrix trickery with Bcoder
dfM_adj = eps0 * dfM;
dfE_adj = eps0 * dfE;
MSM = SSM / dfM_adj;
MSE = SSE / dfE_adj;
F = MSM / MSE;
p = teg_fsig(F, dfM_adj, dfE_adj);
eta = SSM / (SSM + SSE);

function thiscorrvec = inner_corr_pre_gr(y, contvec, B, Bcoder)
Bdummy = teg_B_to_BX(B, Bcoder);
tmp = 1;
for n = 2:size(Bdummy, 2),
    tmp = [tmp; 10.^(n-1)];
end;
BGr = Bdummy * tmp;
uB = unique(BGr);
thiscorrvec = [];
for iu = 1:length(uB),
    f = find(BGr == uB(iu));
    cf = contvec(f);
    yf = y(f);
    c = corr(cf, yf);
    thiscorrvec = [thiscorrvec; c];
end;
