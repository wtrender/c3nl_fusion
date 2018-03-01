function O = teg_RMA(varargin)

% function O = teg_RMA(M, levels, varnames, Betw_group, Betw_labels, Cont, contLabels, factorLabelsOfInterest, argcode)
%
% Requires:     addpath 'D:\Dropbox\Dropbox\Work\Common\teg_basic_funcs'
%
% M is observation x variable-combination matrix.
% levels is vector of levels per factor.
% varnames is a cell array of strings.
%
% Betw_group: Matrix with columns containing categorical between-subject
% factors.
% Betw_labels: cell array of strings, matching the columns of Betw-group.
%
% Cont: Matrix with columns containing continuous between-subject
% factors.
% contLabels: cell array of strings, matching the columns of Cont.
%
% factorLabelsOfInterest: Only analyze effects involving exactly this
% combination of factor labels
%
% argcode: string containing various special characters to adjust the program.
%   'p': Use permutation testing
%   's': Silent, no output to screen
%   'a': Print all tests, regardless of significance
%   'd': Used for recursive analysis of interactions
%   'c': Correct p criterion for the number of effects being tested
%
% Greenhouse-Geisser correction is applied for parametric tests.
%
% Thomas E. Gladwin.
%
% Last update: 16/07/2017.

% Default parameters
betw_subj_interaction_depth = 1;
pCritForFurther = 0.05; % p-value for further post-hoc analyses
pStar = 0.05; % p-value to be printed with an asterisk
verbose0 = 2; % -1 to remove output, 1 or 2 to provide more extensive output
digprec = 3; % digital precision of printouts.
perm_test = 0;
nIts_perm = 2500;
all_str_offset = '';

% Parse other input arguments
[M, NM, levels, varnames, B, Bvarnames, ContRaw, contRawLabels, factorLabelsOfInterest] = teg_RMA_parse_arguments(varargin);

% Parameters based on input arguments
[nSubj, nVar] = size(M);
pStar3 = 0.05 / prod(levels); % p-value criterion to be printed with an asterisk

% Prepare continuous variables
if ~isempty(ContRaw),
    [Cont, contLabels, Cont_vars_involved] = teg_RMA_make_betw_cont(ContRaw, contRawLabels, betw_subj_interaction_depth);
else
    Cont = {};
end;

% Transform levels into dummy predictor matrix X
[X, factorStarts, nColsFactor, labels, cellsets, factors] = teg_RMA_create_ANOVA_dummy(levels, varnames, betw_subj_interaction_depth);

% Prepare dummy variables for Between-group effects
if isempty(B),
    bfactorStarts1 = [];
else,
    bLevels = [];
    for iB = 1:size(B, 2);
        u = unique(B(:, iB));
        bLevels(iB) = length(u);
    end;
    [bX1, bfactorStarts1, bnColsfactor1, blabels, bcellsets, bfactors] = teg_RMA_create_ANOVA_dummy(bLevels, Bvarnames);
end;

% Parse argcode if provided
argcode = '';
if ischar(varargin{end}),
    argcode = varargin{end};
    if ~isempty(findstr(argcode, 'p')),
        perm_test = 1;
    end;
    if ~isempty(findstr(argcode, 's')),
        verbose0 = -1;
    end;
    if ~isempty(findstr(argcode, 'a')),
        perm_test = 0;
        pCritForFurther = 1;
    end;
    if ~isempty(findstr(argcode, 'd')), % Step-down
        pCritForFurther = 3 * pCritForFurther;
        if length(varargin{2}) == 1,
            verbose0 = 2;
        else,
            verbose0 = 0;
        end;
        fd = length(findstr(argcode, 'd'));
        all_str_offset = '';
        for id = 1:length(fd),
            all_str_offset = [all_str_offset '> '];
        end;
    end;
    if ~isempty(findstr(argcode, 'c')),
        pCritForFurther = pCritForFurther / (length(factors) * (1 + size(Cont, 2)) * (1 + length(bfactors)));
    end;
end;

% Set up output structure with R(esults) and labels for tests
O.R = [];
O.labels = {};

% Run over all main effects and interaction effects; 0 is the overall average.
for iEffect = 0:length(factorStarts),

    % Get dummy variables for current effect
    if iEffect == 0,
        X0 = [];
        withinLabel = 'Subject-mean';
    else,
        predcols = factorStarts(iEffect):(factorStarts(iEffect) + nColsFactor(iEffect) - 1);
        X0 = X(:, predcols);
        withinLabel = labels{iEffect};
    end;
    
    % Loop over continuous and between-group factors
    continuousVarsInteractions = size(Cont, 2);
    for iCont = 0:continuousVarsInteractions,
        if iCont == 0,
            contvec = [];
            contstr = '';
        else,
            contvec = Cont(:, iCont);
            contstr = [' x ' contLabels{iCont}];
        end;
        
        for iBetwInt = 0:length(bfactorStarts1),
            betwStr = '';
            if iBetwInt == 0,
                Bcoder = [];
            else,
                a = bfactorStarts1(iBetwInt);
                b = a + bnColsfactor1(iBetwInt) - 1;
                Bcoder = bX1(:, a:b);
                betwStr = [' x ' blabels{iBetwInt}];
            end;
            
            % Create factor labels cell array
            factorsLabels_this = {};
            factorsLabels_this_except_last = {};
            if iEffect > 0,
                for iFactor = 1:length(factors{iEffect}),
                    factorsLabels_this{end + 1} = varnames{factors{iEffect}(iFactor)};
                    if iFactor < length(factors{iEffect}),
                        factorsLabels_this_except_last{end + 1} = varnames{factors{iEffect}(iFactor)};
                    end;
                end;
            end;
            factorsLabels_this{end + 1} = contstr;
            factorsLabels_this_except_last{end + 1} = contstr;
            factorsLabels_this{end + 1} = betwStr;
            factorsLabels_this_except_last{end + 1} = betwStr;
            
            % If in recursive call for exploring interactions: Check whether the effect concerns factors of interest
            if ~isempty(factorLabelsOfInterest),
                if length(factorLabelsOfInterest) ~= length(factorsLabels_this),
                    continue;
                end;
                mismatch0 = 1;
                for iFactor2 = 1:length(factorsLabels_this),
                    if isempty(factorsLabels_this{iFactor2}),
                        ind_to_test = iFactor2;
                    else,
                        ind_to_test = 1:length(factorLabelsOfInterest);
                    end;
                    foundit = 0;
                    for iFactor = ind_to_test,
                        if strcmp(factorLabelsOfInterest{iFactor}, factorsLabels_this{iFactor2}),
                            foundit = 1;
                        end;
                    end;
                    mismatch0 = mismatch0 * foundit;
                end;
                if mismatch0 == 0,
                    continue;
                end;
            end;
            
            %
            % Run the desired statistical analysis on the current
            % within-subject, continuous, and between-group factors.
            %
            % Here: An implementation of Repeated Measures ANOVA.
            %
            [F, df1, df2, p, SSM, SSE, MSM, MSE, eta, eps0, bModel] = teg_RMA_ANOVA(M, X0, B, Bcoder, contvec, perm_test, nIts_perm);
            O.R = [O.R; F df1 df2 p MSM MSE];
            O.labels{end + 1} = [num2str(size(O.R, 1)) ': ' withinLabel betwStr contstr];
            
            %
            % If p is below the criterion for further processing:
            %
            
            if p <= pCritForFurther,
                
                Btmp = [];
                Btypes = [];
                if iBetwInt > 0,
                    [Btmp, Btypes] = teg_RMA_B_to_codes(B, bfactors{iBetwInt});
                end;
                
                if iEffect == 0,
                    M_reduced = mean(M, 2);
                    teg_RMA_report('', X0, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, M_reduced, 2, O.labels{end}, Btmp, Btypes, contvec);
                else,
                    % Get reduced-size dependent-variable matrix.
                    % This only has one column per combination of
                    % relevant levels.
                    M_reduced = teg_RMA_inner_recode_raw(M, NM, cellsets, iEffect);
                    
                    if length(factors{iEffect}) > 1,
                        teg_RMA_report('', X0, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, M_reduced, 0, O.labels{end}, Btmp, Btypes, contvec);
                        % Use M_reduced, factors, nColsFactor
                        lastFactor = factors{iEffect}(end);
                        nLevelsLastFactor = nColsFactor(lastFactor) + 1;
                        posthoc_levels = levels(factors{iEffect}(1:(end - 1)));
                        posthoc_varnames = {};
                        currModStr = 'factors:';
                        for iFactor = 1:(length(factors{iEffect}) - 1),
                            posthoc_varnames{end + 1} = varnames{factors{iEffect}(iFactor)};
                            currModStr = [currModStr ' ' posthoc_varnames{end}];
                        end;
                        currModStr = [currModStr ' ' factorsLabels_this{end - 1} ' ' factorsLabels_this{end}];
                        fprintf([all_str_offset '# START Post-hoc analyses start of model with ' currModStr '\n']);
                        for iLevel = 1:nLevelsLastFactor,
                            posthoc_M = M_reduced(:, iLevel:nLevelsLastFactor:end);
                            fprintf([all_str_offset 'Level ' num2str(iLevel) ' of factor ' varnames{factors{iEffect}(end)} ', model with ' currModStr '\n']);
                            if ~isempty(factorsLabels_this_except_last{end}),
                                breakpoint = 1;
                            end;
                            O_level_down = teg_RMA(posthoc_M, posthoc_levels, posthoc_varnames, B, Bvarnames, ContRaw, contRawLabels, factorsLabels_this_except_last, [argcode 'd']);
                        end;
                        fprintf([all_str_offset '# END Post-hoc analyses end of model with ' currModStr '\n\n']);
                    else,
                        % If not an interactions: Provide additional
                        % descriptives.
                        teg_RMA_report('', X0, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, M_reduced, 2, O.labels{end}, Btmp, Btypes, contvec);
                    end;
                end;
            end;
        end;
    end;
end;
