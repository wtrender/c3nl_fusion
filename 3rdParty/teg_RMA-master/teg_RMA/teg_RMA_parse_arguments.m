function [M, NM, levels, varnames, B, Bvarnames, ContRaw, contRawLabels, factorLabelsOfInterest] = teg_RMA_parse_arguments(vararginPassed)

M = vararginPassed{1};
levels = vararginPassed{2};
% The number of observations per cell can be added as the second half of M.
% This will then be used to remove difference between results due to 
% splitting the per-subject data differently over factors.
if size(M, 2) > prod(levels), 
    b = size(M, 2) / 2;
    NM = M(:, (b + 1):end);
    M = M(:, 1:b);
else
    NM = ones(size(M));
end;

% Remove and store subject effects.
% Take different cell counts into account here if provided.
subjEffect = [];
for iSubj = 1:size(M, 1),
    fNonNaN = ~isnan(M(iSubj, :));
    subjEffect(iSubj, 1) = sum(M(iSubj, fNonNaN) .* NM(iSubj, fNonNaN), 2) ./ sum(NM(iSubj, fNonNaN));
end;
% M_subj = subjEffect * ones(1, size(M, 2));
% M_raw = M;
% M = M - M_subj;

varnames = vararginPassed{3};
if isempty(varnames),
    varnames = cellstr(num2str((1:length(levels))'));
end;

if length(vararginPassed) > 3,
    B = vararginPassed{4};
    Bvarnames = vararginPassed{5};
else,
    B = [];
    Bvarnames = {};
end;

if length(vararginPassed) > 5,
    ContRaw = vararginPassed{6};
    contRawLabels = vararginPassed{7};
else,
    ContRaw = [];
    contRawLabels = {};
end;

if length(vararginPassed) > 8, % There must be an argcode
    factorLabelsOfInterest = vararginPassed{8};
else,
    factorLabelsOfInterest = {};
end;
