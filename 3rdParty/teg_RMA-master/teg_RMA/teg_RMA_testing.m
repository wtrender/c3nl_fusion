function sim_results = teg_RMA_testing

nIts = 1;
nSubj = 320;
withinGroupEffect = 1;
betweenGroupFactor = 1;
betweenGroupEffect = 0;
continuousFactor = 0;
continuousEffect = 0;
noiseFactor = 0.005;

sim_results = [];
for iIt = 1:nIts,
    fprintf(['Sim ' num2str(iIt) '\n']);
    
    levels = [2, 3]; 
    labels = {'W1', 'W2'}; 
    M = randn(nSubj, prod(levels)); 

    % Within-subject effect of W1 and interaction W1 x W2
    if withinGroupEffect == 1,
        M = M + ones(size(M, 1), 1) * [-1 1 -1 1 -1 1];
    end;
    
    % create within-subject dependence
    Noise = randn(size(M));
    for iCol = 2:size(M, 2),
        for iCol2 = 1:size(M, 2),
            if iCol2 ~= iCol,
                Noise(:, iCol) = (Noise(:, iCol) + Noise(:, iCol2)) ./ 2;
            end;
        end;
    end;
    M = M + noiseFactor * Noise;

    % Between-group factors
    B = [];
    betwLabels = {};
    if betweenGroupFactor == 1,
        % Create hugely different group size
        B = [floor(10 * rand(nSubj, 2))];
        betwLabels = {'b1', 'b2'};
        bcol = B(:, 1);
        f1 = find(bcol == 0);
        f2 = find(bcol >1);
        if betweenGroupEffect == 1,
            M(f1, :) = -M(f1, :);
        end;
    end;
    
    % Continuous variable
    Cont = []; 
    contLabels = {};
    if continuousFactor == 1,
        Cont = randn(nSubj, 1); 
        Cont = teg_demean(Cont);
        contLabels = {'c1'};
        if continuousEffect == 1,
            M(:, 1:(size(M, 2) / levels(1))) = M(:, 1:(size(M, 2) / levels(1))) .* (Cont(:, 1) * ones(1, (size(M, 2) / levels(1))));
        end;
    end;

    O = teg_RMA(M, levels, labels, B, betwLabels, Cont, contLabels);
    p_values = O.R(:, 4);
    sim_results = [sim_results; p_values(:)'];
end;

fprintf(['Proportion of tests with p below 0.05:\n']);
disp(mean(sim_results < 0.05, 1));
