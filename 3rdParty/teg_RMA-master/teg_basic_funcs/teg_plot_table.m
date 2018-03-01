function teg_plot_table(varargin)

% function teg_plot_table(m, se0[, minyzero, betweenclusterlabels, withinclusterlabels, title0])
%
% Nesting is columns in rows.
% Rows are plotted on the main X axis.
% Columns are plotted within the clusters.
%
% That is: the values of each row will be grouped together as a cluster of bars.
%
% Example:
% teg_plot_table([1 2 3 2; 3 2 1 -2; 7 8 9 9], [1 2 1 1; 3 4 4 1; 1 1 1 2], 1, {'A', 'B', 'C'}, {'1', '2', '3', '4'}, 'Test')
%
% Adjust the second element of set(gca, 'Position', [0.1300 0.200 0.7750
% 0.7150]) if the labels go off-figure

m = varargin{1};
if length(varargin) > 1,
    se0 = varargin{2};
else,
    se0 = [];
end;
if length(varargin) > 2,
    minyzero = varargin{3};
else,
    minyzero = 1;
end;
if length(varargin) > 3,
    betwlab = varargin{4};
else,
    betwlab = {};
end;
if length(varargin) > 4,
    withinlab = varargin{5};
else,
    withinlab = {};
end;
if length(varargin) > 5,
    title0 = varargin{6};
else,
    title0 = [];
end;
if size(m, 1) == 1,
    m = m';
    se0 = se0';
end;

barSpacing = 0.66;
barW = 0.66 * barSpacing * 1 / size(m, 2);
lineW = barW;

xminall = Inf;
xmaxall = -Inf;
yminall = Inf;
ymaxall = -Inf;

midClusterX = [];
barXv = [];
barLab = {};
allyv = [];

for iR = 1:size(m, 1),
    clusterBaseX = iR;
    
    clusterXv = [];
    for iC = 1:size(m, 2),
        barBaseX = clusterBaseX + barSpacing * iC / size(m, 2);
        h = m(iR, iC);
        % Bar
        xv = [barBaseX barBaseX barBaseX + barW barBaseX + barW];
        yv = [0 h h 0];
        fill0 = fill(xv, yv, 'k');
        hold on;
        allyv = [allyv; h];
        % SE lines
        lX = barBaseX + barW / 2;
        if ~isempty(se0),
            hl = sign(h) * se0(iR, iC);
            xv = [lX lX];
            yv = [h h + hl];
            l0 = line(xv, yv); set(l0, 'Color', 'k');
            xv = [lX - lineW/2, lX + lineW/2];
            yv = [h + hl h + hl];
            l0 = line(xv, yv); set(l0, 'Color', 'k');
            allyv = [allyv; h + hl];
        end;
        clusterXv = [clusterXv lX];
        barXv = [barXv lX];
        if ~isempty(withinlab),
            thislab = withinlab{iC};
            if ~isempty(betwlab),
                thislab = [betwlab{iR} thislab];
            end;
            barLab{end + 1} = thislab;
        end;
    end;
    
    midClusterX(iR) = mean(clusterXv);

end;

xl0 = xlim;
miny = min(allyv);
maxy = max(allyv);
buffy = (maxy - miny) * 0.05;
maxy = maxy + buffy;
miny = miny - buffy;
if minyzero == 1,
    miny = min([0; miny]);
end;
ylim([miny, maxy]);

set(gca, 'XTick', midClusterX);
if ~isempty(betwlab) && isempty(withinlab),
    set(gca, 'XTickLabel', betwlab);
    set(gca, 'XTickLabelRotation', 45);
end;
if ~isempty(withinlab),
    set(gca, 'XTick', []);
    for iL = 1:length(barXv),
        t0 = text(barXv(iL), miny, barLab{iL});
        set(t0, 'Rotation', 45);
        ext0 = get(t0, 'Extent');
        pos0 = get(t0, 'Position');
        pos0(1) = pos0(1) - ext0(3) + barW / 2;
        pos0(2) = pos0(2) - ext0(4) * 1.2;
        set(t0, 'Position', pos0);
    end;
    set(gca, 'Position', [0.1300 0.200 0.7750 0.7150]);
end;

if ~isempty(title0),
    title(title0);
end;
set(gcf,'PaperPositionMode','auto');

fprintf(['print -dtiff -r600 plotfilename\n']);
