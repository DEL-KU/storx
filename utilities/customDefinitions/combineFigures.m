%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by:                                          %
% - Amir M. Mirzendehdel, Aerospace Engineering Department, KU              %
% - Krishnan Suresh, Mechanical Engineering Department, UW-Madison          %
%                                                                           %
% Please send your comments to: amirzend@ku.edu                             %
%                                                                           %
% The code is intended for educational purposes and theoretical details     %
% are discussed in the textbook:                                            %
% Introduction to Shape and Topology Optimization using MATLAB              %
%                                                                           %
% Disclaimer:                                                               %
% The authors reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, we shall not be liable in any event        %
% caused by the use of the program.                                         %
%                                                                           %
% License:                                                                  %
% This software is used, copied and distributed under the licensing         %
% agreement contained in the file LICENSE in the top directory of           %
% the distribution.                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combineFigures(fig_title)
disp('Combining figures, this might take a few seconds ...');
% Get a list of all open figures
pltId = PlotId;
figs = findall(0, 'Type', 'figure');

% Get figure numbers
figNumbers = arrayfun(@(x) x.Number, figs);

% Sort figure numbers
[~, sortIndex] = sort(figNumbers,'descend');

% Sort figure handles based on sorted figure numbers
figs = figs(sortIndex);

numFigs = length(figs);

if numFigs < 1
    disp('No figures to combine');
    return;
end

% Create a new figure for the combined subplots
newFig = figure;
set(newFig, 'Name', fig_title);

% Maximize the size of the figure (set to full screen size)
set(newFig, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Set the number of columns to 3 and calculate the required number of rows
numCols = 3;
numRows = ceil(numFigs / numCols);
tiledlayout(numRows, numCols, 'TileSpacing', 'compact', 'Padding', 'compact');

% Loop over each figure in reverse order (oldest first)
for i = 1:numFigs
    id = numFigs - i + 1;
    newAxes = nexttile;

    % Get the axes of the original figure
    oldAxes = findall(figs(id), 'Type', 'axes');
    % Set the figure name as the title of the subplot
    figName = get(figs(id), 'Name');
    if isempty(figName)
        figName = ['Figure ', num2str(figs(id).Number)];
    end

    if (figs(id).Number == pltId.convergence)
        plotConvergence(figs(id),newAxes,figName);
        continue;
    end

    % Copy the children (like plot lines, images, etc.) to the new subplot
    for j = 1:length(oldAxes)

        oldChildren = get(oldAxes(j), 'Children');

        % Copy only the visible children (i.e., the actual plot objects) to the new axes
        visibleChildren = findobj(oldChildren, 'Visible', 'on');
        copyobj(visibleChildren, newAxes);

        if (figs(id).Number == pltId.design)
            % Assume you have a figure with an image and you've copied the image object handle
            h = findobj(gca, 'Type', 'image'); % Find the image object in the current axes

            % Flip the image data vertically
            h.CData = flipud(h.CData);
        end

        % Preserve axes limits, labels, titles, etc.
        set(newAxes, 'XLim', get(oldAxes(j), 'XLim'));
        set(newAxes, 'YLim', get(oldAxes(j), 'YLim'));
        set(newAxes, 'ZLim', get(oldAxes(j), 'ZLim'));

        zticklabels = get(oldAxes(j),'ZTickLabel');
        set(newAxes,'ZTickLabel',zticklabels)

        % Copy aspect ratio
        set(newAxes, 'PlotBoxAspectRatio', get(oldAxes(j), 'PlotBoxAspectRatio'));

        % Preserve colormap
        colormap(newAxes, colormap(figs(id)));

        % Copy labels and titles
        xlabel(newAxes, get(oldAxes(j), 'XLabel').String);
        ylabel(newAxes, get(oldAxes(j), 'YLabel').String);
        zlabel(newAxes, get(oldAxes(j), 'ZLabel').String);

        % Copy colorbar if it exists
        oldColorbar = findall(figs(id), 'Type', 'colorbar');
        if ~isempty(oldColorbar)
            colorbar(newAxes);  % Create colorbar in the new axes
        end

        % Copy color limits (clim)
        caxis(newAxes, caxis(oldAxes(j)));  % Copy color limits

        % Copy the same view (for 2D or 3D view)
        view(newAxes, get(oldAxes(j), 'View'));

        % Copy legend
        plbc = PlotBC;
        legend_fields = {};
        legend_labels = {};
        hLegend = findobj(figs(id), 'Type', 'legend'); % Find the legend in the original figure
        if ~isempty(hLegend)
            hold on;
            labels = hLegend.String;
            for k = 1:numel(labels)
                label = labels{k};
                if strcmp(label,'fixed $T$')
                    fixed_T = plot(NaN,NaN, ...
                        plbc.fixed_T.marker,'MarkerEdgeColor', ...
                        plbc.fixed_T.color,'MarkerSize',10); hold on
                    % legend
                    legend_fields = [legend_fields;fixed_T ]; %#ok
                    legend_labels = [legend_labels,'fixed $T$']; %#ok
                elseif strcmp(label,'heat flux')
                    flux = plot(NaN, NaN,plbc.flux.marker, ...
                        'MarkerEdgeColor',plbc.flux.color, ...
                        'MarkerFaceColor', plbc.flux.color); hold on
                    legend_fields = [legend_fields;flux ]; %#ok
                    legend_labels = [legend_labels,'heat flux']; %#ok

                elseif strcmp(label,'internal heat')
                    internal_heat = plot(NaN, NaN,plbc.internal_heat.marker, ...
                        'MarkerEdgeColor',plbc.internal_heat.color, ...
                        'MarkerFaceColor', plbc.internal_heat.color); hold on
                    legend_fields = [legend_fields;internal_heat ]; %#ok
                    legend_labels = [legend_labels,'internal heat']; %#ok
                elseif strcmp(label,'fixed $u$')
                    fixed_U = plot(NaN, NaN,plbc.fixed_U.marker, ...
                        'MarkerEdgeColor', plbc.fixed_U.color);
                    hold on;
                    % legend
                    legend_fields = [legend_fields;fixed_U ]; %#ok
                    legend_labels = [legend_labels,'fixed $u$']; %#ok
                elseif strcmp(label,'fixed $v$')
                    fixed_V = plot(NaN, NaN,plbc.fixed_V.marker, ...
                        'MarkerEdgeColor', plbc.fixed_V.color);
                    hold on;
                    % legend
                    legend_fields = [legend_fields;fixed_V ]; %#ok
                    legend_labels = [legend_labels,'fixed $v$']; %#ok
                elseif strcmp(label,'force')
                    force = plot(NaN, NaN,plbc.force.marker, ...
                        'MarkerEdgeColor',plbc.force.color, ...
                        'MarkerFaceColor', plbc.force.color);
                    hold on;
                    legend_fields = [legend_fields;force ]; %#ok
                    legend_labels = [legend_labels,'force']; %#ok
                elseif strcmp(label,'body force')
                    acceleration = plot(NaN, NaN,plbc.acceleration.marker, ...
                        'MarkerEdgeColor',plbc.acceleration.color, ...
                        'MarkerFaceColor', plbc.acceleration.color);
                    hold on;
                    legend_fields = [legend_fields;acceleration ]; %#ok
                    legend_labels = [legend_labels,'body force']; %#ok
                elseif strcmp(label,'flow $u$')
                    flow_U = plot(NaN, NaN, plbc.flow_U.marker, ...
                        'MarkerEdgeColor', plbc.flow_U.color, ...
                        'MarkerFaceColor', plbc.flow_U.color);
                    hold on;
                    legend_fields = [legend_fields;flow_U ]; %#ok
                    legend_labels = [legend_labels,'flow $u$']; %#ok
                elseif strcmp(label,'no-slip $u$')
                    noSlip_U = plot(NaN, NaN, ...
                        plbc.noSlip_U.marker, ...
                        'MarkerEdgeColor', plbc.noSlip_U.color);
                    hold on;
                    legend_fields = [legend_fields;noSlip_U ]; %#ok
                    legend_labels = [legend_labels,'no-slip $u$']; %#ok
                elseif strcmp(label,'flow $v$')
                    flow_V = plot(NaN, NaN, ...
                        plbc.flow_V.marker, ...
                        'MarkerEdgeColor',plbc.flow_V.color, ...
                        'MarkerFaceColor', plbc.flow_V.color);
                    hold on;
                    legend_fields = [legend_fields;flow_V ]; %#ok
                    legend_labels = [legend_labels,'flow $v$']; %#ok
                elseif strcmp(label,'no-slip $v$')
                    noSlip_V = plot(NaN, NaN, plbc.noSlip_V.marker, ...
                        'MarkerEdgeColor', plbc.noSlip_V.color);
                    hold on;
                    legend_fields = [legend_fields;noSlip_V ]; %#ok
                    legend_labels = [legend_labels,'no-slip $v$']; %#ok
                elseif strcmp(label,'fixed $p$')
                    fixed_P = plot(NaN, NaN,plbc.fixed_P.marker, ...
                        'MarkerFaceColor',plbc.fixed_P.color);
                    hold on;
                    legend_fields = [legend_fields;fixed_P ]; %#ok
                    legend_labels = [legend_labels,'fixed $p$']; %#ok
                elseif strcmp(label,'active design domain')
                    activeDomain = plot(NaN, NaN,plbc.activeDomain.marker, ...
                        'MarkerFaceColor',plbc.activeDomain.color);
                    hold on;
                    legend_fields = [legend_fields;activeDomain ]; %#ok
                    legend_labels = [legend_labels,'active design domain']; %#ok
                end
            end
            if ~isempty(legend_labels)
                legend(legend_fields,legend_labels, ...
                    'Location', hLegend.Location);
            end
        end
    end

    title(newAxes, figName);
end

% Close all original figures after copying them
close(figs);
disp('Combining figures successfully completed!');
end

function plotConvergence(figHandle,newAxes,figName)
originalAxes = findall(figHandle, 'type', 'axes');
children = allchild(originalAxes); % Get all child objects of the axes
plot(children(2).XData,children(2).YData, ...
    'LineStyle', children(2).LineStyle,  ...
    'Color', children(2).Color,  ...
    'LineWidth',children(2).LineWidth);
ylabel(newAxes, get(get(originalAxes.YAxis(1), 'Label'), 'String'));

yyaxis right
plot(children(1).XData,children(1).YData, ...
    'LineStyle', children(1).LineStyle,  ...
    'Color', children(1).Color,  ...
    'LineWidth',children(1).LineWidth);

ylabel(newAxes, get(get(originalAxes.YAxis(2), 'Label'), 'String'));

% Preserve axes limits, labels, titles, etc.
set(newAxes, 'XLim', get(originalAxes, 'XLim'));
set(newAxes, 'YLim', get(originalAxes, 'YLim'));
set(newAxes, 'ZLim', get(originalAxes, 'ZLim'));

% Copy aspect ratio
set(newAxes, 'PlotBoxAspectRatio', get(originalAxes, 'PlotBoxAspectRatio'));
title(newAxes, figName);
end

