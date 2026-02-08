function plotSTL(filename)
TR = stlread(filename);

pltId = PlotId;
fig = figure(pltId.geom_stl); clf(fig, 'reset');
set(gcf, 'Name', 'STL Model');
trisurf(TR, 'FaceColor', 'cyan', 'EdgeColor','none','FaceAlpha',0.75);
axis equal; view(3); axis off
camlight; lighting gouraud; material metal;

end