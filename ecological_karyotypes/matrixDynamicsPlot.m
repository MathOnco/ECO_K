function [] = matrixDynamicsPlot(matrix, cloneIDs,origin)
 
% Define the output directory where the .mat files are saved
OUTDIR = 'Results/'; % Adjust this path if necessary

% Remove 'Origin-' prefix from origin to match the filenames
originWithoutPrefix = strrep(origin, 'Origin-', '');

% Construct the .mat filename based on the origin
matFileName = ['clone_colors_', originWithoutPrefix, '.mat'];

% Check if the .mat file exists
if exist(matFileName, 'file')
    % Load the cloneColorTable variable from the .mat file
    load(matFileName, 'cloneColorTable');
    % Filter by clones in this Replicate Group
    idx = ismember(cloneColorTable.CloneID, cloneIDs);
    cloneColorTable = cloneColorTable(idx, :);
    % Extract the RGB color values into nodecols
    nodecols = [cloneColorTable.R, cloneColorTable.G, cloneColorTable.B];
else
    error(['Color table file not found for origin ', originWithoutPrefix]);
end

edgecols={"#8C000F"	,"#054907"};
matrix=matrix';
G = digraph(matrix);
h = plot(G,'LineWidth',0.1,'EdgeColor','white', 'NodeLabel', cloneIDs, 'NodeColor', nodecols,'MarkerSize',16, 'NodeFontSize', 16, 'ArrowSize', 16, 'EdgeAlpha', 0.75);
for x = matrix(matrix~=0)'
    A_=matrix;
    A_(A_~=x) = 0;
    G_sub = digraph(A_);
    highlight(h,G_sub,'EdgeColor',edgecols{1+(x>0)},'LineWidth',10*abs(x))
end

end