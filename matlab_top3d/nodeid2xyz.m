function [x, y, z] = nodeid2xyz(node_id, nelx, nely)
    % forward formula: node_id = z*(nelx+1)*(nely+1)+x*(nely+1)+(nely+1-y);
    y_rem = rem(node_id, nely+1);
    y = mod(-y_rem, nely+1);
    x = rem((node_id - (nely+1-y))/(nely+1), nelx+1);
    z = (node_id-x*(nely+1)-(nely+1-y))/((nelx+1)*(nely+1));
end