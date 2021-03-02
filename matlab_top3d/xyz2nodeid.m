function [node_id] = xyz2nodeid(x, y, z, nelx, nely)
    node_id = z*(nelx+1)*(nely+1)+x*(nely+1)+(nely+1-y);
end