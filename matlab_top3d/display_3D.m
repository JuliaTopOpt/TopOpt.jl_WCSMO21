% === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho, loadnid, fixednid)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1; % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5)  % User-defined display density threshold
                vert = [x y z;    x y-hx z;    x+hx y-hx z;    x+hx y z; ...
                        x y z+hx; x y-hx z+hx; x+hx y-hx z+hx; x+hx y z+hx];
                % ? y-z axis swapped
                vert(:,[2 3]) = vert(:,[3 2]); 
                vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',...
                [0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
% draw load vector
for s=1:size(loadnid,1)
%     node_id = fix(load_dofs(s)/3)+1;
    node_id = loadnid(s);
    [i,j,k] = nodeid2xyz(node_id, nelx, nely, 0);
    x = (i+0)*hx;
    y = (j+0)*hy;
    z = (k+0)*hz;
    % hardcode load vector for now
    fprintf("Load | nodeid %d, (%d, %d, %d), x: %f, y: %f, z: %f\n", node_id, i,j,k, x, y, z)
    quiver3(x,y,z,0,-1.0,0,'r','LineWidth',1);
    hold on;
end
% draw fixities vector
for s=1:size(fixednid,1)
%     dir_id = rem(fixeddof(s),3);
%     if dir_id == 0
%         dir_id = 3;
%     end
%     fix_vec = zeros(3,1);
%     fix_vec(dir_id) = 1;
%     node_id = fix(fixeddof(s)/3)+1;
    node_id = fixednid(s);
    [i,j,k] = nodeid2xyz(node_id, nelx, nely, 0);
    x = (i)*hx;
    y = (j)*hy;
    z = (k-1)*hz;
    if x > 0
        fprintf("Fix | nodeid %d, (%d, %d, %d), x: %f, y: %f, z: %f\n", node_id, i,j,k, x, y, z)
    end
%     quiver3(x,y,z,fix_vec(1),fix_vec(2),fix_vec(3),'b','LineWidth',4);
    vec = eye(3);
    for t=1:3
        quiver3(x,y,z,vec(t,1),vec(t,2),vec(t,3),'b','LineWidth',4);
        hold on;    
    end
end
xlabel("x");ylabel("z");zlabel("y");
axis equal; axis tight; axis on; box on; view([30,30]); pause(1e-6);
end