% === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho, loaddof, fixeddof, force)
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
                % ? y axis mirrored (for no reason?)
                % vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',...
                [0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
% * draw load vector
for s=1:size(loaddof,1)
    dir_id = mod(loaddof(s),3);
    if dir_id == 0
        dir_id = 3;
    end
    fix_vec = zeros(3,1);
    fix_vec(dir_id) = 1;
    fix_vec = fix_vec .* force;
    node_id = floor(loaddof(s)/3)+1;
    [xx,yy,zz] = nodeid2xyz(node_id, nelx, nely);
    xcoord = [xx,yy,zz] .* [hx, hy, hz];
    quiver3(xcoord(1), xcoord(3), xcoord(2), fix_vec(1),fix_vec(3),fix_vec(2), 'r','LineWidth', 1);
    hold on;
end
% * draw fixities vector
for s=1:size(fixeddof,1)
    dir_id = mod(fixeddof(s),3);
    if dir_id == 0
        dir_id = 3;
    end
    fix_vec = zeros(3,1);
    fix_vec(dir_id) = 1;
    node_id = floor(fixeddof(s)/3)+1;
    [xx, yy, zz] = nodeid2xyz(node_id, nelx, nely);
    xcoord = [xx, yy, zz] .* [hx, hy, hz];
    quiver3(xcoord(1), xcoord(3), xcoord(2), fix_vec(1),fix_vec(3),fix_vec(2),'b','LineWidth', 1);
    hold on;    
end
xlabel("x");ylabel("z");zlabel("y");
axis equal; axis tight; axis on; box on; view([30,30]); pause(1e-6);
end