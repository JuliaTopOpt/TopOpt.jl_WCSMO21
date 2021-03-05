% AN 169 LINE 3D TOPOLOGY OPITMIZATION CODE BY LIU AND TOVAR (JUL 2013)
% https://github.com/Top3dAPP/Top3d
% https://www.top3d.app/tutorials/
% https://link.springer.com/article/10.1007%252Fs00158-014-1107-x
% ------
% * TopOpt.jl comparison
% With modifications from Yijiang Huang (yijiangh@mit.edu)
% TODO
% - MMA, OC, fmincon
% - density & sensitivity filter

% MMA obj 2034
% OC obj 1473

% run_top3d(20, 10, 2, 0.3, 3.0, 2.0)
% run_top3d(6, 4, 2, 0.3, 3.0, 2.0)
nelx = 20;
nely = 10;
nelz = 2;
volfrac = 0.3;
penal = 3.0;
rmin = 2.0;
% 'OC', 'MMA', 'fmincon'
optimizer = 'OC'; 
% filter_type = 'density'; % 'sensitivity'
filter_type = 'sensitivity'; % ''
run_top3d(nelx,nely,nelz,volfrac,penal,rmin,optimizer,filter_type);

function run_top3d(nelx,nely,nelz,volfrac,penal,rmin,optimizer,filter_type)
close all
clc
% record time
tic

verbose = true; % printout iterations
fprintf('Optimizer: %s, Filter: %s\n', optimizer, filter_type)

if strcmp(optimizer, 'MMA')
    % load MMA code
    addpath('GGP-Matlab')
end

% USER-DEFINED LOOP PARAMETERS
maxloop = 400;    % Maximum number of iterations
tolx = 1e-6;      % Terminarion criterion
% tolx = 0.001;      % Terminarion criterion
displayflag = 0;  % Display intermediate structure flag

% USER-DEFINED MATERIAL PROPERTIES
E0 = 1.0;         % Young's modulus of solid material
% Young's modulus of void-like material
Emin = 1e-3; % 0.001
nu = 0.3;         % Poisson's ratio

% USER-DEFINED LOAD DOFs
il = nelx; jl = 0; kl = 0:nelz;                         % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
force = -1.0;

% USER-DEFINED SUPPORT FIXED DOFs
[jf,kf] = meshgrid(1:nely+1,1:nelz+1);                  % Coordinates
fixednid = (kf-1)*(nely+1)*(nelx+1)+jf;                 % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs

% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
% external load vector
F = sparse(loaddof,1,force,ndof,1);
% displacement vector
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
% element stiffness matrix
KE = lk_H8(nu);

% first grid of nodes in the x-y plane (z=0)
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));

% node IDs of the first node at each element
edofVec = 3*nodeids(:)+1;

% connectivity matrix, each row is node indices
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);

% * PREPARE DENSITY FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

% INITIALIZE ITERATION
x = repmat(volfrac,[nely,nelx,nelz]);
% x = ones(nely,nelx,nelz);
xPhys = x; 

loop = 0; 
change = 1;

if strcmp(optimizer, 'MMA')
    % * INITIALIZE MMA OPTIMIZER
    m     = 1;                % The number of general constraints.
    n     = nele;             % The number of design variables x_j.
    xmin  = zeros(n,1);       % Column vector with the lower bounds for the variables x_j.
    xmax  = ones(n,1);        % Column vector with the upper bounds for the variables x_j.
    xold1 = x(:);             % xval, one iteration ago (provided that iter>1).
    xold2 = x(:);             % xval, two iterations ago (provided that iter>2).
    low   = ones(n,1);        % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
    upp   = ones(n,1);        % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
    a0    = 1;                % The constants a_0 in the term a_0*z.
    a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
    c_MMA = 10000*ones(m,1);  % Column vector with the constants c_i in the terms c_i*y_i.
    d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
end

global ce % Shared between myfun and myHessianFcn, used by fmincon
if ~strcmp(optimizer, 'fmincon')
    % * START MAIN ITERATION
    while change > tolx && loop < maxloop
        loop = loop+1;
        % * FE-ANALYSIS
        sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),24*24*nele,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);

        % * OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        % element-wise compliance
        ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
        c = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)));
        dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
        dv = ones(nely,nelx,nelz);

        % * FILTERING AND MODIFICATION OF SENSITIVITIES
        if strcmp(filter_type, 'density')
            % compliance gradient
            dc(:) = H*(dc(:)./Hs);  
            % volume gradient
            dv(:) = H*(dv(:)./Hs);
        elseif strcmp(filter_type, 'sensitivity')
            dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
        end

        if strcmp(optimizer, 'OC')
            % * OPTIMALITY CRITERIA UPDATE
            l1 = 0; l2 = 1e9; move = 0.2;
            while (l2-l1)/(l1+l2) > 1e-3
                lmid = 0.5*(l2+l1);
                xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));

                if strcmp(filter_type, 'density')
                    xPhys(:) = (H*xnew(:))./Hs;
                elseif strcmp(filter_type, 'sensitivity')
                    xPhys = xnew;
                end

                if sum(xPhys(:)) > volfrac*nele, l1 = lmid; else l2 = lmid; end
            end
        elseif strcmp(optimizer, 'MMA')
            % * METHOD OF MOVING ASYMPTOTES
            xval  = x(:);
            f0val = c;
            df0dx = dc(:);
            fval  = sum(xPhys(:))/(volfrac*nele) - 1;
            dfdx  = dv(:)' / (volfrac*nele);
            [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low,upp] = ...
            mmasub(m, n, loop, xval, xmin, xmax, xold1, xold2, ...
            f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
            % * Update MMA Variables
            xnew     = reshape(xmma,nely,nelx,nelz);

            if strcmp(filter_type, 'density')
                xPhys(:) = (H*xnew(:))./Hs;
            elseif strcmp(filter_type, 'sensitivity')
                xPhys = xnew;
            end

            xold2    = xold1(:);
            xold1    = x(:);
        end

        change = max(abs(xnew(:)-x(:)));
        x = xnew;
        % PRINT RESULTS
        if verbose, fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f delta x.:%7.3f\n',loop,c,mean(xPhys(:)),change); end
        % PLOT DENSITIES
        if displayflag, clf; display_3D(xPhys); end %#ok<UNRCH>
    end
else
    % * fmincon
    assert(strcmp(filter_type, 'density'), 'Sensitivity filter not implemented yet for fmincon.');
    A = [];
    B = [];
    Aeq = [];
    Beq = [];
    LB = zeros(size(x));
    UB = ones(size(x));
    OPTIONS = optimset('TolX',tolx, 'MaxIter',maxloop, 'Algorithm','interior-point',...
    'GradObj','on', 'GradConstr','on', 'Hessian','user-supplied', 'HessFcn', @myHessianFcn,'Display','none', ...
    'OutputFcn', @(x,optimValues,state) myOutputFcn(x,optimValues,state,displayflag,verbose)...
    );
    % 'PlotFcns',@optimplotfval);
    fmincon(@myObjFcn, x, A, B, Aeq, Beq, LB, UB, @myConstrFcn, OPTIONS);
end

if change < tolx
    term_reason = 'Delta x Converged';
else
    term_reason = 'Iteration limits exceeded';
end
fprintf("Topopt total computation time: %0.2f s, terminated due to %s\n", toc, term_reason)
clf; 
display_3D(xPhys, loaddof(:), fixeddof(:), force);

% fmincon function definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, gradf] = myObjFcn(x)
    % TODO sensitivity filter
    xPhys(:) = (H*x(:))./Hs;
    % FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),24*24*nele,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    % FILTERING AND MODIFICATION OF SENSITIVITIES
    % TODO sensitivity filter
    dc(:) = H*(dc(:)./Hs);
    % RETURN
    f = c;
    gradf = dc(:);
end % myfun

function h = myHessianFcn(x, lambda)
    xPhys = reshape(x,nely,nelx,nelz);
    % Compute Hessian of Obj.
    Hessf = 2*(penal*(E0-Emin)*xPhys.^(penal-1)).^2 ./ (E0 + (E0-Emin)*xPhys.^penal) .* ce;
    % TODO sensitivity filter
    Hessf(:) = H*(Hessf(:)./Hs);
    % Compute Hessian of constraints
    Hessc = 0; % Linear constraint
    % Hessian of Lagrange
    h = diag(Hessf(:)) + lambda.ineqnonlin*Hessc;
end % myHessianFcn

function [cneq, ceq, gradc, gradceq] = myConstrFcn(x)
    % TODO sensitivity filter
    xPhys(:) = (H*x(:))./Hs;
    % Non-linear Constraints
    cneq  = sum(xPhys(:)) - volfrac*nele;
    gradc = ones(nele,1);
    % Linear Constraints
    ceq     = [];
    gradceq = [];
end % mycon

function stop = myOutputFcn(x,optimValues,state,displayflag,verbose)
    stop = false;
    switch state
        case 'iter'
            % Make updates to plot or guis as needed
            xPhys = reshape(x, nely, nelx, nelz);
            %% PRINT RESULTS
            if verbose
                fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',optimValues.iteration,optimValues.fval, ...
                    mean(xPhys(:)),optimValues.stepsize);
            end
            %% PLOT DENSITIES
            if displayflag, figure(10); clf; display_3D(xPhys); end
            title([' It.:',sprintf('%5i',optimValues.iteration),...
                ' Obj. = ',sprintf('%11.4f',optimValues.fval),...
                ' ch.:',sprintf('%7.3f',optimValues.stepsize)]);
        case 'init'
            % Setup for plots or guis
            if displayflag
                figure(10)
            end
        case 'done'
            fprintf('fmincon done.\n');
            % Cleanup of plots, guis, or final plot
            % fprintf("Topopt total computation time: %0.2f s", toc)
            % figure(10); clf; display_3D(xPhys);
        otherwise
    end % switch
end % myOutputFcn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % main function

% =========================================================================
% === This code was written by K Liu and A Tovar, Dept. of Mechanical   ===
% === Engineering, Indiana University-Purdue University Indianapolis,   ===
% === Indiana, United States of America                                 ===
% === ----------------------------------------------------------------- ===
% === Please send your suggestions and comments to: kailiu@iupui.edu    ===
% === ----------------------------------------------------------------- ===
% === The code is intended for educational purposes, and the details    ===
% === and extensions can be found in the paper:                         ===
% === K. Liu and A. Tovar, "An efficient 3D topology optimization code  ===
% === written in Matlab", Struct Multidisc Optim, 50(6): 1175-1196, 2014, =
% === doi:10.1007/s00158-014-1107-x                                     ===
% === ----------------------------------------------------------------- ===
% === The code as well as an uncorrected version of the paper can be    ===
% === downloaded from the website: http://www.top3dapp.com/             ===
% === ----------------------------------------------------------------- ===
% === Disclaimer:                                                       ===
% === The authors reserves all rights for the program.                  ===
% === The code may be distributed and used for educational purposes.    ===
% === The authors do not guarantee that the code is free from errors, and =
% === they shall not be liable in any event caused by the use of the code.=
% =========================================================================
