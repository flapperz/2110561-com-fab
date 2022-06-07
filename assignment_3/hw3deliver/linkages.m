function linkages(scene)
if nargin < 1
    scene = 0;
end

% Set up the scene here.
% Note that links don't have to be placed exactly. The first call to
% solveLinkage() will automatically snap the links so that all the
% constraints are satisfied.
linkages = [];
pins = [];
particles = [];
switch scene
    case 0
        % Crank-rocker
        range = [0,Inf];
        % Bottom link
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [-1;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]'; % display vertices
        % Left link
        links(2).angle = pi/2;
        links(2).pos = [-1;0];
        links(2).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]';
        % Right link
        links(3).angle = pi/2;
        links(3).pos = [1;0];
        links(3).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]';
        % Top link
        links(4).angle = 0;
        links(4).pos = [-1;1];
        links(4).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        
        % Which link is grounded?
        grounded = 1;
        % Which link is the driver?
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        % Bottom-left
        pins(1).links = [1,2];
        pins(1).pts = [0,0;0,0]';
        % Bottom-right
        pins(2).links = [1,3];
        pins(2).pts = [2,0;0,0]';
        % Left-top
        pins(3).links = [2,4];
        pins(3).pts = [1,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [1+rand(1),0;2,0]'; % pin location on link3 is randomized
        
        % List of tracer particles for display
        particles(1).link = 4; % which link?
        particles(1).pt = [0.5;0.1]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        particles(2).link = 4;
        particles(2).pt = [2.5;-0.1];
        particles(2).ptsWorld = zeros(2,0);
    case 1
        % Drag-link
        range = [0,Inf];
        % Bottom link
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [-1;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]'; % display vertices
        % Left link
        links(2).angle = pi/2;
        links(2).pos = [-1;0];
        links(2).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        % Right link
        links(3).angle = pi/2;
        links(3).pos = [1;0];
        links(3).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        % Top link
        links(4).angle = 0;
        links(4).pos = [-1;2];
        links(4).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        
        % Which link is grounded?
        grounded = 1;
        % Which link is the driver?
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        rn = rand(1)/3;
        
        % Bottom-left
        pins(1).links = [1,2];
        pins(1).pts = [0,0;0,0]';
        % Bottom-right
        pins(2).links = [1,3];
        pins(2).pts = [1.5,0;0,0]';
        % Left-top
        pins(3).links = [2,4];
        pins(3).pts = [2+rn,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [2+2*rn,0;2+2.5*rn,0]'; % pin location on link3 is randomized
        
        % List of tracer particles for display
        particles(1).link = 4; % which link?
        particles(1).pt = [0;0]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        particles(2).link = 4;
        particles(2).pt = [2+2*rn;0];
        particles(2).ptsWorld = zeros(2,0);
    case 2
        % Double-rocker
        range = [pi/2,2*pi/6];
        % Bottom link
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [-1;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]'; % display vertices
        % Left link
        links(2).angle = pi/2;
        links(2).pos = [-1;0];
        links(2).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        % Right link
        links(3).angle = pi/2;
        links(3).pos = [1;0];
        links(3).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        % Top link
        links(4).angle = 0;
        links(4).pos = [-1;2];
        links(4).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]';
        
        % Which link is grounded?
        grounded = 1;
        % Which link is the driver?
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        rn = rand(1)/3;
        
        % Bottom-left
        pins(1).links = [1,2];
        pins(1).pts = [0,0;0,0]';
        % Bottom-right
        pins(2).links = [1,3];
        pins(2).pts = [2,0;0,0]';
        % Left-top
        pins(3).links = [2,4];
        pins(3).pts = [2+rn,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [2+2*rn,0;1,0]'; % pin location on link3 is randomized
        
        % List of tracer particles for display
        particles(1).link = 4; % which link?
        particles(1).pt = [0;0]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        particles(2).link = 4;
        particles(2).pt = [1;0];
        particles(2).ptsWorld = zeros(2,0);
    case 3
        % Hoekens
        range = [0, Inf];
        % Bottom link
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [-1;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]'; % display vertices
        % Left link
        links(2).angle = pi/2;
        links(2).pos = [-1;0];
        links(2).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]';
        % Right link
        links(3).angle = pi/2;
        links(3).pos = [1;0];
        links(3).verts = [0,-0.1;2.5,-0.1;2.5,0.1;0,0.1]';
        % Top link
        links(4).angle = 0;
        links(4).pos = [-1;2];
        links(4).verts = [0,-0.1;5,-0.1;5,0.1;0,0.1]';
        
        % Which link is grounded?
        grounded = 1;
        % Which link is the driver?
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        rn = rand(1)/3;
        
        % Bottom-left
        pins(1).links = [1,2];
        pins(1).pts = [0,0;0,0]';
        % Bottom-right
        pins(2).links = [1,3];
        pins(2).pts = [2,0;0,0]';
        % Left-top
        pins(3).links = [2,4];
        pins(3).pts = [1,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [2.5,0;2.5,0]'; % pin location on link3 is randomized
        
        % List of tracer particles for display
        particles(1).link = 4; % which link?
        particles(1).pt = [0;0]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        particles(2).link = 4;
        particles(2).pt = [5;0];
        particles(2).ptsWorld = zeros(2,0);
    case 4
        % Peaucellier-Lipkin
        range = [6*pi/8, 2*pi/8];
        
        % ground link
        links(1).angle = pi/2; % rotation from the positive x-axis
        links(1).pos = [0;-1]; % position of the center of rotation
        links(1).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]'; % display vertices
        
        % driver link
        links(2).angle = pi/2; % rotation from the positive x-axis
        links(2).pos = [0;0]; % position of the center of rotation
        links(2).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]'; % display vertices
        
        % A
        links(3).angle = pi; % rotation from the positive x-axis
        links(3).pos = [0;1]; % position of the center of rotation
        links(3).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]'; % display vertices
         
        % B
        links(4).angle = pi/2; % rotation from the positive x-axis
        links(4).pos = [1;1]; % position of the center of rotation
        links(4).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]'; % display vertices
        
        % C
        links(5).angle = 0; % rotation from the positive x-axis
        links(5).pos = [-1;1]; % position of the center of rotation
        links(5).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]'; % display vertices
        
        % D
        links(6).angle = -pi; % rotation from the positive x-axis
        links(6).pos = [1;1]; % position of the center of rotation
        links(6).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]'; % display vertices
        
        % left
        links(7).angle = pi/2; % rotation from the positive x-axis
        links(7).pos = [0;-1]; % position of the center of rotation
        links(7).verts = [0,-0.1;2.5,-0.1;2.5,0.1;0,0.1]'; % display vertices
        
        % right
        links(8).angle = pi/2; % rotation from the positive x-axis
        links(8).pos = [0;-1]; % position of the center of rotation
        links(8).verts = [0,-0.1;2.5,-0.1;2.5,0.1;0,0.1]'; % display vertices
        
        % ground driver
        pins(1).links = [1,2];
        pins(1).pts = [1,0;0,0]';
        
        % driver-A
        pins(2).links = [2,3];
        pins(2).pts = [1,0;0,0]';
        
        % A-B
        pins(3).links = [3,4];
        pins(3).pts = [1,0;0,0]';
        
        % B-C
        pins(4).links = [4,5];
        pins(4).pts = [1,0;0,0]';
        
        % C-D
        pins(5).links = [5,6];
        pins(5).pts = [1,0;0,0]';
        
        % D-A
        pins(6).links = [6,3];
        pins(6).pts = [1,0;0,0]';
        
        % left-B base
        pins(7).links = [7,4];
        pins(7).pts = [2.5,0;0,0]';
        
        % right-D base
        pins(8).links = [8,6];
        pins(8).pts = [2.5,0;0,0]';
        
        % left-ground
        pins(9).links = [7,1];
        pins(9).pts = [0,0;0,0]';
        
        % right-ground
        pins(10).links = [8,1];
        pins(10).pts = [0,0;0,0]';
        
        grounded = 1;
        driver = 2;
        
        % List of tracer particles for display
        particles(1).link = 4; % which link?
        particles(1).pt = [1;0]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        
    case 5
        % Klann
        range = [0,Inf];
        
        grounded = 1;
        driver = 2;    
        
        rockerArmAxle1 = [1.366, 1.366];
        rockerArmAxle2 = [1.009, 0.574];
        crankShaft = [1.599,0.75];
        
        ldriver = 0.268;
        l6 = 0.5176;
        l4 = 0.522;
        l5 = 0.3206;
        l3 = 0.59;
        
        dh = 0.04;
        
        % Ground link
        links(grounded).angle = 0;
        links(grounded).pos = rockerArmAxle2';
        links(grounded).verts = [rockerArmAxle1-rockerArmAxle2; 0,0; crankShaft-rockerArmAxle2]';
        
        % Driver link
        links(driver).angle = 0;
        links(driver).pos = crankShaft';
        links(driver).verts = [0,0; ldriver/2,dh; ldriver,0; ldriver/2,-dh]';
        
        % l3
        links(3).angle = pi;
        links(3).pos = (crankShaft+[ldriver,0])';
        links(3).verts = [0,0; l3,0; 1.099,-0.116]';
        
        % l4 // mistake
        links(4).angle = 0;
        links(4).pos = ([0,0])';
        links(4).verts = [0,0;0,0]';
        
        % l5
        links(5).angle = -2*pi/3;
        links(5).pos = (crankShaft-[l3-ldriver,0])';
        links(5).verts = [0,0;l5/2,dh;l5,0;l5/2,-dh]';
        
        % l6
        links(6).angle = pi/2;
        links(6).pos = rockerArmAxle1';
        links(6).verts = [0,0;l6/2,dh;l6,0;l6/2,-dh]';
        
        % leg
        links(7).angle = -pi;
        links(7).pos = [1,1.732]';
        links(7).verts = [0,0;0.232,0.865;0,1.732]';
        
        % pin grounded-driver
        pins(1).links = [grounded,driver];
        pins(1).pts = [crankShaft-rockerArmAxle2;0,0]';
        
        % pin driver-l3
        pins(2).links = [driver,3];
        pins(2).pts = [ldriver,0;0,0]';
        
        % pin l3-l4
        pins(3).links = [3,4];
        pins(3).pts = [l3,0;0,0]';
        
        % pin l3-l5
        pins(4).links = [3,5];
        pins(4).pts = [l3,0;0,0]';
        
        % pin grounded-l5
        pins(5).links = [grounded,5];
        pins(5).pts = [0,0;l5,0]';
        
        % pin grounded-l6
        pins(6).links = [grounded,6];
        pins(6).pts = [rockerArmAxle1-rockerArmAxle2;0,0]';
        
        % pin leg-l3
        pins(7).links = [7,3];
        pins(7).pts = [0.232,0.866;1.099,-0.116]';
        
        % pin leg-l6
        pins(8).links = [7,6];
        pins(8).pts = [0,0;l6,0]';
        
        % List of tracer particles for display
        particles(1).link = 7; % which link?
        particles(1).pt = [0;1.732]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
    case 10
        % Extra credit!
end

% Initialize
for i = 1 : length(links)
    links(i).grounded = (i == grounded);
    links(i).driver = (i == driver);
    % These target quantities are only used for grounded and driver links
    links(i).angleTarget = links(i).angle;
    links(i).posTarget = links(i).pos;
end
drawScene(links,pins,particles,true);

% lsqnonlin options
if verLessThan('matlab','8.1')
    opt = optimset(...
        'Jacobian','on',...
        'DerivativeCheck','off',...
        'Display','off'); % final-detailed iter-detailed off
else
    opt = optimoptions(@lsqnonlin,...
        'Jacobian','on',...
        'DerivativeCheck','off',...
        'Display','final-detailed'); % final-detailed iter-detailed off
end

% Simulation loop
t = 0; % current time
T = 1.2; % final time
dt = 0.01; % time step
constAngVel = 2*pi; % driver angular velocity
angVel = 0;
angRange = range;
k = 10;
xMid = (angRange(1) + angRange(2))/2;
links(driver).angleTarget = angRange(1);
while t < T
    % Procedurally set the driver angle.
    % Right now, the target angle is being linearly increased, but you may
    
    xNow = links(driver).angleTarget;
    if angRange(2) == Inf
        links(driver).angleTarget = xNow + dt*constAngVel;
    else
        angAcc = -k * (xNow - xMid);
        angVel = angVel + dt * angAcc;
        % want to do something else.
        links(driver).angleTarget = xNow + dt*angVel;
    end
    % Solve for linkage orientations and positions
    [links,feasible] = solveLinkage(links,pins,opt);
    % Update particle positions
    particles = updateParticles(links,particles);
    % Draw scene
    drawScene(links,pins,particles);
    % Quit if over-constrained
    if ~feasible
        break;
    end
    t = t + dt;
end

end

%%
function [R,dR] = rotationMatrix(angle)
c = cos(angle);
s = sin(angle);
% Rotation matrix
R = zeros(2);
R(1,1) = c;
R(1,2) = -s;
R(2,1) = s;
R(2,2) = c;
if nargout >= 2
    % Rotation matrix derivative
    dR = zeros(2);
    dR(1,1) = -s;
    dR(1,2) = -c;
    dR(2,1) = c;
    dR(2,2) = -s;
end
end

%%
function [links,feasible] = solveLinkage(links,pins,opt)
nlinks = length(links);
% Extract the current angles and positions into a vector
angPos0 = zeros(3*nlinks,1);
for i = 1 : nlinks
    link = links(i);
    ii = (i-1)*3+(1:3);
    angPos0(ii(1)) = link.angle;
    angPos0(ii(2:3)) = link.pos;
end
% Limits
lb = -inf(size(angPos0));
ub =  inf(size(angPos0));
% Solve for angles and positions
[angPos,r2] = lsqnonlin(@(angPos)objFun(angPos,links,pins),angPos0,lb,ub,opt);
% If the mechanism is feasible, then the residual should be zero
feasible = true;
if r2 > 1e-6
    fprintf('Mechanism is over constrained!\n');
    feasible = false;
end
% Extract the angles and positions from the values in the vector
for i = 1 : length(links)
    ii = (i-1)*3+(1:3);
    links(i).angle = angPos(ii(1));
    links(i).pos = angPos(ii(2:3));
end
end

%%
function [f,J] = objFun(angPos,links,pins)
nlinks = length(links);
npins = length(pins);
% Temporarily change angles and positions of the links. These changes will
% be undone when exiting this function.
for i = 1 : nlinks
    ii = (i-1)*3+(1:3);
    links(i).angle = angPos(ii(1));
    links(i).pos = angPos(ii(2:3));
end

% Evaluate constraints
ndof = 3*nlinks;
ncon = 3 + 3 + 2*npins; % 3 for ground, 3 for driver, 2*npins for pins
f = zeros(ncon,1);
J = zeros(ncon,ndof);
k = 0;
% Some angles and positions are fixed
for i = 1 : nlinks
    link = links(i);
    if link.grounded || link.driver
        % Grounded and driver links have their angles and positions
        % prescribed.
        f(k+1,    1) = link.angle - link.angleTarget;
        f(k+(2:3),1) = link.pos - link.posTarget;
        % The Jacobian of this constraint is the identity matrix
        J(k+(1:3),k+(1:3)) = eye(3);
        k = k + 3;
    end
end
% Pin constraints
for i = 1 : npins
    pin = pins(i);
    rows = k+(1:2); % row index of this pin constraint
    indLinkA = pin.links(1); % array index of link A
    indLinkB = pin.links(2); % array index of link B
    linkA = links(indLinkA);
    linkB = links(indLinkB);
    [Ra,dRa] = rotationMatrix(linkA.angle);
    [Rb,dRb] = rotationMatrix(linkB.angle);
    % Local positions
    ra = pin.pts(:,1);
    rb = pin.pts(:,2);
    % World positions
    xa = Ra * ra + linkA.pos;
    xb = Rb * rb + linkB.pos;
    p = xa(1:2) - xb(1:2);
    f(rows,1) = p;
    % Column indices for the angles and positions of links A and B
    colAngA = (indLinkA-1)*3 + 1;
    colPosA = (indLinkA-1)*3 + (2:3);
    colAngB = (indLinkB-1)*3 + 1;
    colPosB = (indLinkB-1)*3 + (2:3);
    % The Jacobian of this constraint is the partial derivative of f wrt
    % the angles and positions of the two links.
    J(rows,colAngA) = dRa * ra;
    J(rows,colPosA) = eye(2);
    J(rows,colAngB) = -dRb * rb;
    J(rows,colPosB) = -eye(2);
    k = k + 2;
end
end

%%
function particles = updateParticles(links,particles)
% Transform particle position from local to world
for i = 1 : length(particles)
    particle = particles(i);
    link = links(particle.link);
    R = rotationMatrix(link.angle);
    x = R * particle.pt + link.pos;
    % Append world position to the array (grows indefinitely)
    particles(i).ptsWorld(:,end+1) = x;
end
end

%%
function drawScene(links,pins,particles,initialize)
if nargin < 4
    initialize = false;
end
if initialize
    clf;
else
    cla;
end
hold on;
grid on;
% Draw links
for i = 1 : length(links)
    link = links(i);
    R = rotationMatrix(link.angle);
    % Draw frame
    p = link.pos; % frame origin
    s = 0.1; % frame display size
    px = p + s*R(:,1); % frame x-axis
    py = p + s*R(:,2); % frame y-axis
    plot([p(1),px(1)],[p(2),px(2)],'r','LineWidth',3);
    plot([p(1),py(1)],[p(2),py(2)],'g','LineWidth',3);
    % Draw link geometry
    if link.grounded
        color = [1 0 0];
    elseif link.driver
        color = [0 1 0];
    else
        color = [0 0 1];
    end
    E = [R,link.pos;0,0,1]; % transformation matrix
    vertsLocal = [link.verts;ones(1,size(link.verts,2))];
    vertsWorld = E * vertsLocal;
    plot(vertsWorld(1,[1:end,1]),vertsWorld(2,[1:end,1]),'Color',color);
end
% Draw pins
for i = 1 : length(pins)
    pin = pins(i);
    linkA = links(pin.links(1));
    linkB = links(pin.links(2));
    Ra = rotationMatrix(linkA.angle);
    Rb = rotationMatrix(linkB.angle);
    xa = Ra * pin.pts(:,1) + linkA.pos;
    xb = Rb * pin.pts(:,2) + linkB.pos;
    plot(xa(1),xa(2),'co','MarkerSize',10,'MarkerFaceColor','c');
    plot(xb(1),xb(2),'mx','MarkerSize',10,'LineWidth',2);
end
% Draw particles
for i = 1 : length(particles)
    particle = particles(i);
    if ~isempty(particle.ptsWorld)
        plot(particle.ptsWorld(1,:),particle.ptsWorld(2,:),'k');
        plot(particle.ptsWorld(1,end),particle.ptsWorld(2,end),'ko');
    end
end
axis equal;
drawnow;
end
