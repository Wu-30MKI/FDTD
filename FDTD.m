

% OneD_FDTD_1
% A 1d FDTD simulation of plane wave pulse propagation in vacuum
clear all;
c0 = 3.0e8; % Speed of Light in Vacuum
% Size of the FDTD space
ke = 100;
% Number of Time Steps
nsteps = 2*ke;
% Source Position
ks = 2;
kstart = ks;
% Cell Size & Time Stepping
dt = ke/nsteps;
%Courant Stability Factor
cc = 1.0;             % Courant Condition
dx = c0*dt/cc;  % Courant Condition
epsilon = 1;
cb(kstart:ke/2) = cc./epsilon;  % Enforces Propagation in the positive z-Direction
cb(ke/2+1:ke) = cc./(1.5*epsilon);  % Places an Interface at the Midpoint of the Calculation Space
% Initialize Field Vectors
Ex = zeros(1,ke);
Hy = zeros(1,ke);
% Initial Gaussian Pulse
t0 = 20;
sigma = 10;
omega = 0.5;
% Start Loop
M = moviein(nsteps);
for t = 1:nsteps
    % E Field Loop
    for k = 2:ke-1
        Ex(k) = Ex(k) + cb(k)*(Hy(k-1) - Hy(k));
    end
    % Source Field
    Ex(ks) = Ex(ks) + exp(-((t-t0)/sigma)^2)*cos(omega*(t-t0));
    % H Field Loop
    for k = 1:ke-1
        Hy(k) = Hy(k) + cb(k)*(Ex(k) - Ex(k+1));
    end
    %Ex(2) = Ex(1);          % Left ABC
    Ex(ke) = Ex(ke-1);   % Right ABC
    plot(Ex);axis([1 ke -2 2]);
    xlabel('distance $z$'),ylabel('E Field')
    M(:,t) = getframe   ;
    % input(' ')
end
%movie(M,1);
