% Demonstrates how mixed integer programming can be used to simulate the
% dynamics of a bouncing ball 
% see paper: Mote, Mark, J. Pablo Afman, and Eric Feron. "Robotic trajectory planning through collisional interaction." 2017 IEEE 56th Annual Conference on Decision and Control (CDC). IEEE, 2017.

clear; clc; close all;
tic

%% Initialization

N = 1;

% Specify Initial and Final States
pInit  = [1];              % intial position  [px0,py0]
vInit  = [0];              % initial velocity [vx0,vy0]

% Other
dt = 0.01;                 % timestep
tau = 150; 
Mp = 10e3;                 % "Big-M" - Set slightly larger than size of arena

eN = -1.55;

%% Main
clear l
l = linOpt();

%% Define State Variables
py =   l.define('py',tau,N);
vy =   l.define('vy',tau,N);
zeta = l.define('zeta',tau,N,'B');

%% Initial and Final Conditions
for n = 1
    l.declareAtom({[1,py(n,1)]  ,'=', pInit(1)} );
    l.declareAtom({[1,vy(n,1)]  ,'=', vInit(1)} );
    l.declareAtom({[1,zeta(n,1)],'=', 0});
end

%% Dynamics
for n = 1
    for i = 1:tau-1
        l.declareAtom({ [1, py(n, i+1), -1, py(n,i), -dt, vy(n,i),  Mp, zeta(n,i+1)], '>', -(9.8/2)*dt^2});
        l.declareAtom({ [1, py(n, i+1), -1, py(n,i), -dt, vy(n,i), -Mp, zeta(n,i+1)], '<', -(9.8/2)*dt^2});
        l.declareAtom({ [1, vy(n, i+1), -1, vy(n,i),                Mp, zeta(n,i+1)], '>', -(9.8)*dt});
        l.declareAtom({ [1, vy(n, i+1), -1, vy(n,i),               -Mp, zeta(n,i+1)], '<', -(9.8)*dt});
        
        l.declareAtom({ [1, py(n, i+1), -(1+eN), py(n, i), -(1+eN)*dt, vy(n, i), -Mp, zeta(n,i+1)], '>', -Mp });
        l.declareAtom({ [1, py(n, i+1), -(1+eN), py(n, i), -(1+eN)*dt, vy(n, i),  Mp, zeta(n,i+1)], '<',  Mp });
        
        l.declareAtom({ [1, vy(n, i+1), -1-eN, vy(n, i), -Mp, zeta(n,i+1)], '>', 0-Mp });
        l.declareAtom({ [1, vy(n, i+1), -1-eN, vy(n, i),  Mp, zeta(n,i+1)], '<', 0+Mp });
    end
end

%% Switching Condition
for n = 1
    for i = 2:tau-1
        l.declareAtom({ [ 1, py(n,i-1), dt, vy(n,i-1),  Mp, zeta(n,i)], '>', -(9.8/2)*dt^2 });
        l.declareAtom({ [ 1, py(n,i-1), dt, vy(n,i-1),  Mp, zeta(n,i)], '<', -(9.8/2)*dt^2 + Mp });
    end
end

%% Optimization
% Cost
c =  zeros(1,l.numStates);
Q =  zeros(l.numStates,l.numStates);

% Model
model.A = sparse(l.A) ;
model.rhs = l.b ;
model.sense = l.sense;
model.lb = -(10^6).*ones(1,l.numStates);  %[zeros(1,tau+1) -(10^6).*ones(1,nstates-tau-1)]%-(10^6).*ones(1,nstates);
model.Q = sparse(Q);
model.obj = c;
model.vtype = l.vType;

% Generate Result
result = gurobi(model);

if strcmp('OPTIMAL',result.status)
    search=false;
else
    search=true;
    tau = tau + 1;
    fprintf('Time: %.2f\n\n',tau*dt);
end

%% Plotting
if search == false
    MS = 4;
    for n = 1:N
        %     Xpos(n,:) = result.x(px(n,:));
        %     Xvel(n,:) = result.x(vx(n,:));
        Ypos(n,:) = result.x(py(n,:));
        Yvel(n,:) = result.x(vy(n,:));
        zeta(n,:) = result.x(zeta(n,:));
    end
    
    % figure(1)
    figure('units','normalized','outerposition',[0 0 1 0.5])
    
    subplot(1,2,1)
    hold on
    for i = 1:length(Ypos)
        if zeta(i) == 0
            plot(dt*(i-1), Ypos(i), 'b.','markersize',MS)
        else
            plot(dt*(i-1), Ypos(i), 'r.','markersize',MS)
        end
    end
    grid minor
    xlabel('Time ($t$) [s]','interpreter','latex')
    ylabel('Position ($s_y$) [m]','interpreter','latex')
    ylim([-0.1,1.1])
    
    % figure(2)
    subplot(1,2,2)
    hold on
    for i = 1:length(Ypos)
        if zeta(i) == 0
            plot(dt*(i-1), Yvel(i), 'b.','markersize',MS)
        else
            plot(dt*(i-1), Yvel(i), 'r.','markersize',MS)
        end
    end
    grid minor
    xlabel('Time ($t$) [s]','interpreter','latex')
    ylabel('Velocity ($v_y$) [m/s]','interpreter','latex')
    ylim([-4.1,4.1])
else
    disp(" ")
    disp("Not Feasible")
    disp(" ")
end