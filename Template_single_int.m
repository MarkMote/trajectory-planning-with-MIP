clear; clc; close all; startup; 
tic
%% Initialization

% Number of Agents 
N = 2; 

% Specify Initial and Final States
pInit  = [0 1 ; 1 0];          % intial position  [px0,py0]
vInit  = [0 0 ; 0 0];          % initial velocity [vx0,vy0]
pFinal = [2 1 ;1 2];           % final position   [pxf,pyf]
vFinal = [0 0 ; 0 0];          % final velocity   [vxf,vyf]
vMax =  0.25;                     % max velocity

% Other
startTime = 1;                 % Search starts at this time and increases t until a feasible soln is found
dt = 0.3;                      % timestep
tau = ceil(startTime/dt);      % converts "startTime" into iterations
result.status = 'Not_Optimal'; % search stops when optimal
search = true;                 % used to find optimal time vector

%% Main loop
while search
    clear l 
    l = linOpt();
    
    %% Define State Variables
    px = l.define('px',tau,N);
    py = l.define('py',tau,N);
    vx = l.define('vx',tau,N);
    vy = l.define('vy',tau,N);
    
    %% Initial and Final Conditions
    for n = 1:N
        l.declareAtom({[1,px(n,1)]  ,'=', pInit(1,n)} );
        l.declareAtom({[1,py(n,1)]  ,'=', pInit(2,n)} );
        
        l.declareAtom({[1,px(n,tau)]  ,'=', pFinal(1,n)});
        l.declareAtom({[1,py(n,tau)]  ,'=', pFinal(2,n)});
    end
    
    %% Dynamics
    for n = 1:N
        for i = 1:tau-1
            l.declareAtom({ [1, px(n, i+1), -1, px(n,i), -dt, vx(n,i)], '=', 0});
            l.declareAtom({ [1, py(n, i+1), -1, py(n,i), -dt, vy(n,i)], '=', 0});
        end
    end
    
    %% Saturations
    P = 40; % Sides of poly used to approximate the norm
    for n = 1:N
        for i = 1:tau
            l.normLeq([vx(n,i), vy(n,i)] , vMax , P); % ||[vx,vy||_2 < vMax approximated with P sided poly
        end
    end
    
    %% Boxes
    % Define Box 
    box = [0.75 1.25 0.75 1.25]; 
    
    % Avoid Boxes 
    for n = 1:N
        for i = 1:tau
            l.avoidBox([px(n,i), py(n,i)], box); 
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
end

%% Plotting
Xpos = zeros(size(px)); Ypos = zeros(size(py)); Xvel = zeros(size(vx)); Yvel = zeros(size(vy));

% Plot
for n = 1:N
    Xpos(n,:) = result.x(px(n,:));
    Xvel(n,:) = result.x(vx(n,:));
    Ypos(n,:) = result.x(py(n,:));
    Yvel(n,:) = result.x(vy(n,:));
end
time = (1:tau).*dt;
figure(1)
hold on
for n = 1:N
    plot(Xpos(n,:),Ypos(n,:),'.')
    plot(pInit(n,1),pInit(n,2),'gx')
    plot(pFinal(n,1),pFinal(n,2),'rx')
end
xlabel('x-position (m)')
ylabel('y-position (m)')
toc
plot([box(1), box(2), box(2), box(1), box(1)],[box(3), box(3), box(4), box(4), box(3)])