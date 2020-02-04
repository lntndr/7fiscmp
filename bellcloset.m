function [out]=bellcloset(in)
dflt.time_span = [0 100];
dflt.x_range = [-2 2];
dflt.psi_0_han = @(x,p) exp(-x.^2);
dflt.psi_0_par = 2;
dflt.lattice_step = 5e-2;
dflt.lap_approx = 2; % same as makeh
dflt.boundary_con = 0; % same as makeh
dflt.potential_han = @(x,p) 0*x;
dflt.potential_par = 2;
dflt.spacing_han = @(N,p) N/N+p;      % spacing function
dflt.spacing_par = 1;
dflt.evolution_method = 'exp'; % or exp or syp
%ODE only
dflt.ode_solver = @ode113;
dflt.relative_tol = 1e-5;
%EXP+SYMP only
dflt.time_step = 5e-2;

%INPUT MANAGMENT

% Input handling and checks
if nargin == 0
    out = dflt;
    return;
end

% Fill all missing fields from default
fname = fieldnames(dflt);
for i = 1:length(fieldnames(dflt))
    if ~isfield(in,char(fname(i)))
        in.(char(fname(i))) = dflt.(char(fname(i)));
    end
end

tspan=in.time_span;
xrange=in.x_range;
psi0_f=in.psi_0_han;
psi0_p=in.psi_0_par;
a=in.lattice_step;
V=in.potential_han;
p=in.potential_par;
evom=in.evolution_method;
%ODE only
reltol=in.relative_tol;
odesolver=in.ode_solver;
%EXP only
dt=in.time_step;

%%%CONSISTENCY CHECK

%main process

% Number of points in the 1D lattice
N  = (xrange(2)-xrange(1))/a + 1;
if ~isinteger(N)
    N = floor(N);
    a = (xrange(2)-xrange(1))/(N-1); % Adjustment of lattice step if not
end                                  % compatible with the x range
% 1D  grid
x = xrange(1):a:xrange(2);

% Initial wave function
psi0 = psi0_f(x,psi0_p);

if strcmpi(evom,'syp')
    H=[];
    if in.boundary_con
        k = fftshift((2*pi/N)*(-floor(N/2):floor((N-1)/2))'/a);
        Utilde = ifft(diag(exp(-1i*dt*k.^2/2))*fft(eye(N)));
    else
        k = (pi/a/(N+1))*(1:N)';
        Utilde = sinft(diag(exp(-1i*dt*k.^2/2))*sinft(eye(N)));
    end
    U = diag(exp(-1i*dt*V(x,p)/2))*Utilde*diag(exp(-1i*dt*V(x,p)/2));
    psi = psi0';
else
    % CALL MAKEH for H
    in.lattice_points=N;
    makeh_res=makeh(in);
    H=makeh_res.H;
    if strcmpi(evom,'ode') %convert directly in consisetncy check if not exp
        U = [];
        opts = odeset('reltol',reltol);
        [dt,psi] = odesolver(@schreqn,tspan,psi0,opts);
    elseif strcmpi(evom,'exp')
        U = expm(-1i*H*dt);
        psi = psi0';
    end
end

out.in=in;
out.H=H;
out.U=U;
out.x=x;
out.psi=psi;
out.dt=dt;

%%% nested

    function ydot = schreqn(~,y)
        ydot = -1i*H*y;
    end

end