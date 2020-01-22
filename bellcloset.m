function [out] = bellcloset(in)
%BELLCLOSET is a basic 1D quantum mechanics simulator
dflt.time_span = 10;
dflt.x_range = [-1 1];
dflt.psi_0 = @(x)exp(-(x-2).^2) + exp(-(x+2).^2);
dflt.lattice_step = 0.01;
dflt.time_step = 0.001;
dflt.lap_approx = 0;
dflt.boundary_condition = true;
dflt.odesolver = 0;
dflt.rel_tol = eps;
dflt.potential_han = @(x,p)x.^2;
dflt.potential_par = 1;
dflt.t_ev_method = 'ode'; % or 'symp' or 'exp'

% input handling and checks
if nargin == 0
    out = dflt;
    return;
end

% fill all missing fields from default
fname = fieldnames(dflt);
for jname = 1:length(fname)
    if ~isfield(in,fname{jname})
        in.(fname{jname}) = dflt.(fname{jname});
    end
end

% time handling

if in.time_step>0 && in.time_span>=0
    dt = in.time_step;
    tspan = ceil(in.time_span/dt);
else
    error('Invalid time span or time step');
end

if length(in.x_range)==2
    xrng = in.x_range;
    dx = in.lattice_step; % LAVORARE SU SPAZIATURA COME FUNZIONE di N,p
    N = abs(diff(xrng)/dx);
    if (xrng(2)-xrng(1))>0
        x = xrng(1):dx:xrng(2)-dx;
    end
else
    error('The selected box boundaries are not valid');
end

if isa(in.potential_han,'function_handle')
    V = in.potential_han;
else
    error('The given potential is not a valid function handle');
end

p = in.potential_par;

if isa(in.psi_0,'function_handle')
    psi0 = in.psi_0;
else
    error('The given psi_0 is not a valid function handle');
end

if % all except symp
    mkh=makeh(in);
    H=mkh.H;
    if % ode
        opts = odeset('reltol',in.reltol);
        [times,psi] = odesolver(@schreqn,tspan,psi0,opts);
    else % exp
        U = expm(-1i*H*dt);
        psi = psi_0(x)';
    end
else % symp
    if in.boundary_condition;
        n = floor(nd/2);
        nn = floor((nd-1)/2);
        k = (2*pi/nd)*(-n:nn);
        k = fftshift(k/dx);
        stepper = @pbcstep;
    else
        k = (pi/dx/(nd+1))*(1:nd);
        stepper = @dbcstep;
    end
end

disp("strawberry");




