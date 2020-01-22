function [out] = bellcloset1(in)
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
dflt.potential = @(x)x.^2;

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

% short variable names and basic consistency check

if length(in.x_range)==2
    xr = in.x_range;
    dx = in.lattice_step; % LAVORARE SU SPAZIATURA COME FUNZIONE di N,p
    nd = abs(diff(xr)/dx);
    if (xr(2)-xr(1))>0
        x = xr(1):dx:xr(2)-dx;
    end
else
    error('The selected box boundaries are not valid');
end

% time handling

if in.time_step>0 && in.time_span>=0
    dt = in.time_step;
    ns = ceil(in.time_span/dt);
else
    error('Invalid time span or time step');
end

% deploy boundary conditions

bc = in.boundary_condition;

if bc %periodic
    n = floor(nd/2);
    nn = floor((nd-1)/2);
    k = (2*pi/nd)*(-n:nn);
    k = fftshift(k/dx);
    stepper = @pbcstep;
else
    k = (pi/dx/(nd+1))*(1:nd);
    stepper = @dbcstep;
end

% functions 

if isa(in.potential,'function_handle')
    V = in.potential;
else
    error('The given potential is not a valid function handle');
end

if isa(in.psi_0,'function_handle')
    psi = in.psi_0;
else
    error('The given psi_0 is not a valid function handle');
end

% generate psi0 plot

y = feval(psi,x);
h = plot(x,abs(y).^2,'.-');
xlim(xr);
ylim([0,max(y)]);

c=0;
while c<=ns
    y = exp(-1i*dt*feval(V,x)/2).*y;
    y = stepper(y,k,dt);
    set(h,'YData',abs(y).^2);
    title(strcat('t= ',num2str(c*dt),' s'));
    drawnow limitrate;
    c = c+1;
end

function y = dbcstep(x,k,dt)
    y = sinft(exp(-1i*dt*k.^2/2).*sinft(x));

function y=pbcstep(x,k,dt)
    y = ifft(exp(-1i*dt*k.^2/2).*fft(x));
