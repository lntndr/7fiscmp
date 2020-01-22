function [out] = bellcloset7(in)
%BELLCLOSET is a basic 1D quantum mechanics simulator
dflt.box_boundaries = [-4 4];
dflt.boundary_condition = 'pbc';    % dbc or Dirichlet for Dirichlet B. C.
                                    % pbc or Periodic for Periodic B. C.
                                    % Case insensitive
dflt.nodes = 2048;
dflt.time_step = 0.001;
dflt.time_span = 10;
dflt.potential = @(x)x.^2;
dflt.psi_0 = @(x)exp(-(x-2).^2) + exp(-(x+2).^2);

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

if in.nodes>0 && floor(in.nodes)==in.nodes
    nd = in.nodes;
else
    error('It is needed at least one node');
end

if length(in.box_boundaries)==2
    bb = in.box_boundaries;
    if (bb(2)-bb(1))>0
        dx = (bb(2)-bb(1))/nd;
        x = (bb(1)+(dx/2)):dx:(bb(2)-(dx/2));
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

if strcmpi(bc,'dbc') || strcmpi(bc,'dirichlet')
    k = (pi/dx/(nd+1))*(1:nd);
    stepper = @dbcstep;
elseif strcmpi(bc,'pbc') || strcmpi(bc,'periodic')
    n = floor(nd/2);
    nn = floor((nd-1)/2);
    k = (2*pi/nd)*(-n:nn);
    k = fftshift(k/dx);
    stepper = @pbcstep;
else
    error('The selected boundary condition is not valid');
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
xlim(bb);
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
