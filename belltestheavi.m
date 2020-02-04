% testing script
mhin.time_span = [0 1.5];
mhin.x_range = [-5 5];

mhin.psi_0_han = @(x,p) exp(-x.^2);
mhin.psi_0_par = [0 10];

mhin.lattice_step = 0.05;

mhin.lap_approx = 0;
mhin.boundary_con = 0;

mhin.potential_han = @(x,p)-atan(x)+heaviside(x+1)+pi/2;

mhin.potential_par = [0.1 0.2 5];

mhin.spacing_han = @(N,p) N*0+p;      % spacing function
mhin.spacing_par = 0.05;

%ODE ONLY
mhin.odesolver = @ode113;
mhin.relative_tol = 1e-6;

%EXP+SYMP only
mhin.time_step = 5e-4;

methods=["ode","exp","syp"];

for j=1:3
mhin.evolution_method=char(methods(j));
disp(methods(j));
tic;
mhot=bellcloset(mhin);
toc;
mhot.sxylim1=[0 2*pi];
mhot.sxylim2=[0 1];
mhot.is_wave_squared=true;
bellplot(mhot);
end
