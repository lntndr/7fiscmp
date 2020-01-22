function harmostest


in.boundary_con = true;                     % 1 period, 0 Dirichlet
in.potential_han = @(x,p) 0.5*p(1)^2*x.^2;     %
in.potential_par = 1;
in.lattice_points = 2^10;                   % numero punti
in.spacing = @(N,p)((pi^2/2)/(0.5*(1*N/2)^2))^(1/4);
in.compute_eig = true;

hold on

for j=0:2:6
    in.lap_approx=j;
    out=makeh(in);
    plot(out.E);
end

