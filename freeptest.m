function freeptest

in.lattice_points = 2^10; 
in.boundary_con = true;                     % 1 period, 0 Dirichlet
in.potential_han = @(x,p)0;     %
in.potential_par = 1;
in.spacing = @(N,p)1;
in.compute_eig = true;

hold on

for j=0:2:6
    in.lap_approx=j;
    out=makeh(in);
    plot(out.E);
end
