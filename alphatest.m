function alphatest

in.lattice_points = 2^10;                           % numero punti
in.spacing =@(N,p)((pi^2/2)/(0.5*N^p))^(1/(p+2));   %
in.potential_par = 2;                               % 0 fft else diag
in.boundary_con = true;                             % 1 period, 0 Dirichlet
in.potential_han = @(x,p) 0.5*abs(x).^p(1);           %

hold on

for j=0:2:6
    in.lap_approx=j;
    out=makeh(in);
    plot(out.E);
end

