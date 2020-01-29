function [out]=makeh(in)

dflt.lattice_points = 2^10;         % numero punti
dflt.spacing_han = @(N,p) N/N+0*p;      % spacing function
dflt.spacing_par = 1;
dflt.lap_approx = 0;                % 0 fourier, else diagonal p value
dflt.boundary_con = true;           % 1 periodic, else dirichlet
dflt.potential_han = @(x,p) 0*x*p;  % V function handle
dflt.potential_par = 1;             % V parameters

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

% short variables
N = in.lattice_points;
la = in.lap_approx;
bnd = in.boundary_con;
pp = in.potential_par;
sp = in.spacing_par;
dx = in.spacing_han(N,sp);
x = (-(N-1)/2:(N-1)/2)'*dx;
V = in.potential_han(x,pp);

if la                       
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% diagonal laplace
    % trova alpha
    q = la-1;
    A = repmat(0:q,la,1);
    A(1,:) = 2;
    A(1,1) = 1;
    A(2:end,2:end) = A(2:end,2:end).^repmat(2*(1:(q))',1,q);
    b = zeros(la,1);
    b(2,1) = 1;
    alpha = A\b;
    % calcola diagonali
    D=-alpha(1)*diag(ones(N,1))./2;
    for c=2:la
        D=D+(-alpha(c)*diag(ones(N+1-c,1),c-1));    
    end
    H0 = (1/dx^2/2)*D;
    
    if bnd==false %DBC, else PBC
        alpha=alpha*0;
    end 
    
    for j=2:length(alpha)
        for c=1:j-1
            H0(c,N-(c-1)) = -alpha(j);
        end
    end
    
    H0 = H0 + H0';
    H = H0 + diag(V);
    
else                       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fourier laplace
    k = (2*pi/N)*(-floor(N/2):floor((N-1)/2))';
    k = fftshift(k/dx);
    if bnd                  % PBC
        T = ifft(diag(k.^2/2)*fft(eye(N)));
    else                    % DBC
        k = (pi/dx/(N+1))*(1:N)';
        T = sinft(diag(k.^2/2)*sinft(eye(N)));
    end
    T = real(T+T')/2;
    H = T + diag(V);
    
end

% Numerical computation of the eigenvalues

out.H = H;
out.in = in;