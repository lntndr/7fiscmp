function sf = sinft(f)
% columnwise sine fourier transform without doubling;
% N.B.: f(0,:)=f(nr+1,:)=0 is implicit but NEITHER are needed;
% hence the array those transform is delta_(n,1) is 
% sin(pi*(1:nr)'/(nr+1));

    persistent auxsin

    [nr,nc] = size(f);
    N = nr+1; 
    %treat special case of a row vector
    if nr == 1
        N = nc+1;
        f = f.';
        nc = 1;
    end
    
    if ~isequal(size(f),size(auxsin))
        auxsin = sin(pi*(1:N-1)'/N)*ones(1,nc);
    end

    g = zeros(N,nc);
    g(1,:) = 0;
    ff = flipud(f);
    g(2:end,:) = auxsin.*(f + ff) + 0.5*(f - ff);
    g = g/sqrt(2*N);
    
    g = fft(g);
    
    n = floor(N/2);
    p = N-2*n;
    fg = flipud(g);
    sf = zeros(N-1,nc);
    sf(1,:) = g(1,:);
    sf(2:2:end,:) = 1i*(g(2:n+p,:) - fg(1:n+p-1,:));
    sf(3:2:end,:) = ones(n-1,1)*sf(1,:) + cumsum(g(2:n,:) + fg(1:n-1,:));
    
    if nr == 1
        sf = sf.';
    end
    
