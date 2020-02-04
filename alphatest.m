function out=alphatest(in)

dflt.alpha = [0.5 5];
dflt.expNmatrix =  9:13;
dflt.reltol = 1e-4;

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

%short vars
alp=in.alpha;
n=in.expNmatrix;
rtl=in.reltol;

hin.spacing_han =@(N,p)((pi^2/2)/(0.5*(max(N))^p(1)))^(1/(p(1)+2)); % spacing function
hin.spacing_par = alp;
hin.boundary_con = true;                                
hin.potential_han = @(x,p) 0.5*abs(x).^p(1);      
hin.potential_par = alp;


for j=1:length(alp)
    figure;
    sgtitle(['\alpha = ' num2str(alp(j)) ', Rel. tol. = ' num2str(rtl)]);
    dx=subplot(1,2,1);
    sx=subplot(1,2,2);
    bohrs=bohrsomm(2^n(end),alp(j));
    hin.potential_par=alp(j);
    hold(dx,'on');
    hold(sx,'on');
    trng=zeros(length(n),2);
    for k=1:length(n)
        %preliminare
        nbohrs=bohrs(1:2^n(k));
        hin.lattice_points=2^n(k);
        %principale
        hout=makeh(hin);
        hres=sort(eig(hout.H));
        %tolleranze
        hrat=hres./nbohrs;
        tol=find((abs(hrat-1))<rtl);
        if isempty(tol)
            trng(k,:)=1;
        else
            trng(k,1)=tol(1);
            trng(k,2)=tol(end);
        end
        %grafica
        plot(dx,hres,'DisplayName',sprintf('n=2^{%d}',n(k)));
    end
    
    plot(dx,bohrs,'DisplayName','B.-S. eig.');
    set(dx,'YScale', 'log');
    set(dx,'XScale', 'log');
    axis(dx,[0 2^(n(end)+1) 0 inf]);
    legend(dx,'Location','southeast');
    title(dx,'Eigenvalues');
    xlabel(dx,"Eigenvalues' index");
    ylabel(dx,'Eigenvalues');
    hold(dx,'off');
    
    plot(sx,n,trng,'.-');
    set(sx,'YScale', 'log');
    axis(sx,[(n(1)-1) n(end)+1 1 inf]);
    legend(sx,['Lower limit';'Upper limit'],'Location','southeast');
    title(sx,'Range of compatible eigens');
    xlabel(sx,'n-th power of 2');
    ylabel(sx,"Eigenvalues' index");
    set(sx,'XTick',n);
    grid minor;
    hold(sx,'off');
    
end

end

function E = bohrsomm(N,alp)

b = 2*alp/(alp + 2);
E = 0.5*(pi*alp/2/beta(3/2,1/alp))^b * (1/2:N-1/2)'.^b;

end