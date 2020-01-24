function bellplot(in)
dflt.in=[];
dflt.U=[];
dflt.psi=[];

% fill all missing fields from default
fname = fieldnames(dflt);
for jname = 1:length(fname)
    if ~isfield(in,fname{jname})
        in.(fname{jname}) = dflt.(fname{jname});
    end
end

U=in.U;
x=in.x;
psi=in.psi;
V=in.in.potential_han(x,in.in.potential_par);

if isempty(U)
    y=psi(1,:);
    kfin=size(psi,1);
    stepper=@ulesstep;
else
    y=U*psi;
    kfin=ceil(diff(in.in.time_span)/in.in.time_step);
    stepper=@ustep;
end

plot(x,V);
hold on
ylim([0,max(abs(y).^2)]);
h = plot(x,abs(y).^2,'.-');

% eval for limit

for k=2:kfin
    set(h,'YData',abs(y).^2);
    drawnow nocallbacks;
    y=stepper(k,y);
end

    function y=ustep(~,yo)
        y=U*yo;
    end
    
    function y=ulesstep(k,~)
        y=psi(k,:);
    end

end