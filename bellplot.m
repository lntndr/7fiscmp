function bellplot(in)
    dflt.in=[];
    dflt.U=[];
    dflt.psi=[];
    dflt.plot_wave_skip=24; % if ==0, <0 or >0 
    dflt.is_wave_squared=true;

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
    pw=in.plot_wave_skip;
    
    if in.is_wave_squared
        wav=@sqwave;
    else
        wav=@wave;
    end
    
    if isempty(U)
        y=psi(1,:);
        kfin=size(psi,1);
        stepper=@ulesstep;
    else
        y=U*psi;
        kfin=ceil(diff(in.in.time_span)/in.in.time_step);
        stepper=@ustep;
    end
    
    stop = uicontrol('style','toggle','string','stop',...
            'units','normalized','position',[.45 .01 .1 .05]);
       
    plot(x,V);
    hold on
    %ylim come diff tra estremi
    ylim([-0.2*max(wav(y)),1.2*max(wav(y))]);
    h = plot(x,wav(y),'.-');

    % eval for limit
    for k=2:kfin
        if ~mod(k,pw)
            set(h,'YData',wav(y));
            drawnow;
        end
        y=stepper(k,y);
        if get(stop,'value')
            break
        end
    end
    
    hold off
    
    %NESTED

    function y=ustep(~,yo)
        y=U*yo;
    end

    function y=ulesstep(k,~)
        y=psi(k,:);
    end

end

function y=sqwave(y0)
    y=abs(y0).^2;
end

function y=wave(y0)
    y=real(y0);
end
