function bellplot(in)
    dflt.in=[];
    dflt.U=[];
    dflt.H=[];
    dflt.psi=[];
    dflt.plot_wave_skip=1; % if ==0, <0 or >0 
    dflt.plot_E_skip=1;
    dflt.is_wave_squared=true;

    % fill all missing fields from default
    fname = fieldnames(dflt);
    for jname = 1:length(fname)
        if ~isfield(in,fname{jname})
            in.(fname{jname}) = dflt.(fname{jname});
        end
    end

    U=in.U;
    H=in.H;
    x=in.x;
    psi=in.psi;
    V=in.in.potential_han(x,in.in.potential_par);
    pw=in.plot_wave_skip;
    pe=in.plot_E_skip;
    
    if in.is_wave_squared
        wav=@sqwave;
    else
        wav=@wave;
    end
    
    if isempty(U)
        y=psi(1,:);
        kfin=size(psi,1);
        stepper=@ulesstep;
        estep=@ulesestep;
    else
        y=U*psi;
        kfin=ceil(diff(in.in.time_span)/in.in.time_step);
        stepper=@ustep;
        estep=@uestep;
    end
    
    %MAIN
    
    bfig=figure;
    
    sx = subplot(1, 2, 1, 'Parent', bfig);
    dx = subplot(1, 2, 2, 'Parent', bfig);
    
    stop = uicontrol('style','toggle','string','stop',...
            'units','normalized','position',[.45 .01 .1 .05]);
    
    hold(sx, 'on');
    plot(sx,x,V);
    wvplot = plot(sx,x,wav(y),'.-');
    eplot = animatedline(dx);
    
    % eval for limit
    for k=2:kfin
        
        if ~mod(k,pw)
            set(wvplot,'YData',wav(y));
        end
        
        if ~mod(k,pe)
            addpoints(eplot,k,estep(k));
        end
        
        drawnow;
        
        y=stepper(k,y);
        
        if get(stop,'value')
            break
        end
    end
    
    hold off

    function y=ustep(~,yo)
        y=U*yo;
    end

    function y=ulesstep(k,~)
        y=psi(k,:);
    end

    function y=uestep(~)
        y=abs(psi'*(U*psi))/(psi'*psi);
    end

    function y=ulesestep(k)
        y=abs(psi(k,:)*(H*psi(k,:)'))/(psi(k,:)*psi(k,:)');
    end

end

function y=sqwave(y0)
    y=abs(y0).^2;
end

function y=wave(y0)
    y=real(y0);
end
