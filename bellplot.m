function bellplot(in)

    dflt.in=[];
    dflt.U=[];
    dflt.H=[];
    dflt.x=[];
    dflt.psi=[];
    dflt.plot_wave_skip=1; % if ==0, <0 or >0 
    dflt.plot_E_skip=1;
    dflt.is_wave_squared=true;
    dflt.sxylim1=[0 inf];
    dflt.sxylim2=[0 inf];

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
    sxylim1=in.sxylim1;
    sxylim2=in.sxylim2;
    dt=in.dt;
    
    if in.is_wave_squared
        wav=@sqwave;
        sylabel='Probability';
    else
        wav=@wave;
        sylabel='\psi(0)';
    end
    
    if isempty(U)
        y=psi(1,:);
        kfin=size(psi,1);
        stepper=@psi_ulesstep;
        yestep=@ene_ulesstep;
        xestep=@tim_ulesstep;
    else
        y=U*psi;
        kfin=ceil(diff(in.in.time_span)/in.in.time_step);
        stepper=@psi_ustep;
        yestep=@ene_ustep;
        xestep=@tim_ustep;
    end
    
    %MAIN
    
    bfig=figure;
    stop = uicontrol('style','toggle','string','stop',...
            'units','normalized','position',[.45 .01 .1 .05]);
    
    sx = subplot(1, 2, 1, 'Parent', bfig);
    
    dx = subplot(1, 2, 2, 'Parent', bfig);
    
    hold(sx, 'on');
    ylabel(sx,'Energy');
    plot(sx,x,V,'color','k');
    axis(sx,[x(1) x(end) sxylim1(1) sxylim1(2)]);
    yyaxis(sx,'right');
    axis(sx,[x(1) x(end) sxylim2(1) sxylim2(2)]);
    ylabel(sx,sylabel);
    xlabel(sx,'Position');
    ylabel(dx,'Energy');
    xlabel(dx,'Time');
    wvplot = plot(sx,x,wav(y),'.-');
    
    eplot = animatedline(dx);
    
    % eval for limit
    for k=2:kfin
        sgtitle(['Method: ',in.in.evolution_method,...
            '; t=',num2str(round(xestep(k),2))])
        
        if ~mod(k,pw)
            set(wvplot,'YData',wav(y));
        end
        
        if ~mod(k,pe)
            addpoints(eplot,xestep(k),yestep(k));
        end
        
        drawnow limitrate;
        
        y=stepper(k,y);
        
        if get(stop,'value')
            break
        end
    end
    
    hold off

    function y=psi_ustep(~,yo)
        y=U*yo;
    end

    function y=psi_ulesstep(k,~)
        y=psi(k,:);
    end

    function y=ene_ustep(~)
        y=abs(psi'*(U*psi))/(psi'*psi);
    end

    function y=ene_ulesstep(k)
        y=abs(psi(k,:)*(H*psi(k,:)'))/(psi(k,:)*psi(k,:)');
    end

    function y=tim_ustep(~)
        y=dt*k;
    end

    function y=tim_ulesstep(k)
        y=dt(k);
    end

end

function y=sqwave(y0)
    y=abs(y0).^2;
end

function y=wave(y0)
    y=real(y0);
end
