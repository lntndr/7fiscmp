function [out]=bellcloset(in)
dflt.time_span = [0 100];
dflt.x_range = [-2 2];
dflt.psi_0_han = @(x,p) exp(-x.^2);
dflt.psi_0_par = 2;
dflt.lattice_step = 5e-2;
dflt.lap_approx = 2; % same as makeh
dflt.boundary_con = 0; % same as makeh
dflt.potential_han = @(x,p) 0*x;
dflt.potential_par = 2;
dflt.evolution_method = 'exp'; % or exp or syp
%ODE only
dflt.ode_solver = @ode113;
dflt.relative_tol = 1e-5;
%EXP only
dflt.time_step = 5e-2;
%GRAPHICS
dflt.plot_options = 'd'; % 'd': dynamic, 'f': final wave function, 
                         % 'e': energy plot, 'n': no plot
dflt.wave_function_or_probability = 'w'; % 'w': wave function, 'p': probability
dflt.plotframe_skips = 0;

%INPUT MANAGMENT

% Input handling and checks
if nargin == 0
    out = dflt;
    return;
end

% Fill all missing fields from default
fname = fieldnames(dflt);
for i = 1:length(fieldnames(dflt))
    if ~isfield(in,char(fname(i)))
        in.(char(fname(i))) = dflt.(char(fname(i)));
    end
end

tspan=in.time_span;
xrange=in.x_range;
psi0_f=in.psi_0_han;
a=in.lattice_step;
% approx=in.lap_approx;
% bc=in.boundary_con;
V=in.potential_han;
evom=in.evolution_method;
%ODE only
reltol=in.relative_tol;
odesolver=in.ode_solver;
%EXP only
dt=in.time_step;
%GRAPH only
plot_options = in.plot_options;
plot_skips = in.plotframe_skips;
wave_prob = in.wave_function_or_probability;

%%%CONSISTENCY CHECK

%main process

% Number of points in the 1D lattice
N  = (xrange(2)-xrange(1))/a + 1;
if ~isinteger(N)
    N = floor(N);
    a = (xrange(2)-xrange(1))/(N-1); % Adjustment of lattice step if not
end                                  % compatible with the x range
% 1D  grid
x = xrange(1):a:xrange(2);

% Initial wave function
psi0 = psi0_f(x);

% CALL MAKEH for H
in.lattice_points=N;
in.spacing=@(N,p)a;
makeh_res=makeh(in);
H=makeh_res.H;
if strcmpi(evom,'ode') %convert directly in consisetncy check if not exp
    opts = odeset('reltol',reltol);
    [times,psi] = odesolver(@schreqn,tspan,psi0,opts);
    %%%%%%%%%%COPIED PART FOR TESTING%%%%%%%%%DO NOT USE%%%%%%%%%%%%%%%
    legend_word = ["wave function","probability"];
    if wave_prob == 'p'
        psi0 = (abs(psi0)).^2;
        psi = (abs(psi)).^2;
        legend_word = char(legend_word(2));
    else
        legend_word = char(legend_word(1));
    end

    % Plot of the final wave function
    if plot_options == 'f' 
        clf
        yyaxis right
        maxV = max(V(x));
        energy = psi0*(H*psi0')/(psi0*psi0');
        maxE = max(maxV, energy);
        plot(x,V(x))
        axis([x(1) x(end) -9/8*maxE 9/8*maxE])
        hold on
        yline(energy,'--','Color',[0.85 0.33 0.10]);
        ylabel('energy')
        yyaxis left
        plot(x, psi0,'--','Color',[0 0.45 0.74])
        plot(x, real(psi(end,:)),'-')
        axis([x(1) x(end) -max(psi0) max(psi0)])
        title(['Evolution with ', char(odesolver)]);
        legend(['Initial ',legend_word,' (t = ',...
            num2str(tspan(1)),')'],['Final ',legend_word,' (t = ',...
            num2str(tspan(2)),')'],'Potential','Initial energy','Location','southwest');
        ylabel(legend_word)

    % Dynamic plot
    elseif plot_options == 'd' 
        clf
        yyaxis right
        maxV = max(V(x));
        energy = psi0*(H*psi0')/(psi0*psi0');
        maxE = max(maxV, energy);
        plot(x,V(x))
        axis([x(1) x(end) -9/8*maxE 9/8*maxE])
        hold on
        yyaxis left
        h = plot(x, real(psi0),'-');
        N_times = size(psi,1);
        if plot_skips == 0
            plot_skips = 1;
        end
        plot_index = 1:plot_skips:N_times;
        title(sprintf('t = %-8.1f',0))
        axis([x(1) x(end) -max(psi0) max(psi0)])
        set(gcf, 'resize','off');
        stop = uicontrol('style','toggle','string','stop',...
            'units','normalized','position',[.45 .01 .1 .05]);   
        for i = 1:N_times/plot_skips
            delete(h);
            h = plot(x, real(psi(plot_index(i),:)),'-');
            axis([x(1) x(end) -max(psi0) max(psi0)])
            title(sprintf('t = %-8.3f',times(plot_index(i))));
            pause(0.01)
            if get(stop,'value')
                break
            end
        end
        plot(x, real(psi0),'--','Color',[0 0.45 0.74]);
        ylabel(legend_word)
        yyaxis right 
        yline(energy,'--','Color',[0.85 0.33 0.10]);
        ylabel('energy')
        legend(['Final ',legend_word],['Initial ',legend_word],'Potential',...
            'Initial Energy','Location','southwest');

    % Plot of energies (energy vs time)
    elseif plot_options == 'e'
        energy = zeros(size(psi,1),1);
        for i = 1:size(psi,1)
            energy(i) = abs(psi(i,:)*(H*psi(i,:)'))/(psi(i,:)*psi(i,:)');
        end
        plot(times,energy);
        title({'Energy vs time',['(Evolution with ',char(odesolver),', relative '...
        'tolerance = ',num2str(in.relative_tolerance),')']})
        xlabel('time')
        ylabel('energy')

    % Plot of the imaginary part of the computed energies
    elseif plot_options == 'i'
        energy = zeros(size(psi,1),1);
        for i = 1:size(psi,1)
            energy(i) = psi(i,:)*(H*psi(i,:)')/(psi(i,:)*psi(i,:)');
        end
        plot(imag(energy));
        title({'Imaginary part of complex energies',['(Evolution with ',...
            char(odesolver),'relative tolerance = ',num2str(in.relative_tolerance),')']})
        xlabel('Number of computed energies')
        ylabel('Imaginary part')

    % 'no plot' option; in this case the output is the matrix of the solutions 
    % for each time step
    elseif plot_options == 'n' 
        out = [times psi];     
    end
    %%%%%%%%%%END COPIED PART%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(evom,'exp')
    U = expm(-1i*H*dt);
    psi = psi0';
    %%%%%%%%%%COPIED PART FOR TESTING%%%%%%%%%DO NOT USE%%%%%%%%%%%%%%%
    N_times = floor((tspan(2)-tspan(1))/dt);

    legend_word = ["wave function","probability"];
    if plot_options == 'f' % Plot of the final wave function
        clf
        yyaxis right
        maxV = max(V(x));
        energy = psi0*(H*psi0')/(psi0*psi0');
        maxE = max(maxV, energy);
        plot(x,V(x))
        axis([x(1) x(end) -9/8*maxE 9/8*maxE])
        ylabel('energy')
        hold on
        yline(energy,'--','Color',[0.85 0.33 0.10]);
        yyaxis left
        if wave_prob == 'p'
            maxY = max(abs(psi0).^2);
        else
            maxY = max(psi0);    
        end
        axis([x(1) x(end) -maxY maxY])
        if wave_prob == 'p'
            plot(x,abs(psi0).^2,'--','Color',[0 0.45 0.74]);
            legend_word = char(legend_word(2));
        else
            plot(x,real(psi0),'--','Color',[0 0.45 0.74]);
            legend_word = char(legend_word(1));
        end
        ylabel(legend_word)
        for i = 1:N_times
            psi = U*psi;
        end
        if wave_prob == 'p'
            plot(x,abs(psi).^2,'-');
        else
            plot(x, real(psi),'-');
        end
        title('Evolution with matrix exponentiation');
        legend(['Initial ',legend_word,' (t = ',num2str(tspan(1)),')'],... 
            ['Final',legend_word,' (t = ',num2str(tspan(2)),')'],...
            'Potential','Initial energy','Location','best');

    elseif plot_options == 'd' % Dynamic plot
        clf   
        time = tspan(1);
        stop = uicontrol('style','toggle','string','stop',...
                'units','normalized','position',[.45 .01 .1 .05]);
        yyaxis right
        maxV = max(V(x));
        energy = psi0*(H*psi0')/(psi0*psi0');
        maxE = max(maxV, energy);
        plot(x,V(x))
        axis([x(1) x(end) -9/8*maxE 9/8*maxE])
        hold on
        yyaxis left
        if plot_skips == 0
            plot_skips = 1;
        end
        if wave_prob == 'p'
            maxY = max(abs(psi0).^2);
        else
            maxY = max(psi0);    
        end
        h =  plot(x, real(psi),'-');

        for i = 1:floor(N_times/plot_skips)
            for j = 1:plot_skips
                psi = U*psi;
            end
            time = time + plot_skips*dt;
            pause(0.01)
            if wave_prob == 'p'
                delete(h);
                h = plot(x,abs(psi).^2,'-');
            else
                delete(h);
                h = plot(x, real(psi),'-');            
            end
            title(sprintf('t = %-8.3f',time));
            set(gcf, 'resize','off');   
            axis([x(1) x(end) -maxY maxY])
            pause(0.01)
            if get(stop,'value')
                break
            end
        end
        hold on;
        if wave_prob == 'p'
            legend_word = char(legend_word(2));
            plot(x, abs(psi0).^2,'--','Color',[0 0.45 0.74]);
        else
            legend_word = char(legend_word(1));
            plot(x, real(psi0),'--','Color',[0 0.45 0.74]);
        end

        ylabel(legend_word)
        yyaxis right 
        yline(energy,'--','Color',[0.85 0.33 0.10]);
        ylabel('energy')
        legend(['Final ',legend_word],['Initial ',legend_word],'Potential',...
            'Initial energy','Location','southwest');

    % Plot of energies (energy vs time)
    elseif plot_options == 'e'
        energy = zeros(N_times,1);
        energy(1) = psi0*(H*psi0')/(psi0*psi0');
        for i = 1:N_times
            psi = U*psi;
            energy(i+1) = abs(psi'*(H*psi))/(psi'*psi);
        end    
        times_vector = tspan(1):dt:tspan(2);
        plot(times_vector,energy);
        title({'Energy vs time','(Evolution with matrix exponentiation)'})
        xlabel('time')
        ylabel('energy')

    % Plot of the imaginary part of the computed energies    
    elseif plot_options == 'i'
        energy = zeros(N_times,1);
        energy(1) = psi0*(H*psi0')/(psi0*psi0');
        for i = 1:N_times
            psi = U*psi;
            energy(i+1) = psi'*(H*psi)/(psi'*psi);
        end    
        plot(imag(energy));
        title({'Imaginary part of complex energies',['(Evolution with matrix '...
            'exponentiation)']})
        xlabel('Number of computed energies')
        ylabel('Imaginary part')

    % 'no plot' option; in this case the output is the matrix of the solutions 
    % for each time step
    elseif plot_options == 'n'          
        times = tspan(1):dt:tspan(2);   
        out = [psi'; zeros(N_times,size(psi,1))];
        for i = 1:N_times
            psi = U*psi;
            if wave_prob == 'p'
                out(i+1,:) = abs(psi').^2;
            else
                out(i+1,:) = psi';
            end
        end
        out = [times' out];
    end
    %%%%%%%%%%END COPIED PART%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%% nested

    function ydot = schreqn(~,y)
        ydot = -1i*H*y;
    end 

end