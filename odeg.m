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