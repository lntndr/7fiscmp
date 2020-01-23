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