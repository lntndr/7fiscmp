function harmostest

in.lattice_points = 2^11;                              % numero punti
in.spacing_han =@(N,p)((pi^2/2)/(0.5*(1*N/2)^2))^(1/4); % spacing function
in.boundary_con = true;                                % 1 periodic,
in.potential_han = @(x,p) 0.5*p(1)^2*x.^2;             % V function handle

%service vars
N=in.lattice_points;
max_p=9;
trials=15;

%times vars 
fou_times=zeros(trials,1);
dia_times=zeros(trials,max_p-2);

%eig vars
fou_eig=zeros(N,1);

%legend vars

%fourier
in.lap_approx=0;
tic;out=makeh(in);fou_times(1,1)=toc;
fou_eig(:,1)=sort(eig(out.H));
dia_leg=string();

%multidiag
for j=2:max_p
    in.lap_approx=j;
    tic;out=makeh(in);dia_times(1,j-1)=toc;
    dia_eig(:,j-1)=sort(eig(out.H));
    dia_leg(j-1,:)=sprintf('p=%d, %d-diag',j,(j*2)-1); 
end
%analytical
ana_eig = ((0:N-1) + 1/2)';

%first plot
eig_f=figure;
eig_x=axes(eig_f);
xlim([0 2048]);
title("Harm. oscillator's eigenvalues evaluated with different methods");
xlabel("Eigenvalues' index");
ylabel("Eigenvalues");
hold(eig_x,'on');
di_p=plot(eig_x,dia_eig);
set(di_p,{'color' 'DisplayName'}, ...
    [num2cell(jet(max_p-1),2) num2cell(dia_leg)]);
fo_p=plot(eig_x,fou_eig,'g','LineWidth',1, ...
    'DisplayName','Fourier method');
an_p=plot(eig_x,ana_eig,'k','LineWidth',1.5, ...
    'DisplayName','Analytical eigenvalues');
legend('Location','northwest');
hold(eig_x,'off');

%second plot
fe_diff=abs(fou_eig-ana_eig);
ae_diff=abs(dia_eig-repmat(ana_eig,1,max_p-1));
dif_f=figure;
dif_x=axes(dif_f);
xlim([0 2048]);
title("Difference between analytical and numerical eigenvalues");
xlabel("Eigenvalues' index");
ylabel("|numerical-analytical|");
set(dif_x, 'YScale', 'log')
hold(dif_x,'on');
fd_p=plot(dif_x,fe_diff,'g','LineWidth',1,'DisplayName','Fourier method');
di_p=plot(dif_x,ae_diff);
set(di_p,{'color' 'DisplayName'}, ...
    [num2cell(jet(max_p-1),2) num2cell(dia_leg)]);
% set(di_p,{'color'}, num2cell(jet(max_p-1),2));
ep_p=yline(eps,'.k','DisplayName','eps');
legend('Location','northwest');
hold(dif_x,'off');

for k=2:trials
    in.lap_approx=0;
    tic;[~]=makeh(in);fou_times(k,1)=toc;
    for j=2:max_p
        in.lap_approx=j;
        tic;[~]=makeh(in);dia_times(k,j-1)=toc;
    end
end

%time plots
mft=mean(fou_times);
mdt=mean(dia_times);
eft=std(fou_times)./sqrt(trials);
edt=std(dia_times)./sqrt(trials);

tim_f=figure;
tim_x=axes(tim_f);
title("Time needed to set up H matrix");
xlabel("p value");
ylabel("s");
xlim([-0.25,max_p+0.25])
hold(tim_x,'on');
errorbar(2:max_p,mdt,edt,'.g','DisplayName','n-diagonal method');
errorbar(0,mft,eft,'.r','DisplayName','Fourier method');

set(tim_x, 'YScale', 'log')
legend('Location','northeast')
hold(tim_x,'off');
