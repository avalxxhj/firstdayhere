run go.m;

fpath=strcat(pwd,'/functions');
fpaths=[{''},{'/save2word'},{'/skt'}];
for f=1:length(fpaths)
    addpath(strcat(fpath,cell2mat(fpaths(f))));
end

file=strcat(pwd,'\exhibits\week2.doc');
delete(file);

%geometric Brownian motion
%parameters
T=1; %number of time periods (years)
m=10; %number of simulations
sigma=.18;
mu=.11;
S0=3830; %value of the S&P 500 on Feb. 3 2021

%loop over n
n=[252/21 252 252*10];
for i=1:length(n)
    dt=T/n(i);
    Smc=repmat(S0,m,n(i)+1); %grid for values of S at each time step
    %simulations (note looping twice, by n and by m)
    BM=sigma*sqrt(dt)*normrnd(0,1,m,n(i));  %Brownian motion using Monte-Carlo draws at each time step
    Smc(:,2:end)=S0*exp(cumsum((mu-.5*sigma^2)*dt+BM,2)); %values of S in the grid
    
    %plot paths
    h=figure('Visible','off');
    hold on
    axl=gca;
    plot(dt*(0:n(i)),Smc);
    scatter(0,S0,50,'k','d','filled');
    hline=refline([0 S0]);
    set(hline,'Color','k','Linestyle','--');
    xlim([0 T]);
    title({strcat('Simulations of Geometric Brownian Motion: n=',int2str(n(i)))},'Fontsize',9,'FontName','Garamond'); 
    ylabel(strcat('S&P 500 (Price)'),'Fontsize',8);
    xlabel({strcat('Time (0 to 1 year)')},'Fontsize',8);
    !taskkill /f /im winword.exe
    save2word(file);   
    hold off;
end


%test of loops
T=1;
n=2520;
dt=T/n;
m=1000; %number of paths/simulations
%with loops (SLOW!)
tic;
Smc_loops=repmat(S0,m,n+1);
for i=1:m %loop over simulations/paths 
   for j=1:n %loop over time steps (sequentially)
       Smc_loops(i,j+1)=Smc_loops(i,j)*exp((mu-.5*sigma^2)*dt+sigma*sqrt(dt)*normrnd(0,1));
   end
end
toc;

%without loops (FAST!)
tic;
Smc_NOloops=repmat(S0,m,n+1);
BM=sigma*sqrt(dt)*normrnd(0,1,m,n);  %Brownian motion using Monte-Carlo draws at each time step (IID)
Smc_NOloops(:,2:end)=S0*exp(cumsum((mu-.5*sigma^2)*dt+BM,2)); %values of S in the grid
toc;




