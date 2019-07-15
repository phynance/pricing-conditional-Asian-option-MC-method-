clear
format long
tic
% Monte Carlo on an arithmetic-average fixed-strike conditional Asian put option

%%%%%% Problem and method parameters %%%%%%%%%
sigma = 0.2;  % implied volatility
r = 0.05;     % risk-free interest rate
T = 1.;       % expiry
t = 0.;       % initial time
tau = T-t;    % time-to-maturity        
Dt = 0.001;   % size of time step in MC simulations 
M = 1E5;      % number of asset paths in MC simulations
K = 3.;       % fixed strike    
B = 1.5 ;     % 50% observation barrier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sarray = [ 0 : 0.1 : 4];
II = 0.; % initial value of continuous-arithmetic sum of asset prices above B (it is non-zero with initial time is not 0 )
J = 0. ; % initial value of continuous-arithmetic sum of time whose asset prices above B (it is non-zero with initial time is not 0 )

numS = size(Sarray,2);

%%%%%%  initial pairs of arrays for antithetic variate and lists for storing results %%%%%%
transMean1= zeros(M,1); transMean2= zeros(M,1);
Sfinal = zeros(M,1); Sfinal2 = zeros(M,1);
Parith = zeros(M,1);Parith2 = zeros(M,1); 
Y = zeros(M,1); Y2 = zeros(M,1);
Pmean = zeros(numS,1); Zstd = zeros(numS,1);
stderr = zeros(numS,1); confmc = zeros(numS,2); 
N = round(tau/Dt);  
Spath = zeros( 1,N+1 ) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : numS
    S = Sarray(i);
    parfor n = 1:M   %%%%    execute for-loop iterations in parallel
    %for n = 1:M     %%%%    plain for-loop
        ranVec = randn(1,N);

        % Standard Monte Carlo
            Spath = S*cumprod(1+ r*Dt + sigma*sqrt(Dt)*ranVec + sigma^2 *0.5* Dt* (ranVec .* ranVec - 1) );
            Spath = [S Spath];
            OverBpath1 = Spath(Spath >= B); % filter asset prices above B
            transMean1(n,1) = ( II + sum(OverBpath1)* Dt) / ( J +(size(OverBpath1,2))*Dt ) ;       
            Parith(n,1) = exp(-r*tau)*max( K -transMean1(n,1),0);

        % antithetic path 
            Spath2 = S*cumprod(1+ r*Dt - sigma*sqrt(Dt)*ranVec + sigma^2 *0.5* Dt* (-ranVec .* -ranVec - 1) );
            Spath2 = [S Spath2];
            OverBpath2 = Spath2(Spath2 >= B); % filter asset prices above B
            transMean2(n,1) = ( II + sum(OverBpath2)* Dt) / ( J +(size(OverBpath2, 2))*Dt ) ; 
            Parith2(n,1) = exp(-r*tau)*max( K -transMean2(n,1),0);
   end
   Parith3 = ( Parith + Parith2 )/2;
   Pmean(i,1) = mean(Parith3)                   % simulated option price for each asset price S
   Zstd(i,1) = std(Parith3);                    % standard deviation for each simulated option price
   stderr(i,1) = Zstd(i,1) / sqrt(M);           % standard error for each simulated option price
   confmc(i,1) = Pmean(i,1)-1.96*stderr(i,1);   % lower confidence level
   confmc(i,2) = Pmean(i,1)+1.96*stderr(i,1);   % upper confidence level
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot( Sarray ,  Pmean)

toc
