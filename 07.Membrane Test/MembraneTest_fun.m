function [Cm_fast,Cm_slow,Rm,Rs,Rsquare,Rsquare_adj,RMSE,tau_fast,tau_slow] = MembraneTest_fun(MembraneTest,baseFileName,T)

% function used to calculate Cm and monitor Rs
I_Raw = MembraneTest;

% check if Current unit is pA
if all(I_Raw(1:20,1)<0.01) %if the first all 20 elements < 0.01
    I_pA = I_Raw*1000;
    fprintf(1,'Unit of %s\ncurrent is nA, switched to pA\n',baseFileName)
else
    I_pA = I_Raw;
end

% Global definitions
fontsize = 10;

% Intrinsic properties
size_I_Raw = size(I_Raw,1);

% Time axes
t = 0:T:(size_I_Raw-1) * T;
% show time axes in work space
t_trans = t';

% t_Index = t/T + 1
% t = (Index - 1)*T = Index*T - 1

I_pA_max = max(I_pA);
I_pA_min = min(I_pA);

figure
title('Membrane test','FontSize',fontsize * 1.2)
Ilim = [I_pA_min-25 I_pA_max+25];
plot(t,I_pA)
ylim(Ilim)


% command voltage -5 mV
cmd_V = -5; 

% find negative peak in "I_Raw"
% "loc" means t value, "i" means index
[I_neg_pk_rev,I_neg_loc] = findpeaks(-I_pA(1:(10/T+1)),'MinPeakProminence',10); % only find peak in first 10ms 
I_neg_peak = -I_neg_pk_rev;

% calculate Series Resistance
Rs = cmd_V/I_neg_peak*1000;

% define stable state current as 15-20ms mean value
I_ss = mean(I_pA(15/T+1:20/T+1));

% calculate Membrane Resistance
Rm = cmd_V/I_ss*1000-Rs;

% exponential fit for I_pA between negative peak to 9ms
x = t(1,I_neg_loc:(15/T+1))'; % transverse x to column vector use '
y = I_pA(I_neg_loc:(15/T+1),1); 

%% Fit: 'fitted curve'.
%  Data for 'fitted curve' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Trust-Region';
opts.Display = 'Off';
opts.Normalize = 'on';
opts.Robust = 'LAR';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

%obtain the goodness of fitting
Rsquare = gof.rsquare;
Rsquare_adj = gof.adjrsquare;
RMSE = gof.rmse;

if Rsquare < 0.9
    fprintf(1,'Fitting of %s\nis not successful, please reset the parameter',Name)
    cftool(x,y)% open curve fit tool
end

% Plot fit with data.
figure( 'Name', 'fitted curve' );
h = plot( fitresult, xData, yData );
text(14.5,-50,sprintf('R-sq = %g\nAdjR-sq = %g\nRMSE = %g',gof.rsquare,gof.adjrsquare,gof.rmse))
legend( h, 'Current(pA) vs. Time(ms)', 'fitted curve', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Time(ms)', 'Interpreter', 'none' );
ylabel( 'Current(pA)', 'Interpreter', 'none' );
grid on

% extract two time constants and adjrsquare
MyCoeffs = coeffvalues(fitresult);
Kfast = -MyCoeffs(1,2);
Kslow = -MyCoeffs(1,4);
tau_fast = 1/Kfast;
tau_slow = 1/Kslow;

%% calculate membrane capacitance
Cm_fast = tau_fast*(1/Rs+1/Rm)*1000;
Cm_slow = tau_slow*(1/Rs+1/Rm)*1000;

end