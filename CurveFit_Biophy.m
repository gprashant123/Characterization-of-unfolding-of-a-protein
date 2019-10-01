function CurveFit_Biophy()
%{
The function imports dataset Unf3.dat and permforms nonlinear regression by
making use of the lsqcurvefit() function.

It makes use of another function myExp.m which specifies the nonlinear
function, determined by temperature dependent thermodynamic parameters.

Finally it prints te estimate and the errors of the parameters
%}


Unf3=importdata('Unf3.dat');  %importing data from Unf3.dat

x=Unf3(:,1);   %Temperature Values
y=Unf3(:,2);    %Signal Values

fx=gradient(y);     %Gradient of the signals
minima=min(fx);     %minimum of gradient gives the point of inflection
indx_min=find(fx==minima);
Tm0=x(indx_min);     %Initializing Tm corresponding to minimum gradient

N=50;
Hm0=2.9*N*10^3;   %Initializing Hm in J/mol         
Cp0=50*N;           %initializing Cp in J/K

XFol=x(42:52);      %choosing points in the folded region
YFol=y(42:52);      
XFol=[ones(11,1) XFol];
paramsFol=inv(XFol'*XFol)*XFol'*YFol;  %Finding linear relation between Folded Signal and Temperature
% paramsFol = XFol\YFol;

XUnf=x(86:101);  %choosing points in the unfolded region
YUnf=y(86:101);
XUnf=[ones(16,1) XUnf];
paramsUnf=inv(XUnf'*XUnf)*XUnf'*YUnf;   %Finding linear relation between Unfolded Signal and Temperature

init=[Tm0 Hm0 Cp0 paramsFol(2) paramsFol(1) paramsUnf(2) paramsUnf(1)];  %Vector of initialized parameters

% Visualise Input Data
figure(1); plot(x,y,'o','MarkerSize',6)


lb = [-inf -inf -inf -inf -inf -inf -inf]; % Lower bounds
ub = [inf inf inf inf inf inf inf]; % Upper bounds

% options for curve fitting
options.Display = 'final';
options.MaxFunEvals = 1000;
options.TolFun = 1e-6;
options.TolX = 1e-6;

% Fit using the defined function
[xfit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit('myExp',init,x,y,lb,ub,options);
ci = nlparci(xfit,residual,'jacobian',jacobian,'alpha',0.318);  %for calculating coonfidence intervals for 1 std dev

S_error= (ci(:,2)-ci(:,1))/2;

% p=length(init);    %number of parameters 
% n=101;  %number of samples
% sigma_cap_sq=resnorm/(n-p);  
% cov_betacap = inv(jacobian'*jacobian + 10^-13*diag(diag(jacobian)))*sigma_cap_sq;  %Calculating Covariance Matrix
% S_error = sqrt(diag(cov_betacap));      %Calculating fitting errors of all the parameters

% Calculate the fit using the calculated parameters

[yfit,pu,pf,H,S,G] = myExp(xfit,x);   

% Plot the fit over the data points
hold on;
plot(x,yfit,'r');
hold on
% plot(XFol,XFol*paramsFol,'b');
StFol = @(x) xfit(5)+xfit(4)*x;   %Folding baseline
StUnf = @(x) xfit(7)+ xfit(6)*x;   %Unfolding baseline
Lim1=[240 360];
Lim2=[345 424];
fplot(StFol,Lim1,'-k');
fplot(StUnf,Lim2,'-m');

% Plot title
title('Signal vs Temperature of a protein having 50 residues','FontWeight','Bold');

legend('Data points', 'Fitted curve' , 'Folding baseline', 'Unfolding baseline');
xlabel('Temperature (K)');
ylabel('Signal');

%Plotting Probability vs Temperature
figure(2)
plot(x,pf,'b');
hold on
plot(x,pu,'k');
legend('Folded population', 'Unfolded population');
title('Probability of folded/unfolded population vs Temperature');
ylabel('Probability');
xlabel('Temperature');

%Plotting Enthalpy vs Temperature
figure(3)
plot(x,H,'g');
hold off
plot(x,H/1000,'g');
txt1='Relation between unfolding enthalpy with temperature';
xlabel('Temperature (K)');
ylabel('Unfolding enthalpy (kJ/mol)');
title(txt1);

%Plotting Gibb's Free Energy vs Temperature
figure(4)
plot(x,G/1000,'k');
txt2="Relation between stability (Gibb's Free Energy change) with Temperature";
title(txt2);
xlabel('Temperature (K)');
ylabel("Gibb's Free energy change (kJ/mol)");

%Plotting Entropy vs Temperature
figure(5)
plot(x,S,'m');
txt3='Relation between unfolding entropy with temperature';
title(txt3);
xlabel('Temperature (K)');
ylabel("Unfolding Entropy (J/mol)");

%Printing Estimated Values
fprintf('Estimated Thermodynamic Parameters:\n');
fprintf('Tm = %f' , xfit(1)); fprintf(' K\n');
fprintf('Hm = %f' ,xfit(2)); fprintf(' J/mol\n');
fprintf('Cp = %f' ,xfit(3)); fprintf(' J/K\n');
fprintf('Sm = Hm/Tm = %f' ,xfit(2)/xfit(1)); fprintf(' J/K\n\n');

%Printing Error values
S_error=full(S_error);      
fprintf('Fitting errors of estimated thermodynamic parameters:\n');
fprintf('For Tm : %ld' ,S_error(1)); + fprintf(' K\n');
fprintf('For Hm : %ld', S_error(2)); + fprintf(' J/mol\n');
fprintf('For Cp : %ld', S_error(3)); + fprintf(' J/K mol\n');
fprintf('For S : %ld', (xfit(2)/xfit(1))*sqrt((S_error(1)/xfit(1))^2 + (S_error(2)/xfit(2))^2)); + fprintf(' J/K mol\n');

end

