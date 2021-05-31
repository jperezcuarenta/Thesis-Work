% This code pertains to Figures 2.3 and 2.4 in thesis manuscript

% This function returns a plot of 
% 1) convolution for increasing values of t.

% CasoN parameter
% if CasoN = 1, graph of non normalized derivatives
% else, graph of normalized derivatives 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% casoN: Scalar value of either 1 or 2
%
%
% Output
% Figures:
% If casoN == 1, show non-normalized derivative figures
% If casoN == 2, show normalized derivative figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function weback2(casoN)
close all
T = 0.1:0.5:2;
M = length(T);
N = 1;    
% sampling interval constant
s_int = 0.1;
% interval for function 'f(t)'
x1 = -1; x2=1;
%  x1 = -4; x2 = 4;
t = [ x1:s_int:x2 ];
% definition of function 'f(t)'
f = 1*ones(1, length(t)); 
%  f = sin(t);
% interval for function 'go(t1)'
t1 = [-8:s_int:8];
% definition of function 'go(t1)'
syms x tSym
symKt = 1/(2*sqrt(tSym))*exp(-pi/(2*tSym)*x^2);
symLaplacian = diff(symKt,x,2);
Laplacian = matlabFunction(symLaplacian);

if casoN == 1
    for jj = 1:length(T)
        tSelect = T(jj);
        go = Laplacian(tSelect,t1);
        % convolve: note the multiplation by the sampling interval
        c = s_int * conv(f, go);
        % flip 'go(t1)' for the graphical convolutions g = go(-t1)
        g = fliplr(go);
        tf = fliplr(-t1);
        % slide range of 'g' to discard non-ovelapping areas with 'f' in the convolution
        tf = tf + ( min(t)-max(tf) );
        % get the range of function 'c' which is the convolution of 'f(t)' and 'go(t1)'
        tc = [ tf t(2:end)];
        tc = tc+max(t1);
        ded1 = -8:x1;
        ded2 = x2:8;
        subplot(M,1,jj)
        plot(tc,c);
        ylim([min(c)-0.2 max(c)+0.2])
        xlim([-3 3])
    end
elseif casoN == 2
    for jj = 1:length(T)
        tSelect = T(jj);
        go = Laplacian(tSelect,t1);
        % convolve: note the multiplation by the sampling interval
        c = tSelect.*s_int * conv(f, go);
        % flip 'go(t1)' for the graphical convolutions g = go(-t1)
        g = fliplr(go);
        tf = fliplr(-t1);
        % slide range of 'g' to discard non-ovelapping areas with 'f' in the convolution
        tf = tf + ( min(t)-max(tf) );
        % get the range of function 'c' which is the convolution of 'f(t)' and 'go(t1)'
        tc = [ tf t(2:end)];
        tc = tc+max(t1);
        ded1 = -8:x1;
        ded2 = x2:8;

        % Plot
        subplot(M,1,jj)
        plot(tc,c);
        ylim([min(c)-0.2 max(c)+0.2])
        xlim([-3 3])
    end
end
 
return  
 

% Ignore
syms x t 

Kt(x,t) = (1/sqrt(2*t))*exp(-(pi*x.^2)/(2*t));
LapKer(x,t) = diff(Kt,x,2);
DxLapKer = diff(LapKer,x,1);
eq1 = DxLapKer == 0;
xCritical = solve(eq1,x);
late
maxLapKer1 = LapKer(xCritical(2),t);
maxLapKer2 = LapKer(xCritical(3),t);

xx1 = linspace(-20,-8);
xx2 = linspace(-8,8);
xx3 = linspace(8,20);

figure;
hold on
plot(xx1,0*xx1,'k','LineWidth',1.5);
plot(xx2,0*xx2+1,'k','LineWidth',1.5);
plot(xx3,0*xx3,'k','LineWidth',1.5);
xline(-8,'LineWidth',1.5);
xline(8,'LineWidth',1.5);

