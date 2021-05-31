% This code pertains to Figure 2.2 in thesis manuscript

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% T:    Scalar value greater than zero
%
% Output
% Figures:
% 1) characteristic function \chi_{[-1,1]}(x)
% 2) laplacian of gaussian kernel
% 3) convolution of 1) and 2)
%
% Notes:
% Example: weback1(0.3)
% Original code modified from here: 
% % https://www.mathworks.com/matlabcentral/fileexchange/4616-animated-convolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function weback1(T)
close all

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
 tSelect = T;
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
  figure;
%  sgtitle(sprintf( 't=%f', tSelect ) )
  subplot(3,1,1)
  %  yline(0,[-8 -2],'b')
  ded1 = -8:x1;
  ded2 = x2:8;
  plot(t,f, 'b');
  hold on
  plot(ded1,0*ded1,'b',ded2,0*ded2,'b');
  title('$\chi(x)$','Interpreter','Latex');
  line([x1 x1],[0 1]);
  line([x2 x2], [0 1]);
  ylim([-0.5 1.5])
  subplot(3,1,2);
  plot( t1, go, 'r');
  ylim([-10 5])
  title('$\partial_{t} K_{t} = \bigtriangleup K_{t}$','Interpreter','Latex')
  
  subplot(3,1,3)
  plot(tc,c);
  title('$\partial_{t} K_{t}*\chi(x)$','Interpreter','Latex')
  ylim([min(c)-0.2 max(c)+0.2])
  xlim([-8 8])
  tSelect
 
  sgtitle('Characteristic Function, Laplacian of Gaussian, and Convolution','Interpreter','Latex','FontSize',14)
  %save gcf 
  txT = '1_LaplacianResponse1.png';
  print(gcf,txT,'-dpng','-r700');
return  
 
% Ignore:
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

