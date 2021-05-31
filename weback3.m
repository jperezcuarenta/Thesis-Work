% This code pertains to Figures 2.5 and 2.6 in thesis manuscript

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% None
%
% Output:
% 1) Toy image with labeled locations where wavelet coefficients
% are computed
% 2) Several evaluations of $|Sf|$ with varying parameter $\alpha$
%
% Notes:
% Requires yellowSquare.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code will read binary square
% then compute scales for 
% large, medium, and small sized square regions 
% 1st row: alpha = -0.4, 
% 2nd row: change location, same alphas
% 3rd row: change location, same alphas

close all
clear all

f = yellowSquare(1);
alpha0 = [-0.6];
alphaVec = [alpha0];

% [x,y] position for insertShape command
pos1 = [60, 56]; % large 
pos2 = [203, 76]; % small 
pos3 = [204, 204]; % med
pos = [pos1; pos3; pos2];
posShape = [pos 0*zeros(3,1)+5];

fMarked = insertShape(f,'FilledCircle',posShape,'Color','red','Opacity',1);
% include border:
fMarked(1,:,:) = 0;
fMarked(end,:,:) = 0;
fMarked(:,1,:) = 0;
fMarked(:,end,:) = 0;

newPos = [pos(1,1)-30, pos(1,2);...
    pos(2,1)-35, pos(2,2);...
    pos(3,1)-30, pos(3,2)];

% adding text
TextVec = {'$P_1$','$P_3$','$P_2$'};
close all
figure; imshow(fMarked);
shift_Variable = 30;
text(newPos(:,1),newPos(:,2),TextVec,'Color','red','Interpreter','Latex',...
    'FontSize',12);
title('Binary Input Image $f$','Interpreter','Latex','FontSize',14)

%%
fMarked_Text = insertText(fMarked,[pos(:,1), pos(:,2)+1],TextVec,...
    'TextColor','red','BoxOpacity',0);
figure;
imshow(fMarked_Text);
title('Binary Input Image $f$','Interpreter','Latex','FontSize',14)

% proceed with comparison of alpha on small, medium, large 
% extrema 

%%
for ii = 2:7
    newAlpha = alpha0 + 0.2*(ii-1);
    alphaVec = [alphaVec, newAlpha];
end

% subplot dimension 
% 3xlength(alphaVec)
M = length(alphaVec);
N = 3;
figure;
%FigH = figure('Position', get(0, 'Screensize'));
%sgtitle('Oscillatory Level for Various $\alpha$','Interpreter','Latex');
subplotCount = [1:3];

for kk = 1:length(alphaVec)
    alphaIn=alphaVec(kk);
    subplotCountIn = subplotCount+3*(kk-1); 
    for ii=1:3
        posInput = pos(ii,:);
%        subplot(M,N,subplotCountIn(ii));

% this works
%        subaxis(M,N,subplotCountIn(ii), 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.1);

% my attempts
        subaxis(M,N,subplotCountIn(ii), 'sh', 0.05,'sv', 0.03, 'Padding', 0, 'Margin', 0.1);
        Sf = localScale(f,alphaIn);
        plot(Sf(:,posInput(2),posInput(1)));
        xticks([0 255]); xlim([0 255]);
        if kk == 4
            mainTitle = strcat('$','\alpha',{' '},...
            {' = '},'0','$');
            title(mainTitle,'Interpreter','Latex');
            set(gca,'fontsize',8)
        else
            mainTitle = strcat('$','\alpha',{' '},...
            {' = '},num2str(alphaIn),'$');
            title(mainTitle,'Interpreter','Latex');
            set(gca,'fontsize',8)
        end
    end
end

return

% Ignore
% font size 8 works somehow 
sgtitle('Illustration of $|Sf|$ for Various Values of $\alpha$','Interpreter','Latex');
set(gcf, 'Color', 'w');
export_fig('1_NumericalAlphaSquare','-eps');