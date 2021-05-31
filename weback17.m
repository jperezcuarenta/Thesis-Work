% This code pertains to Figures 4.1-4.4 in thesis manuscript

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% None
%
% Output:
% Figures
% 1) Fig 4.1: Grayscale image f to be studied
% 2) Fig 4.2: Analysis of |Sf| wrt several angles
% 3) Fig 4.2: Analysis of |Sf| at coarse region
% 4) Fig 4.3: Analysis of |Sf| at vertical edge
%
% Notes:
% Requires 'texture.mat' file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Print Fig 4.1
TEXTURE = load('texture.mat');
f = TEXTURE.f;
loc_matrix = [120 128; 48 43; 184 230];
Text_Matrix = {'Pixel 2', 'Pixel 1', 'PIxel 3'};
posShape = [loc_matrix 0*zeros(3,1)+5];

fig = figure;
fMarked = insertShape(f,'FilledCircle',posShape,'Color','red','Opacity',1);
imshow(fMarked);
colormap(gray);
axis off;

%% Graphs Fig 4.1-4.3

% parameters we won't modify as much
alpha_Vec = -0.4:0.2:0.4;
theta_Resolution = 32;
%theta_Resolution = 12;
tau_axis = 0:255;
theta_axis = linspace(0,pi,theta_Resolution);

% pixels to visit
% first is oscillatory region
% second is coarse region
% border region
loc_matrix = [120 128; 48 43; 184 230];

% leave as zero
animate_option = 0;

% overall size of subaxis (M rows and N columns)
M = length(alpha_Vec); N = 1;
% Option N = 2 is for both 3D and 2D view
% N = 2;

% here we will control horizontal alignment when looping 
axis_position = @(jj) [1, 2]+2*(jj-1);

z_storage = zeros(theta_Resolution,length(tau_axis),length(alpha_Vec),3);
% need storage for less computing time 
for mm = 1:length(alpha_Vec)
    current_alpha = alpha_Vec(mm);
    Sf = weback14(f,current_alpha,theta_Resolution,animate_option);        
    for nn=1:3
    z_storage(:,:,mm,nn) = Sf(:,:,loc_matrix(nn,2),loc_matrix(nn,1))';
    end
end

%% plotting
% is there a proper angle for viewing without white lines ?
close all

zmax_Vec = [7.5,1.5,0.3];
caz = 36.9404; cel = 21.8257;
        
fig_Titles = ["4_NumericalOscillatory",...
    "4_NumericalCoarse","4_NumericalBoundary"];

spacing_val = 0.08;

% k_sel selects oscillatory (1), coarse(2) or boundary (3) region
k_sel = 3;

for kk = k_sel:k_sel
    zmax = zmax_Vec(kk);
    current_loc = loc_matrix(kk,:);
    fig = figure;
    for jj = 1:length(alpha_Vec)
        % pick alpha and choose subfigure        
        current_alpha = alpha_Vec(jj);        

        % these are only needed for 
        % evaluating 3D and 2D graphs 
%        ax_pos = axis_position(jj);
%        ax_pos1 = ax_pos(1); ax_pos2 = ax_pos(2);

        z = z_storage(:,:,jj,kk);
        % compute scales
%        Sf = weback14(f,current_alpha,theta_Resolution,animate_option);        
        % store scales at chosen location
%        z = Sf(:,:,current_loc(2),current_loc(1))';

%        plot 3d
%        ax1 = subaxis(M,N,ax_pos1,'Spacing', spacing_val, 'Padding', 0, 'Margin', 0.1);
%        hAxes = axes;
%        H = surf(tau_axis,theta_axis,z);
%        set(H,'LineStyle','none');
%        ax = gca;
%        ax.SortMethod = 'childorder';
%        view(caz,cel);
%        xlim([0 255]);
%        xticks([0 250])
%        xticklabels({'0','250'})
%        xlabel('$\tau$','Interpreter','Latex');
%        yticks([0 pi]);
%        yticklabels({'0', '\pi'})
%        ylabel('$\theta$','Interpreter','Latex');
%        hcb = colorbar;
%        caxis([0 zmax])
%        zlim([0 zmax])
%        ax1.YDir = 'normal';
%        title_Text = 'Prominent Features in 3D';
%        final_Text = strcat(title_Text,{' '},'for $\alpha = ',{' '},num2str(current_alpha),'$');
%        title(final_Text,'Interpreter','Latex')

        % plot 2d
 %       ax2 = subaxis(M,N,ax_pos2,'Spacing', spacing_val, 'Padding', 0, 'Margin', 0.1);
        ax2 = subaxis(M,N,jj,'Spacing', spacing_val, 'Padding', 0, 'Margin', 0.1);        
        imagesc(tau_axis,theta_axis,z);
%        title_Text = 'Prominent Features in 2D';
        title_Text = 'Projection of $|Sf|$ in 2D';
        final_Text = strcat(title_Text,{' '},'for $\alpha = ',{' '},num2str(current_alpha),'$');
        title(final_Text,'Interpreter','Latex');
        ax2.YDir = 'normal';
        xlim([0 255]);
        xticks([0 50 100 150 200 250])
        xticklabels({'0','50','100','150','200','250'})
        xlabel('$\tau$','Interpreter','Latex');
        yticks([0 pi/4 pi/2 3*pi/4 pi]);
        yticklabels({'0','\pi/4','\pi/2','3\pi/4', '\pi'})
        ylabel('$\theta$','Interpreter','Latex');
        hcb = colorbar;
        caxis([0 zmax])
        zlim([0 zmax])
    end
    overall_Text = 'Visualization of $|Sf|$ with Varying $\alpha$ Values';
    sgtitle(overall_Text,'Interpreter','Latex');
end

return

%% Ignore
% if weback17.m is maximized 
% strecth until utf-8 
current_fig_Text = fig_Titles(k_sel);
set(gcf, 'Color', 'w');
export_fig(current_fig_Text,'-eps');
