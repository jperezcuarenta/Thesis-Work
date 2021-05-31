% This code pertains to Figures 3.5-3.8 in thesis manuscript

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% f: Gray scale image 
%
% Output:
% Figures
% 1) Fig 3.5: 3D view of f
% 2) Fig 3.6: Suppression of low scales wrt blurring
% 3) Fig 3.7: Suppression of medium scales wrt blurring
% 4) Fig 3.8: Suppression of large scales wrt blurring
%
% Notes:
% Requires yellowSquare.m
% Requires export_fig.m only for printing figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

if exist('storage_Vec','var') == 0
    % Blurring 
    f = yellowSquare(1);
    f_double = double(f);
    alphaIn = 0;
    
    [M,N] = size(f);
    %sigma_Vec = 0.5:0.05:5;
    sigma_Vec = [0.5:0.05:0.7, 16:20, 61:8:93];
    loop_cnt = length(squeeze(sigma_Vec));

    % Location for rectangles to study 
    pos1 = [60, 56]; % large 
    pos2 = [203, 76]; % small 
    pos3 = [204, 204]; % med
    pos = [pos1; pos3; pos2];
    posShape = [pos 0*zeros(3,1)+5];

    % This vector will store 
    % Sf curves for each pixel and sigma 
    % 1st Dimension: 3 (1 for each pixel)
    % 2nd Dimension: 256 ( to store Sf curve)
    % 3rd Dimension: length(sigma_Vec) (to compute Sf for each sigma)
    storage_Vec = zeros(3,256,loop_cnt);
    f_blur_Vec = zeros(M,N,loop_cnt);
    for jj =1:loop_cnt
        sigma_current = sigma_Vec(jj);
        f_blur = imgaussfilt(f_double,sigma_current);
        f_blur_Vec(:,:,jj) = f_blur;
        Sf = localScale(f_blur,alphaIn);   
    %    h = surf(f_blur);
    %    zlim([0 1])
    %    set(gca,'xtick',[]);     set(gca,'ytick',[])
    %    txt = strcat('\sigma =','{ }',num2str(sigma_current));
    %    set(h,'LineStyle','none'); title(txt);
    %    colorbar;
    %    pause(0.01/2)
        for kk = 1:3
            posInput = pos(kk,:);
            storage_Vec(kk,:,jj) = Sf(:,posInput(2),posInput(1));
        end
    end
else 
    % Pick sigma_Vec adequate for each pixel
    % Pixel 1 is large square (need final entries of sigma_Vec) 
    % Pixel 2 is medium square (need "medium entries" of sigma_Vec)
    % Pixel 3 is small square (need earlier entries of sigma_Vec)
    pixel = 3:-1:1;
    visit_Vec = [1:5; 6:10; 11:15];
    % loop through pixels
    cnt = 1;
    for kk =1:3
        pixel_Sel = pixel(kk);

        % position of sigmas 
        temp_Vec = visit_Vec(kk,:);
        % sampling sigma values 
        sample_sigma_Vec = sigma_Vec(temp_Vec);
        % loop through all sigmas
        loop_cnt2 = length(squeeze(sample_sigma_Vec));
        figure(kk);
            for jj=1:loop_cnt2
    %            figure(cnt);
    %            cnt=cnt+1;
                % Careful here
                % storage_Vec reads the image with less blur 
                % so we can fix the maximum y value
                tempMax = storage_Vec(pixel_Sel,:,temp_Vec(1));
                maxY = max(tempMax(:));
                % read Sf curve at a given pixel and sigma 
                tempSf = storage_Vec(pixel_Sel,:,temp_Vec(jj));
                subaxis(5,1,jj, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.1,'Holdaxis');
                plot(tempSf); xlim([0 255]); ylim([0 maxY])
            end
    end
    
    cnt = 1;
    for kk=4:6
        figure(kk)
        txt = ['Blurring Effects of $K_t * f$']; 
        for jj=1:5
            subaxis(5,1,jj,'Spacing', 0.05, 'Padding', 0, 'Margin', 0.1)
            H = surf(f_blur_Vec(:,:,cnt));            
            set(H,'LineStyle','none'); 

            % Trying something new 
            set(gca, 'Color', 'w');
            ax = gca;
            ax.SortMethod = 'childorder';
            view(caz,cel);
            
            zlim([0 1.2]); colorbar;
            xlim([0 255]); ylim([0 255]);
%            set(gca,'xtick',[]); set(gca,'ytick',[]);
%            txt = strcat('$sigma$ =','{ }',num2str(sigma_current));
            set(H,'LineStyle','none'); 
            if mod(cnt,5) == 1
                title(txt,'Interpreter','Latex','FontSize',14);
            end
            colorbar;
            cnt=cnt+1;
        end
    end    
end 

return 

% window view for storing in thesis 
caz = 29.8643;
cel = 35.2893;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following code will save 
% EPS files for blurred f 
% eps files need to be opened prior to running the following code 
figure_title_Vec = {'3_BlurringEffects_SmallSquare', ...
    '3_BlurringEffects_MediumSquare', '3_BlurringEffects_LargeSquare'};
for jj=4:6
    current_txt = figure_title_Vec(jj-3);
    gc = figure(jj);
%    ax = gca;
%    ax.SortMethod = 'childorder';
%    view(caz,cel);
    set(gc, 'Color', 'w');
    export_fig(string(current_txt),'-eps');
end


% use tabular for Fig 3.3 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to adjust figure vertically before saving 
% the following code will save the original image f in 3D by using 
% the surf command 
% get camera angle and use on all graphs 
H1 = figure;
H = surf(f);
% try surfc
zlim([0 1.2]); colorbar;
xlim([0 255]); ylim([0 255]);
set(H,'LineStyle','none'); 
set(H1, 'Color', 'w');

% cannot graph with these but they erase lines 
%ax = gca;
%ax.SortMethod = 'childorder';

% set camera view for all figures 
view(caz,cel);

%figure;
%fGraph = flipud(255*f);
%M = stdfilt(fGraph);
%meshCanopy(fGraph,M,@spring);
title('Three-dimensional Interpretation of $f$','Interpreter','Latex','FontSize',14)
title('Three-dimensional Interpretation of $f$','Interpreter','Latex')
export_fig('3_3dVisualizationOfImage','-eps');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to print squares with different blur 
% do it in 3D so effect of blur is clear 
fMarked = insertShape(f,'FilledCircle',posShape,'Color','red');
% include border for printing image in manuscript 
fMarked(1,:,:) = 0;
fMarked(end,:,:) = 0;
fMarked(:,1,:) = 0;
fMarked(:,end,:) = 0;