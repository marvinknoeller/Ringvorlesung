% Clear the workspace and figures
clear;
clc;
close all;
% load data
load('cqdata.mat')
extend = double(extend);
num_planes = 1;
if slice_to_0 == 'x'
    [y, z] = meshgrid(linspace(-extend, extend, n1), linspace(-extend, extend, n2));
    x = zeros(size(y));
elseif slice_to_0 == 'y'
    [x, z] = meshgrid(linspace(-extend, extend, n1), linspace(-extend, extend, n2));
    y = zeros(size(x));
elseif slice_to_0 == 'z'
    [x, y] = meshgrid(linspace(-extend, extend, n1), linspace(-extend, extend, n2));
    z = zeros(size(x));
elseif slice_to_0 == "yz"
    [x1, z1] = meshgrid(linspace(-extend, extend, n1), linspace(-extend, extend, n2));
    y1 = zeros(size(x1));
    [x2, y2] = meshgrid(linspace(-extend, extend, n1), linspace(-extend, extend, n2));
    z2 = zeros(size(x2));
    num_planes = 2;
end


if num_planes == 1
    figure;
    for tt = 1 :1: M+1
        % Sample scalar field function
        scalar_f = @(x, y, z) reshape(ui_on_plane(:,tt)+real(us_on_plane(:,tt)),n1,n2);
        surfcolor = scalar_f(x, y, z);
        % Create the surface plot
        surf(x, y, z, surfcolor, 'EdgeColor', 'none');
        % light('Position', [0 0 1]);
        colormap(jet)
        colorbar;
        hold on;
        sphere;
        tri = msh.TRIANGLES;
        pos = msh.POS;
        m = trisurf(tri(:,1:3),pos(:,1),pos(:,2),pos(:,3),(0*ui_on_scat(:,tt)),...c
            'EdgeColor', 'none','FaceAlpha',1);
        material(m, 'shiny');
        axis equal;
        % Set axis limits
        xlim([-extend, extend]);
        ylim([-extend, extend]);
        zlim([-extend, extend]);
        clim([-.4,.4])
        hold off;
        view(3); % Set the view to 3D
        axis off
        ax = gca;
        title(strcat("Time t = ", num2str(double(tt)/(double(M)+1)*double(T))))
        drawnow
    end

else
    figure;
    vI1 = ~isnan(ui_on_plane(1:n1*n2,1));
    vI2 = ~isnan(ui_on_plane(n1*n2+1:end,1));
    vals1 = ui_on_plane(1:n1*n2,:)+real(us_on_plane(1:n1*n2,:));
    vals2 = ui_on_plane(n1*n2+1:end,:)+real(us_on_plane(n1*n2+1:end,:));
    maxvals = max(max(max(vals1(vI1,:)),max(vals2(vI2,:))));
    minvals = min(min(min(vals1(vI1,:)),min(vals2(vI2,:))));
    for tt = 1 :1: M+1
        % Sample scalar field function
        scalar_f1 = @(x, y, z) reshape(ui_on_plane(1:n1*n2,tt)+real(us_on_plane(1:n1*n2,tt)),n1,n2);
        scalar_f2 = @(x, y, z) reshape(ui_on_plane(n1*n2+1:end,tt)+real(us_on_plane(n1*n2+1:end,tt)),n1,n2);
       
        surfcolor1 = scalar_f1(x1, y1, z1);
        surfcolor2 = scalar_f2(x2, y2, z2);
        % Create the surface plot
        surf(x1, y1, z1, surfcolor1, 'EdgeColor', 'none');
        hold on
        surf(x2, y2, z2, surfcolor2, 'EdgeColor', 'none');
        light('Position', [0 -1 0]);
        light('Position', [0 0 1]);
        colormap(jet)
        colorbar;
        hold on;
        sphere;
        tri = msh.TRIANGLES;
        pos = msh.POS;
        m = trisurf(tri(:,1:3),pos(:,1),pos(:,2),pos(:,3),(0*ui_on_scat(:,tt)),...c
            'EdgeColor', 'none','FaceAlpha',1);
        material(m, 'shiny');
        axis equal;
        % Set axis limits
        xlim([-extend, extend]);
        ylim([-extend, extend]);
        zlim([-extend, extend]);
        % clim([-.3678,.3678])
        clim([minvals,maxvals])
        hold off;
        view(3); % Set the view to 3D
        axis off
        ax = gca;
        title(strcat("Time t = ", num2str(double(tt-1)/(double(M))*double(T))))
        drawnow
    end
end
