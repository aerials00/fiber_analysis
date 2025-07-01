% Code to analyze fibers in rheoconfocal images
clear
clc
close all

level = 2^8-1; %8 bit image

% Import images from the directory
testcase = 'CTL cells'; %Example of folder name
mainDir = 'maindirectory';
dataDir = [mainDir, testcase, '\'];
codeDir = 'directory where codes are stored';
saveDirectory = [dataDir, 'save_plots' , '\'];
if ~isfolder(saveDirectory)
    mkdir(saveDirectory); %Create folder where the final plots are stored
end

% Enter the number of images to be analyzed
num_images = 5; %To change accordingly

% The threshold list for individual images (Check this by running them
% separately prior to the entire for loop)
thresh_list = [threshold reflection]; 
thresh_GFP_list = [threshold GFP];

fiberOrientation_Radial_all = [];
distance_cell_edge_all = [];
dist_filtered_orientation_all = []; % orientation of the filtered (2um)distance

% Calibration
micron2pixel = 9.1055;
% Threshold for distance: 2 micron from spheroid edge 
dist_thresh = 2*micron2pixel;

visu = 0; %parameter to change between 0 and 1. 
% if it is 0 it processes all the data. If it is 1 we analyse the single experiments. 

% Processing begins here
for k = [] %1:num_images alternatively to process all images. if visu = 1, change k according to the experiment you want to analyse 
    fname_cell = [sprintf('rect%d', k), '_GFP.tif']; % Load the cell GFP image
    fname_coll = [sprintf('rect%d', k), '_refl.tif']; % Load the fiber image
    fname_orientation = ['Orientation_Mask_', fname_coll]; % Load the orientation mask obtained from ImageJ

    img_cell = double(imread([dataDir,   fname_cell]))/level;
    img_coll = double(imread([dataDir,'collagen_test\', fname_coll]))/level;


    
    %% Processing the cell image
    se1 = strel('disk', 5);
    se2 = strel('disk', 1);
    % Binarize
    thresh_GFP = thresh_GFP_list(k);
    mask_cell = imbinarize(img_cell, thresh_GFP); % Threshold the GFP image 
    mask_cell = imclose(mask_cell, se1); % Morpholgical closing
    mask_cell = bwfill(mask_cell, 'holes'); % Fill holes
    mask_cell = bwareafilt(mask_cell, [500, 1e5]); % Area threshold 
    mask_cell_2 = mask_cell; % this mask includes other cell
    mask_cell = imclearborder(mask_cell); % this excludes other cells from mask
  
    %% Find centroid area and boundaries
    s = regionprops(mask_cell, 'Centroid', 'Area');
    cell_centroid = s.Centroid;
    cell_area = s.Area;
    cell_radius = sqrt(cell_area/pi); % Find the equivalent radius
    Bcell = bwboundaries(mask_cell); % Boundary of the cell
    mask_cell = imcomplement(mask_cell);
    mask_cell_2 = imcomplement(mask_cell_2);
    if visu == 1
        figure(1);  clf
        imshowpair(img_cell, img_coll); axis equal xy
        hold on
        caxis([0 0.1])
        cellfun(@(x) plot(x(:,2),x(:,1)), Bcell, 'UniformOutput',0);
        plot(cell_centroid(1), cell_centroid(2), 'o')
        viscircles([cell_centroid(1), cell_centroid(2)], cell_radius)
    end
    %% Process the collagen image
    cd(codeDir)
    % Preprocess the image to enhance edges (optional)
    J = img_coll;
    thresh = thresh_list(k);
    Bw = imbinarize(J, thresh);
    Bw1 = bwareafilt(Bw, [10, 1e5]); % fiber mask

    % Load the orientation mask
    orientation_mask = imread([dataDir, 'Orientation_Masks\', fname_orientation]);

    % Mask the orientation image with fiber mask
    fiber_orientation = Bw1.*orientation_mask .*mask_cell_2; 
    
    if visu ==1
        figure(2); clf
        subplot(1,2,1)
        imshow(J, []);
        subplot(1,2,2);
        imagesc(fiber_orientation); axis equal
        axis tight
    end
 

    %% Orientation analysis

    fiber_orientation(fiber_orientation==0) = NaN; % Getting rid of zeros from multiplication with mask
    [fiberY0, fiberX0] = find(~isnan(fiber_orientation));

    % Converting subscripts [row,col] to linear indices
    ind = sub2ind(size(fiber_orientation), fiberY0, fiberX0);

    % Extract non-zero pixel values and corresponding x and y coordinates
    orient = fiber_orientation(ind);

    % Centroid of cell
    x_center = cell_centroid(1);
    y_center = cell_centroid(2);
    
    % Reference the fibers to cell centroid
    fiberX_CAF = fiberX0-x_center;
    fiberY_CAF = fiberY0-y_center;

    % distance from centroid to the fiber
    dist_center = sqrt(fiberX_CAF.^2+fiberY_CAF.^2);

    % Convert from Cartesian to polar co-ordinates
    [TH, R] = cart2pol(fiberY_CAF,fiberX_CAF);
    TH = -TH;

    ind_q1 = TH<=pi/2 & TH>=0;
    ind_q2 = TH<0 & TH>=-pi/2;
    ind_q3 = TH<-pi/2;
    ind_q4 = TH>=pi/2;

    THq1 = TH(ind_q1);
    THq2 = TH(ind_q2);
    THq3 = TH(ind_q3);
    THq4 = TH(ind_q4);

    % This is required to convert the ImageJ orientation consistent with
    % cart2pol from MATLAB
    TH(ind_q1)  = pi/2 - TH(ind_q1); % Concerns bottom left of image
    TH(ind_q2)  = -(pi/2+TH(ind_q2)); % Concerns bottom right of image
    TH(ind_q3) = -(pi/2+TH(ind_q3)); % Concerns top right of image
    TH(ind_q4) = -(pi/2 -(pi-TH(ind_q4))); % Concerns top left of image

    m1 = tand(orient);
    m2 = tan(TH);

    fiberOrientation_Radial = atand((m2-m1)./(1+m1.*m2));
    distance_cell_edge = dist_center-cell_radius;    %%distance of the pixel from the cell edge 

    ind_filt = find(distance_cell_edge<dist_thresh & distance_cell_edge>=0);   %%filtering for distances from the cell edge below the threshold
    dist_filtered_orientation = fiberOrientation_Radial(ind_filt);


    fiberOrientation_Radial_all = [fiberOrientation_Radial_all; fiberOrientation_Radial(:)]; %%radial orientation of all pixel with respect to the centroid
    distance_cell_edge_all = [distance_cell_edge_all; distance_cell_edge(:)]; %%distance of all pixels from the centroid
    dist_filtered_orientation_all = [dist_filtered_orientation_all; dist_filtered_orientation(:)]; %%orientation of all pixels in the filtered distance with respect to the centroid 


    figure(3);
    % Histogram of relative orientation of all fibers
    subplot(1, num_images, k)
    bin_edges = linspace(-90, 90, 10); % Generate # bin edges from -90 to 90
    h = histogram((fiberOrientation_Radial),bin_edges);
    xlabel('Relative orientation $(\phi)$ [-]', 'Interpreter','latex')
    ylabel('PDF [-]', 'Interpreter','latex')
    h.Normalization = 'probability';
    axis([-90 90 0 0.2]);
    xticks(-90:20:90); % Set x-axis ticks from -90 to 90 with intervals of 10
    axis square
    set(gca, 'FontSize', 18);
    title('All fibers')


    figure(4);
    % Histogram of relative orientation with distance threshold
    subplot(1, num_images, k)
    h = histogram((dist_filtered_orientation),bin_edges);
    xlabel('Relative orientation $(\phi)$ [-]', 'Interpreter','latex')
    ylabel('PDF [-]', 'Interpreter','latex')
    h.Normalization = 'probability';
    axis([-90 90 0 0.2]);
    xticks(-90:20:90); % Set x-axis ticks from -90 to 90 with intervals of 10
    axis square
    set(gca, 'FontSize', 18);

  
%%
     D1 = distance_cell_edge/micron2pixel;
    ind2 = find(D1>=0);
    xedges = linspace(0, max(D1), 11);
    yedges = linspace(-90,90, 11);
    figure(5);
    subplot(1, num_images, k)
    h = bivariate_dist(D1(ind2), (fiberOrientation_Radial(ind2)), 'distance', 'angle', 1, xedges, yedges);
    ax = gca; % Get current axis handle
    ylabel(ax, 'Relative Orientation [$\phi$, $\deg$]', 'Interpreter', 'latex');
    xlabel(ax, 'Distance from R [$\mu$m]', 'Interpreter', 'latex');
    axis([1 max(D1) -90 90]);
    
 
end

%% Plot the distance filtered histogram for all test images
figure(1); clf
h = histogram((dist_filtered_orientation_all),bin_edges);
% title(sprintf('Average Angles: %0.1fhrs', i*delta_T));
xlabel(ax, '$(\phi)$ [$^{\circ}$]', 'Interpreter', 'latex');
ylabel('PDF [-]', 'Interpreter','latex')
h.Normalization = 'probability';
axis([-90 90 0 0.35]);
xticks(-90:20:90); % Set x-axis ticks from -90 to 90 with intervals of 10
ax = gca;
xlabel(ax, '$\phi$ [$^{\circ}$]', 'Interpreter', 'latex');
axis square
set(gca, 'FontSize', 18);
saveas(gcf, [saveDirectory, testcase, '_orientation_distfilt_all.fig'], 'fig');
print(gcf, fullfile(saveDirectory, [testcase, '_orientation_distfilt_all.png']), '-dpng', '-r300');

%% Plot the bivariate distribution for all test images
D = distance_cell_edge_all/micron2pixel;
ind2 = find(D>=0);
figure(2); clf
% xedges = linspace(0, 20, 11);
% yedges = linspace(-90,90, 11);
xedges = linspace(0, 16, 9);  % Define edges for 10 bins in x-axis
yedges = linspace(-90, 90, 10);  % Define edges for 9 bins in y-axis

[h, xe, ye] = bivariate_dist(D(ind2), fiberOrientation_Radial_all(ind2), 'distance', 'angle', 1, xedges, yedges);
ax = gca; % Get current axis handle
set(gcf, 'Position', [100, 100, 800, 700]);  %set figure size
set(gca, 'Position', [0.15, 0.15, 0.7, 0.7]); %adjust axis position
set(gca, 'LineWidth', 1);  % Make axis lines and ticks thicker
set(gca, 'FontSize', 22); 
ylabel(ax, '$\phi$ [$^{\circ}$]', 'Interpreter', 'latex', 'FontSize', 32);
xlabel(ax, '$d$ [$\mu$m]', 'Interpreter', 'latex', 'FontSize', 32);
xticks(0:2:20);  % Set x-axis ticks
yticks(-90:20:90);  % Set y-axis ticks
saveas(gcf, [saveDirectory, testcase, '_bivariate_phi_distance_all.fig'], 'fig');
print(gcf, [saveDirectory, testcase, '_bivariate_phi_distance_all.png'], '-dpng', '-r300');