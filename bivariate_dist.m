function [h, edgesX, edgesY] = bivariate_dist(x,y, str1, str2, m, xedges, yedges)
% x = Dc; % Random data for demonstration, replace with your data
% y = AR; % Random data for demonstration, replace with your data

% Number of bins for the histogram
% % numBins = 10;
% 
% % Create a 2D histogram
[counts, edgesX, edgesY] = histcounts2(x, y, xedges, yedges, 'Normalization','count');
% 
% % Normalize by column of x
% sumcol = sum(counts, 2);
% normat = repmat(sumcol, [1,numBins]);
% counts = counts./normat;
% % Get the bin centers for plotting
% centerX = edgesX(1:end-1) + diff(edgesX)/2;
% centerY = edgesY(1:end-1) + diff(edgesY)/2;

 % Normalize by the sum of counts in each x-bin
    sumcol = sum(counts, 2);  % Sum over y-bins for each x-bin
    normat = repmat(sumcol, [1, size(counts, 2)]);  % Expand to match counts size
    
    % Avoid division by zero
    normat(normat == 0) = 1;
    counts = counts ./ normat;

    % Compute bin centers
    centerX = edgesX(1:end-1) + diff(edgesX)/2;
    centerY = edgesY(1:end-1) + diff(edgesY)/2;

% Plot the bivariate histogram with log scale for the color
if m==0
    figure;
    h = imagesc(centerX, centerY, log10(counts'));
    set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
    colorbarHandle = colorbar; % Create the colorbar and get its handle
    colorbarHandle.Label.String = 'log(PDF [-])'; % Set the title
    colormap(flipud(cbrewer2('Spectral', 256))); % Choose colormap (optional)
    % title('Bivariate distribution with log-scaled counts');
    xlabel(str1, 'Interpreter','latex');
    ylabel(str2, 'Interpreter','latex');
    clim([0, max(log10(counts(:)))]); % Set color scale range
    axis square
end
if m==1
    %     figure;
    h = imagesc(centerX, centerY, (counts'));
    set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
    colorbarHandle = colorbar; % Create the colorbar and get its handle
    colorbarHandle.Label.String = 'PDF [-]' % Set the title
    colorbarHandle.Label.Interpreter = 'latex';
    colorbarHandle.Label.FontSize = 28;
    colormap('hot'); % Choose colormap (optional)
%     colormap(redblue);
    % title('Bivariate distribution with linear-scaled counts');
    xlabel(str1, 'Interpreter','latex');
    ylabel(str2, 'Interpreter','latex');
    %     clim([0, max((counts(:)))]); % Set color scale range
    axis square
    % Apply fixed caxis for consistent color range across multiple plots
    caxis([0, 0.3]);
end