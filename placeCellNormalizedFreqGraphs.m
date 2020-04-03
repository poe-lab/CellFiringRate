function placeCellNormalizedFreqGraphs(matFile)


% Variables:
%   b = bin number
%   probabilityB = the probability for occupancy of bin b by the subject
%   firingRateB = the mean firing rate of the cell for bin b
%   meanFiringRate = the average firing rate of the cell for all bins.
%   i = cell number
%   j,k = bin index in the matrix
load(matFile, 'Occ', 'TC');
CellOccupancyMatrix = Occ; % A matrix of occupancy calculations for each bin of the maze
cellTotalSpikeMatrix = TC;
clear Occ TC
% probabilityB = [];
% firingRateB = [];
numberOfCells = size(cellTotalSpikeMatrix, 2); % Code for finding length of array
%meanFiringRate = zeros(1,numberOfCells);
%cellInfoContent = zeros(1,numberOfCells);
gridSize = size(CellOccupancyMatrix);
numberRows = gridSize(1);
numberColumns = gridSize(2);

for i = 1:numberOfCells    
    normFreq = cellTotalSpikeMatrix{i}./CellOccupancyMatrix;
    for j=1:numberRows
        for k=1:numberColumns
            if isequalwithequalnans(normFreq(j,k),NaN)
                if isequal(CellOccupancyMatrix(j,k), 0)
                    CellOccupancyMatrixB(j,k) = NaN;
                    normFreq(j,k) = NaN;
                
                end 
            end
        end
    end
    figure(i)
    axes('FontSize',20,'FontName','Arial')
    surf(normFreq)
    axis([0 numberRows 0 numberColumns 0 150 0 20])
    colorbar('location','eastoutside')
    xlabel('Pixels')
    ylabel('Pixels')
    zlabel('Hz')
    title(['Cell ' num2str(i) ' Average Firing Rate'])
end
 for j=1:numberRows
    for k=1:numberColumns
        if isequal(CellOccupancyMatrix(j,k), 0)
            CellOccupancyMatrix(j,k) = NaN;
        end
    end
 end
                                        
figure(i+1)
axes('FontSize',20,'FontName','Arial')
surf(CellOccupancyMatrix)
axis([0 numberRows 0 numberColumns])
colorbar('location','eastoutside')
xlabel('Pixels')
ylabel('Pixels')
zlabel('Hz')
title(['Occupancy Colormap (s)'])


