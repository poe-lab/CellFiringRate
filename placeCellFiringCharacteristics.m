
function placeCellFiringCharacteristics(Occ, TC, matFile)


% Variables:
%   b = bin number
%   probabilityB = the probability for occupancy of bin b by the subject
%   firingRateB = the mean firing rate of the cell for bin b
%   meanFiringRate = the average firing rate of the cell for all bins.
%   i = cell number
%   j,k = bin index in the matrix

CellOccupancyMatrix = Occ; % A matrix of occupancy calculations for each bin of the maze
CellFiringRateMatrix = TC;
clear Occ TC
% probabilityB = [];
% firingRateB = [];
numberOfCells = size(CellFiringRateMatrix, 2); % Code for finding length of array
meanFiringRate = zeros(1,numberOfCells);
cellInfoContent = zeros(1,numberOfCells);
gridSize = size(CellOccupancyMatrix);
numberRows = gridSize(1);
numberColumns = gridSize(2);

for i = 1:numberOfCells
% Spatial Information Content of Cell Discharge
    meanFiringRate(i) = sum(sum(CellFiringRateMatrix{1,i}))/(numberRows*numberColumns);
    Q = log2(CellFiringRateMatrix{1,i}/meanFiringRate(i));
    for j = 1:numberRows
        for k = 1:numberColumns
            if isequal(Q(j,k), -Inf)  %If the bin has no firing of the cell, the log2 is -Inf.  We set it =0 so that it does not add or subract information.
                Q(j,k) = 0;
            end
        end
    end
    totalRecordingTime = sum(sum(CellOccupancyMatrix));
    probabilityOccupancyMatrix = CellOccupancyMatrix/totalRecordingTime;
    binInfoContent = probabilityOccupancyMatrix .* (CellFiringRateMatrix{1,i}/meanFiringRate(i)) .* Q;
    cellInfoContent(i) = sum(sum(binInfoContent));
    
% Sparsity of Spatial Firing Distribution
    %Indicates relative proportion of maze on which the cell fired.
    sparsity(i) = sum(sum((probabilityOccupancyMatrix .* CellFiringRateMatrix{1,i} .* CellFiringRateMatrix{1,i})/(meanFiringRate(i)^2)));
    
% Selectivity of Firing in the Field    
    [A,b] = max(CellFiringRateMatrix{1,i});
    [maxFiringRate(i),c] = max(A);
    binMaxFiringRate(i, 1:2) = [b(c) c];
    fieldSelectivity(i) = maxFiringRate(i)/meanFiringRate(i); %A measure of the bin that the cell most prefers to fire.
end
numCells = 1:1:numberOfCells;
M = [numCells' meanFiringRate' cellInfoContent' sparsity' binMaxFiringRate fieldSelectivity'];
%d = {'CellNumber', 'AvgFiringR', 'CellInfoCo', 'Sparsity11','xBinMaxFRa', 'yBinMaxFRa', 'FieldSelec'; M};
s = xlswrite([regexprep(matFile, '.mat', '.xls')], M)


% fid = fopen('TT02.xls','a');
% fprintf(fid,'Cell');
% for i = 1:numberOfCells
%     fprintf(fid,'\t');
%     fprintf(fid,num2str(i));
% end
% fprintf(fid,'\n');
% fprintf(fid,'Average Firing Rate\t');
% fprintf(fid,num2str(meanFiringRate));
% fprintf(fid,'\n');
% fprintf(fid,'Cell Info Content\t');
% fprintf(fid,num2str(cellInfoContent));
% fprintf(fid,'\n');
% fprintf(fid,'Sparsity\t');
% fprintf(fid,num2str(sparsity));
% fprintf(fid,'\n');
% fprintf(fid,'x-coordBinMaxFiringRate\t');
% fprintf(fid,num2str(binMaxFiringRate(:,1)));
% fprintf(fid,'\n');
% fprintf(fid,'y-coordBinMaxFiringRate\t');
% fprintf(fid,num2str(binMaxFiringRate(:,2)));
% fprintf(fid,'\n');
% fprintf(fid,'Field Selectivity\t');
% fprintf(fid,num2str(fieldSelectivity));
% fclose(fid);



