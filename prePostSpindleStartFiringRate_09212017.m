function prePostSpindleStartFiringRate_09212017

%% Load tetrode data:
working_dir=pwd;
current_dir='C:\SleepData';
cd(current_dir);

[spikeFile, spikePath] = uigetfile({'*.ntt',...
        'Sorted Neuralynx Tetrode File (*.NTT)'},'Select a Spike Sorted Data File');
if isequal(spikeFile,0) || isequal(spikePath,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    spikeFileName = fullfile(spikePath, spikeFile);
end

[spikeTimeStamps, spikeCellNumbers] = Nlx2MatSpike(spikeFileName, [1 0 1 0 0], 0, 1, [] );
clear spikeFileName spikePath
nonZerosIndex = find(spikeCellNumbers);
spikeCellNumbers = spikeCellNumbers(nonZerosIndex)';
spikeTimeStamps = spikeTimeStamps(nonZerosIndex)';
clear nonZerosIndex
spikeTimeStamps = spikeTimeStamps/1000000;
numOfCells = max(spikeCellNumbers);

%% Select cell to analyze:
cellCheck = 0;
while isequal(cellCheck, 0)
    prompt = {['Choose a cell to be analyzed (max=' num2str(numOfCells) ')']};
    def = {'0'};
    dlgTitle = 'Select Cell Number';
    lineNo = 1;
    answer = inputdlg(prompt,dlgTitle,lineNo,def);
    cellNum = str2num(answer{1,1});
    if isempty(cellNum) || sum(cellNum > numOfCells) > 0
        
    else
        cellCheck = 1;
    end
end
targetIdx = ismember(spikeCellNumbers,cellNum);
clear cellCheck prompt def dlgTitle lineNo answer numOfCells
spikeTimeStamps = spikeTimeStamps(targetIdx);
spikeCellNumbers = spikeCellNumbers(targetIdx);
clear targetIdx

%% Load spindle data from .MAT file:
working_dir=pwd;
current_dir='C:\SleepData';
cd(current_dir);

[spindleFile, spindlePath] = uigetfile({'*.mat',...
        'Detected spindles file (*.MAT)'},'Select the spindle data file:');
if isequal(spindleFile,0) || isequal(spindlePath,0)
    uiwait(errordlg('You need to select a file. Please try again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    spindleFileName = fullfile(spindlePath, spindleFile);
end

load(spindleFileName, '-mat')
targetIdx = EschenkoSpindle.scoring==2 | EschenkoSpindle.scoring==6;
spindleStartTS = TimeStamps(EschenkoSpindle.startIdx(targetIdx));
spindlePeakTime = spindleStartTS + (EschenkoSpindle.duration .* EschenkoSpindle.symmetry');
spindleStartTS = spindlePeakTime;
clear EschenkoSpindle spindlePath spindleFileName CSCFilename scoredFile

%% Calculate the firing rate of cells around spindle start time:
numSpindles = size(spindleStartTS, 1);
numOfCells = length(cellNum);
timeWindowSize = .150; % in seconds
firingRateMatrix = zeros(numSpindles, 40, numOfCells);
for i = 1:numSpindles
    for j = 1:40
        binStartTs = spindleStartTS(i) + timeWindowSize * (j - 21);
        binStopTs = spindleStartTS(i) + timeWindowSize * (j - 20);
        spikesInBin = spikeCellNumbers(spikeTimeStamps >= binStartTs & spikeTimeStamps < binStopTs);
        if ~isempty(spikesInBin)
            for k = 1:numOfCells
                firingRateMatrix(i,j,k) = sum(spikesInBin == cellNum(k));
            end
        end
    end
end    
    
firingRateMatrix = firingRateMatrix ./ timeWindowSize;
averageCellFiringRate = squeeze(mean(firingRateMatrix, 1));
averageFiringRate = mean(averageCellFiringRate, 2);
    