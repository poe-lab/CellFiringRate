function cellFiringRateByState_v06122017
working_dir=pwd;
%% Call the sleep scored file name
current_dir='C:\';
cd(current_dir);
scoredCheck = 0;
while isequal(scoredCheck, 0)
    [scoredFile, scoredPath] = uigetfile({'*.xls','Excel 1997-2003 File (*.xls)'},...
        'Select the Sleep Scored File');
    if isequal(scoredFile,0) || isequal(scoredPath,0)
        uiwait(errordlg('You need to select a file. Please try again',...
            'ERROR','modal'));
    else
        cd(working_dir);
        stageScoredFile= fullfile(scoredPath, scoredFile);
        %Load sleep scored file:
        try
            [numData, stringData] = xlsread(stageScoredFile);
            scoredCheck = 1;
        catch %#ok<*CTCH>
            % If file fails to load, it will notify user and prompt to
            % choose another file.
            uiwait(errordlg('Check if the scored file is saved in Microsoft Excel format.',...
             'ERROR','modal'));
         scoredCheck = 0;
        end

    end
end

%% Detect if states are in number or 2-letter format:
if isequal(size(numData,2),3)
    
    scoredStates = numData(:,2:3);
    clear numData stringData
else
    scoredStates = numData(:,2);
    clear numData
    stringData = stringData(3:end,3);
    [stateNumber] = stateLetter2NumberConverter(stringData);
    scoredStates = [scoredStates stateNumber];
    clear stateNumber stringData
end

%% Load tetrode data:
current_dir='C:\';
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

waithandle= waitbar(0.2,'Loading spikes ..... ');pause(0.2);
[spikeTimeStamps, spikeCellNumbers] = Nlx2MatSpike(spikeFileName, [1 0 1 0 0], 0, 1, [] );
nonZerosIndex = find(spikeCellNumbers);
spikeCellNumbers = spikeCellNumbers(nonZerosIndex)';
spikeTimeStamps = spikeTimeStamps(nonZerosIndex)';
clear nonZerosIndex
spikeTimeStamps = spikeTimeStamps/1000000;
close(waithandle);
uniqueCells = unique(spikeCellNumbers);
numberOfCells = length(uniqueCells);
averageFiringRate = zeros(8, numberOfCells);

%% Label spikes with states:
labeledSpikeStates = zeros(size(spikeTimeStamps,1),1);
stateTargetInterval=[];
lengthScoredStates =length(scoredStates);
lengthScoredSubStates =[]; %#ok<*NASGU>
for i = 1:8
    isoCount = 1;
    lengthScoredSubStates = lengthScoredStates;
    scoredSubStates = scoredStates;
    for j = 1:lengthScoredSubStates
        if isequal(scoredStates(j,2),i)
            scoredSubStates(j,2) = 1;
        else
            scoredSubStates(j,2) = 0;   
        end  
    end
    firstIsoInd = find(scoredSubStates(:,2)==1, 1);
    if isempty(firstIsoInd)
        scoredSubStates = [];
    elseif firstIsoInd > 1 %First of target stage detected is not at index = 1.
        scoredSubStates(1:firstIsoInd-1,:) = [];
    end
    [lengthScoredSubStates, ~] =size(scoredSubStates); %Recalculate the length of the array due to removal of initial rows.
    if lengthScoredSubStates < 2
        stateTargetInterval = [];
    else
        stateTargetInterval(isoCount,1) = scoredSubStates(1,1);  %#ok<*AGROW>
        %The following FOR loop generates isolated intervals based on user-selected states:
        for j = 2:lengthScoredSubStates   
            if isequal(scoredSubStates(j,2),1)
                if isequal(scoredSubStates(j-1,2),0)
                    stateTargetInterval(isoCount,1) = scoredSubStates(j,1); %Looking at time
                end
                if isequal(j, lengthScoredSubStates)
                    stateTargetInterval(isoCount,2) = scoredSubStates(j,1) + 10; %Looking at time. Need to add 10s to catch last epoch boundary
                end
            elseif isequal(scoredSubStates(j,2),0) && isequal(scoredSubStates(j-1,2),1)
                stateTargetInterval(isoCount,2) = scoredSubStates(j,1); %Looking at time
                isoCount = isoCount + 1;
            end
        end
    end

    [lengthIsoArray, ~] = size(stateTargetInterval);
    if isequal(0,lengthIsoArray)
    else
        tempFiringRate = size(lengthIsoArray,numberOfCells);
        for m = 1:lengthIsoArray % Extract all of the sub-intervals for states containing spikes.
           subIntervalIndx = find(spikeTimeStamps(:) >= stateTargetInterval(m,1) & spikeTimeStamps(:) <= stateTargetInterval(m,2));
           if isempty(subIntervalIndx)
           else
               labeledSpikeStates(subIntervalIndx) = i;
               subIntervalCells = spikeCellNumbers(subIntervalIndx);
               for n = 1:numberOfCells
                   tempFiringRate(m,n) = sum(subIntervalCells== uniqueCells(n))/...
                       (stateTargetInterval(m,2)-stateTargetInterval(m,1));  
               end
               clear subIntervalCells
           end
        end
        averageFiringRate(i,:) = mean(tempFiringRate);
        clear tempFiringRate subIntervalIndx
    end
    clear stateTargetInterval
end

[fileName,ignore2] = uiputfile(['stateLabeledSpikesRat_Day_TT_.mat'],'Save state labeled spikes as:');
save(fileName, 'labeledSpikeStates','spikeTimeStamps','spikeCellNumbers','averageFiringRate','uniqueCells');