%FiringRateMap.m  
%computes firing rate vs location on straightened track for place field
%firing

global ORIG_TRACK_COORDS
global ROTATED_TRACK_COORDS
global ROTATED_VT_data
global CW
global VT_data
global SMOOTH_VT_data
global lap_numbers
global total_num_laps
global vt_pts_in_lap
global TRK_LEN_PXLS
    
global LONG_LEN
global TRK_LONG_LEN_PXLS
global TRK_SHORT_LEN_PXLS
global TOTAL_TRK_LEN
   
global SNAPPED_VT_data
global SN_VT_POS_ON_ARM_PXLS
global SN_ARM_IND
  
working_dir = pwd;
DATADIR = 'C:\Sleepdata\8-4\247\Novel 1813 To 3189s';
%DATADIR = 'C:\Rat 182\2004-8-15RUN';
%DATADIR = 'C:\Rat 159\2004-4-22RUN';
track_coords = [];
pf_info = [];
place_cell_info = [];

% read in .cpf, .trk files and Nvt file
cd(DATADIR);
[cpf_file, cpf_path] = uigetfile({'*.cpf','Place Field File (*.cpf)'},'Select a Place Field File');
place_field_file = fullfile(cpf_path, cpf_file);

rotated_track_coords = csvread(place_field_file,3,0,[3,0,6,5]);
pf_info = csvread(place_field_file,9,0,[9,0,9,2]);
place_cell_info = csvread(place_field_file,12,0);

place_field_start_cm = pf_info(1);
place_field_end_cm = pf_info(2);
TOTAL_TRK_LEN = pf_info(3);

[nvt_file, nvt_path] = uigetfile({'*.Nvt',...
            'Neuralynx Video File (*.Nvt)'},'Select a Corresponding Video Position File');
        vt_file = fullfile(nvt_path, nvt_file);
 
[track_file, track_path] = uigetfile({'*.trk','Track Definition File (*.trk)'},'Select a Track Definition File');
trk_file = fullfile(track_path, track_file);
corners = csvread(trk_file,3,0);

cd(working_dir);
first_vt_tstamp = place_cell_info(1,1) - 5000000;
last_vt_tstamp = place_cell_info(length(place_cell_info),1) + 5000000;
[tstamps xpos ypos] = Read_NVT_File(vt_file,first_vt_tstamp,last_vt_tstamp,1);
    
    % check for correct VT file
    while(abs(mean(xpos)) < 100)
        working_dir = pwd;
        %data_dir = horzcat(ROOT_DIR,'SleepData\Datafiles');
        data_dir = DATADIR;
        cd(DATADIR);
        
        [nvt_file, nvt_path] = uigetfile({'*.Nvt',...
            'Neuralynx Video File (*.Nvt)'},sprintf('Select a different Video Position File for %s',place_field_file));
        vt_file = fullfile(nvt_path, nvt_file);
        DATADIR = nvt_path;
        cd(working_dir);
        [tstamps xpos ypos] = Read_NVT_File(vt_file,first_vt_tstamp,last_vt_tstamp,1);
        
    end
    
% put place cell info into cell array place_cell_firing{lap #, i} where
% i=1 timestamp
% i=2 rotated x
% i=3 rotated y
% i=4 position on straightened track (cm)

lap_ctr = place_cell_info(1,5);
lap_numbers = [];
lap_numbers(1) = lap_ctr;
num_laps = 1;
spike_ctr = 1;
spikes_in_lap = [];
for(i=1:length(place_cell_info))
    if(place_cell_info(i,5) == lap_ctr)
        place_cell_firing{lap_ctr,1}(spike_ctr) = place_cell_info(i,1);
        place_cell_firing{lap_ctr,2}(spike_ctr) = place_cell_info(i,2);
        place_cell_firing{lap_ctr,3}(spike_ctr) = place_cell_info(i,3);
        place_cell_firing{lap_ctr,4}(spike_ctr) = place_cell_info(i,4);
        spike_ctr = spike_ctr + 1;
    else
        spikes_in_lap(lap_ctr) = spike_ctr - 1;
        spike_ctr = 1;
        lap_ctr = place_cell_info(i,5);
        num_laps = num_laps + 1;
        lap_numbers(num_laps) = lap_ctr;
        place_cell_firing{lap_ctr,1}(spike_ctr) = place_cell_info(i,1);
        place_cell_firing{lap_ctr,2}(spike_ctr) = place_cell_info(i,2);
        place_cell_firing{lap_ctr,3}(spike_ctr) = place_cell_info(i,3);
        place_cell_firing{lap_ctr,4}(spike_ctr) = place_cell_info(i,4);
        spike_ctr = spike_ctr + 1;
    end
end
spikes_in_lap(lap_ctr) = spike_ctr - 1;
total_num_laps = num_laps;

% find locations of spikes with least and greatest positions on straightened track
min_pos_spike = min(place_cell_firing{lap_numbers(1),4});
max_pos_spike = max(place_cell_firing{lap_numbers(1),4});
for k=1:total_num_laps
    if(min(place_cell_firing{lap_numbers(k),4}) < min_pos_spike)
            min_pos_spike = min(place_cell_firing{lap_numbers(k),4});
    end
    if(max(place_cell_firing{lap_numbers(k),4}) > max_pos_spike)
            max_pos_spike = max(place_cell_firing{lap_numbers(k),4});
    end
end

% put VT info in cell array VT_data{lap #, i} where 
% i=1 timestamp
% i=2 raw x
% i=3 raw y
% include VT data for 1 sec pre first spike and post last spike, unless ...
% position on straightened track is 20 pxls beyond min_pos_spike or 20 pxls
% less than max_pos_spike, then include 2 secs pre first or post last
%j = 1;
%if(place_cell_firing{lap_numbers(j),4}(1) > min_pos_spike+20)
%    lap_start_time = place_cell_firing{lap_numbers(j),1}(1) - 2000000;
%else
%    lap_start_time = place_cell_firing{lap_numbers(j),1}(1) - 1000000;
%end
%if(place_cell_firing{lap_numbers(j),4}(spikes_in_lap(lap_numbers(j))) < max_pos_spike-20)
%    lap_end_time = place_cell_firing{lap_numbers(j),1}(spikes_in_lap(lap_numbers(j))) + 2000000;
%else
%    lap_end_time = place_cell_firing{lap_numbers(j),1}(spikes_in_lap(lap_numbers(j))) + 1000000;
%end
%k = 1;
%vt_pts_in_lap = [];
%for(i=1:length(tstamps))
%    if(tstamps(i) <= lap_end_time)
%        if(tstamps(i) >= lap_start_time)
%            VT_data{lap_numbers(j),1}(k) = tstamps(i);
%            VT_data{lap_numbers(j),2}(k) = xpos(i);
%            VT_data{lap_numbers(j),3}(k) = ypos(i);
%            k = k + 1;
%        end
%    else
%        vt_pts_in_lap(lap_numbers(j)) = k - 1
%        k = 1;
%        j = j + 1;
%        if(j > length(lap_numbers))
%            continue
%        end
%        if(place_cell_firing{lap_numbers(j),4}(1) > min_pos_spike+20)
%            lap_start_time = place_cell_firing{lap_numbers(j),1}(1) - 2000000;
%        else
%            lap_start_time = place_cell_firing{lap_numbers(j),1}(1) - 1000000;
%        end
%        if(place_cell_firing{lap_numbers(j),4}(spikes_in_lap(lap_numbers(j))) < max_pos_spike-20)
%            lap_end_time = place_cell_firing{lap_numbers(j),1}(spikes_in_lap(lap_numbers(j))) + 2000000;
%        else
%            lap_end_time = place_cell_firing{lap_numbers(j),1}(spikes_in_lap(lap_numbers(j))) + 1000000;
%        end
%                
%    end
%end

for k=1:total_num_laps
    if(place_cell_firing{lap_numbers(k),4}(1) > min_pos_spike+20)
        lap_start_time = place_cell_firing{lap_numbers(k),1}(1) - 2500000;
    else
        lap_start_time = place_cell_firing{lap_numbers(k),1}(1) - 1000000;
    end
    %%%%%%%%%%%%%
    %if(k == 11)
    %    lap_start_time = place_cell_firing{lap_numbers(k),1}(1) - 500000;
    %end
    %%%%%%%%%%%%%
    if(place_cell_firing{lap_numbers(k),4}(spikes_in_lap(lap_numbers(k))) < max_pos_spike-20)
        lap_end_time = place_cell_firing{lap_numbers(k),1}(spikes_in_lap(lap_numbers(k))) + 2500000;
    else
        lap_end_time = place_cell_firing{lap_numbers(k),1}(spikes_in_lap(lap_numbers(k))) + 2000000;
    end
    
    iinds = find((tstamps >= lap_start_time) & (tstamps <= lap_end_time));
    VT_data{lap_numbers(k),1} = tstamps(iinds);
    VT_data{lap_numbers(k),2} = xpos(iinds);
    VT_data{lap_numbers(k),3} = ypos(iinds);
    vt_pts_in_lap(lap_numbers(k)) = length(iinds);
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% liinearly interpolate missing (0,0) positions in RAW VT data

for k=1:total_num_laps
    for i=1:vt_pts_in_lap(lap_numbers(k))
        if((VT_data{lap_numbers(k),2}(i) == 0) && (VT_data{lap_numbers(k),3}(i) == 0))  
            prev_ind = i-1;
            next_ind = i+1;
            %while previous position is at (0,0) - find previous index
            while((prev_ind > 1) && (VT_data{lap_numbers(k),2}(prev_ind) == 0) && (VT_data{lap_numbers(k),3}(prev_ind) == 0))
                prev_ind = prev_ind - 1;
            end
            
            %while next position is at (0,0) - find next index
            while((next_ind < vt_pts_in_lap(lap_numbers(k))) && (VT_data{lap_numbers(k),2}(next_ind) == 0) && (VT_data{lap_numbers(k),3}(next_ind) == 0))
                next_ind = next_ind + 1;
            end
            
            %if last point, use previous 2 points for interpolation
            if(i == vt_pts_in_lap(lap_numbers(k)))
                next_ind = prev_ind -1;
            end
            
            %if first point, use next 2 points for interpolation - problem
            %if first few points are zero, then they won't be removed
            if(i == 1)
                prev_ind = next_ind + 1;
            end
            
            ts_diff_tot = VT_data{lap_numbers(k),1}(next_ind) - VT_data{lap_numbers(k),1}(prev_ind);
            for j=prev_ind+1:next_ind-1
                this_ts_diff = VT_data{lap_numbers(k),1}(j) - VT_data{lap_numbers(k),1}(prev_ind);
                frac_1 = this_ts_diff/ts_diff_tot;
                frac_2 = 1 - frac_1;
                VT_data{lap_numbers(k),2}(j) = frac_1 * VT_data{lap_numbers(k),2}(next_ind) + frac_2 * VT_data{lap_numbers(k),2}(prev_ind);
                VT_data{lap_numbers(k),3}(j) = frac_1 * VT_data{lap_numbers(k),3}(next_ind) + frac_2 * VT_data{lap_numbers(k),3}(prev_ind);
            end
        end
    end
end
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth out raw VT data replacing each point with a 15 pt average in a
% moving window centered at the point

for k=1:total_num_laps
    for i=1:vt_pts_in_lap(lap_numbers(k))
        
        sum_start = i - 7;
        if (i <= 8)
            sum_start = 1;
        end
        sum_end = i + 7;
        if (sum_end > vt_pts_in_lap(lap_numbers(k)))
            sum_end = vt_pts_in_lap(lap_numbers(k));
        end
        
        xsum = sum(VT_data{lap_numbers(k),2}(sum_start:sum_end));
        ysum = sum(VT_data{lap_numbers(k),3}(sum_start:sum_end));
        
        %check for zero positions
        zero_ctr = 0;
        if ((~all(VT_data{lap_numbers(k),2}(sum_start:sum_end))) && (~all(VT_data{lap_numbers(k),3}(sum_start:sum_end))))
            for j=sum_start:sum_end
                if ((VT_data{lap_numbers(k),2}(j) == 0) && (VT_data{lap_numbers(k),3}(j) == 0))
                    zero_ctr = zero_ctr + 1;
                end
            end
        end
        
        num_in_avg = sum_end - sum_start + 1 - zero_ctr;
        SMOOTH_VT_data{lap_numbers(k),1}(i) = VT_data{lap_numbers(k),1}(i);
        SMOOTH_VT_data{lap_numbers(k),2}(i) = xsum/num_in_avg;
        SMOOTH_VT_data{lap_numbers(k),3}(i) = ysum/num_in_avg;
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transform VT_data to rotated x,y 

ORIG_TRACK_COORDS{1,1} = corners(:,1)';
ORIG_TRACK_COORDS{1,2} = corners(:,2)';
ORIG_TRACK_COORDS{2,1} = corners(:,3)';
ORIG_TRACK_COORDS{2,2} = corners(:,4)';
ORIG_TRACK_COORDS{3,1} = corners(:,5)';
ORIG_TRACK_COORDS{3,2} = corners(:,6)';
 
    trk_ofst_ang = 0;
    trk_anchor = 0;
    % find the index of the closest corner to the origin.
    min_val = 99999999999999;
    min_ind = 0;
    for i=1:4
        if(ORIG_TRACK_COORDS{2,1}(i).^2+ORIG_TRACK_COORDS{2,2}(i).^2 < min_val)
            min_ind = i;
            min_val = ORIG_TRACK_COORDS{2,1}(i).^2+ORIG_TRACK_COORDS{2,2}(i).^2;
        end
    end
    
    if((min_ind == 1) || (min_ind == 3))
        trk_ofst_ang = atan( (ORIG_TRACK_COORDS{2,2}(min_ind+1) - ORIG_TRACK_COORDS{2,2}(min_ind)) /...
                             (ORIG_TRACK_COORDS{2,1}(min_ind+1) - ORIG_TRACK_COORDS{2,1}(min_ind)) );
        trk_anchor = min_ind;
    else
        trk_ofst_ang = atan( (ORIG_TRACK_COORDS{2,2}(min_ind-1) - ORIG_TRACK_COORDS{2,2}(min_ind)) /...
                             (ORIG_TRACK_COORDS{2,1}(min_ind-1) - ORIG_TRACK_COORDS{2,1}(min_ind)) );
        trk_anchor = min_ind-1;                
    end
    
    %rotate the track
    for i=1:4
        if (i == trk_anchor)
            ROTATED_TRACK_COORDS{2,1}(i) = ORIG_TRACK_COORDS{2,1}(i);
            ROTATED_TRACK_COORDS{2,2}(i) = ORIG_TRACK_COORDS{2,2}(i);
            for track = [1 3]
                zeroedX = ORIG_TRACK_COORDS{track,1}(i)-ORIG_TRACK_COORDS{2,1}(trk_anchor);
                zeroedY = ORIG_TRACK_COORDS{track,2}(i)-ORIG_TRACK_COORDS{2,2}(trk_anchor);
            
                [phase,mag] = cart2pol(zeroedX,zeroedY);
                phase = phase - trk_ofst_ang;
                [x,y] = pol2cart(phase,mag);
                ROTATED_TRACK_COORDS{track,1}(i) = ORIG_TRACK_COORDS{2,1}(trk_anchor) + x;
                ROTATED_TRACK_COORDS{track,2}(i) = ORIG_TRACK_COORDS{2,2}(trk_anchor) + y;
            end 
        else
            for track = [1 2 3]
                zeroedX = ORIG_TRACK_COORDS{track,1}(i)-ORIG_TRACK_COORDS{2,1}(trk_anchor);
                zeroedY = ORIG_TRACK_COORDS{track,2}(i)-ORIG_TRACK_COORDS{2,2}(trk_anchor);
            
                [phase,mag] = cart2pol(zeroedX,zeroedY);
                phase = phase - trk_ofst_ang;
                [x,y] = pol2cart(phase,mag);
                ROTATED_TRACK_COORDS{track,1}(i) = ORIG_TRACK_COORDS{2,1}(trk_anchor) + x;
                ROTATED_TRACK_COORDS{track,2}(i) = ORIG_TRACK_COORDS{2,2}(trk_anchor) + y;
            end
        end
    end
    %plot_rotated_plus(handles);
    
    % now... readjust track!!!  This needs to be done because the points
    % don't exactly straiten after the transformation (rounding error).
    horz_corn = 2;
    vert_corn = 4;
    diag_corn = 3;
    if(trk_anchor == 3)
        horz_corn = 4;
        vert_corn = 2;
        diag_corn = 1;
    end

    for track = 1:3
        %fprintf('\nlocked corner: (%4.1f, %4.1f)\n',ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,1}(trk_anchor),ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,2}(trk_anchor));
        %fprintf('Horizontal corner (pre): (%4.1f, %4.1f)\n',ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,1}(horz_corn),ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,2}(horz_corn));
        ROTATED_TRACK_COORDS{track,2}(horz_corn) = ROTATED_TRACK_COORDS{track,2}(trk_anchor);
        %fprintf('Horizontal corner (post): (%4.1f, %4.1f)\n',ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,1}(horz_corn),ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,2}(horz_corn));
        
        %fprintf('Vertical corner (pre): (%4.1f, %4.1f)\n',ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,1}(vert_corn),ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,2}(vert_corn));
        ROTATED_TRACK_COORDS{track,1}(vert_corn) = ROTATED_TRACK_COORDS{track,1}(trk_anchor);
        %fprintf('Vertical corner (post): (%4.1f, %4.1f)\n',ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,1}(vert_corn),ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,2}(vert_corn));
       
        %fprintf('Diag corner (pre): (%4.1f, %4.1f)\n',ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,1}(diag_corn),ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,2}(diag_corn));
        ROTATED_TRACK_COORDS{track,2}(diag_corn) = ROTATED_TRACK_COORDS{track,2}(vert_corn);
        ROTATED_TRACK_COORDS{track,1}(diag_corn) = ROTATED_TRACK_COORDS{track,1}(horz_corn);
        %fprintf('Diag corner (post): (%4.1f, %4.1f)\n',ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,1}(diag_corn),ROTATED_TRACK_COORDS{ACTIVE_CLUST_NUM,track,2}(diag_corn));
    end
    
    % compute rotated VT data
    for j=1:total_num_laps
        for i=1:vt_pts_in_lap(lap_numbers(j))
            
            zeroedX = SMOOTH_VT_data{lap_numbers(j),2}(i)-ORIG_TRACK_COORDS{2,1}(trk_anchor);
            zeroedY = SMOOTH_VT_data{lap_numbers(j),3}(i)-ORIG_TRACK_COORDS{2,2}(trk_anchor);
            
            [phase,mag] = cart2pol(zeroedX,zeroedY);
            phase = phase - trk_ofst_ang;      
            [x,y] = pol2cart(phase,mag);
            ROTATED_VT_data{lap_numbers(j),1}(i) = ORIG_TRACK_COORDS{2,1}(trk_anchor) + x;
            ROTATED_VT_data{lap_numbers(j),2}(i) = ORIG_TRACK_COORDS{2,2}(trk_anchor) + y;
        end      
    end
    TRK_SHORT_LEN_PXLS = abs(ROTATED_TRACK_COORDS{2,2}(3) - ROTATED_TRACK_COORDS{2,2}(2));
    TRK_LONG_LEN_PXLS = abs(ROTATED_TRACK_COORDS{2,1}(2) - ROTATED_TRACK_COORDS{2,1}(1));
    TRK_LEN_PXLS = 2*(TRK_SHORT_LEN_PXLS + TRK_LONG_LEN_PXLS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot lap vt positions and spike positions for each lap
figure(1);
clf;
plot(ROTATED_VT_data{lap_numbers(1),1},ROTATED_VT_data{lap_numbers(1),2}*(-1),place_cell_firing{lap_numbers(1),2},place_cell_firing{lap_numbers(1),3}*(-1),'kx');
hold on;
for k=2:total_num_laps
    plot(ROTATED_VT_data{lap_numbers(k),1},ROTATED_VT_data{lap_numbers(k),2}*(-1),place_cell_firing{lap_numbers(k),2},place_cell_firing{lap_numbers(k),3}*(-1),'kx');
end
plot( [ROTATED_TRACK_COORDS{1,1} ROTATED_TRACK_COORDS{1,1}(1)],...
      -[ROTATED_TRACK_COORDS{1,2} ROTATED_TRACK_COORDS{1,2}(1)],'k -');
plot( [ROTATED_TRACK_COORDS{2,1} ROTATED_TRACK_COORDS{2,1}(1)],...
      -[ROTATED_TRACK_COORDS{2,2} ROTATED_TRACK_COORDS{2,2}(1)],'r -');
plot( [ROTATED_TRACK_COORDS{3,1} ROTATED_TRACK_COORDS{3,1}(1)],...
      -[ROTATED_TRACK_COORDS{3,2} ROTATED_TRACK_COORDS{3,2}(1)],'k -');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    snap_to_track(handles);           
% snap VT data to track
  
    minX = [];
    maxX = [];
    minY = [];
    maxY = [];
    % find the min + max (x,y) for each leg of the track.
    for j=1:4
        next_ind = j+1;
        if(next_ind > 4)
            next_ind = 1;
        end
        if( ROTATED_TRACK_COORDS{2,1}(j) < ROTATED_TRACK_COORDS{2,1}(next_ind) )
            maxX(j) = ROTATED_TRACK_COORDS{2,1}(next_ind);
            minX(j) = ROTATED_TRACK_COORDS{2,1}(j);
        else
            minX(j) = ROTATED_TRACK_COORDS{2,1}(next_ind);
            maxX(j) = ROTATED_TRACK_COORDS{2,1}(j);
        end
        if( ROTATED_TRACK_COORDS{2,2}(j) < ROTATED_TRACK_COORDS{2,2}(next_ind) )
            maxY(j) = ROTATED_TRACK_COORDS{2,2}(next_ind);
            minY(j) = ROTATED_TRACK_COORDS{2,2}(j);
        else
            minY(j) = ROTATED_TRACK_COORDS{2,2}(next_ind);
            maxY(j) = ROTATED_TRACK_COORDS{2,2}(j);
        end
    end
    
    % find closest leg + snap
for k=1:total_num_laps    
    for i=1:vt_pts_in_lap(lap_numbers(k))
        deltaX = [];
        deltaY = [];
        total_dists = [];
        closestX = 0;
        closestY = 0;
        for j=1:4
            if( ( ROTATED_VT_data{lap_numbers(k),1}(i) > minX(j) ) && ( ROTATED_VT_data{lap_numbers(k),1}(i) < maxX(j) ) )
                % within x-range for leg - x distance = 0
                deltaX(j) = 0;
            else
                % find closer x corner, assign the x difference to be the x
                % distance.
                [deltaX(j),theind] = min([abs(maxX(j) - ROTATED_VT_data{lap_numbers(k),1}(i)), abs(minX(j) - ROTATED_VT_data{lap_numbers(k),1}(i))]);
                if(theind == 1)
                    closestX = maxX(j);
                else
                    closestX = minX(j);
                end
            end
            if( ( ROTATED_VT_data{lap_numbers(k),2}(i) > minY(j) ) && ( ROTATED_VT_data{lap_numbers(k),2}(i) < maxY(j) ) )
                % within y-range for leg - y distance = 0
                deltaY(j) = 0;
            else
                % find closer y corner, assign the y difference to be the y
                % distance.
                [deltaY(j),theind] = min([abs(maxY(j) - ROTATED_VT_data{lap_numbers(k),2}(i)), abs(minY(j) - ROTATED_VT_data{lap_numbers(k),2}(i))]);
                if(theind == 1)
                    closestY = maxY(j);
                else
                    closestY = minY(j);
                end
            end
            dist = (deltaX(j)) + (deltaY(j));
            total_dists(i,j) = dist;
        end
        [val,min_ind] = min(total_dists(i,:));
        
        % min_ind is now the probable arm of the track.  If either deltaX
        % or deltaY is zero, just snap the other value to the track.
        switch min_ind
            case 1      % arm 1 - horizontal - y val same as y of corner 1 & 2
                SNAPPED_VT_data{lap_numbers(k),2}(i) = ROTATED_TRACK_COORDS{2,2}(1);
                if(deltaX(min_ind) == 0)
                    SNAPPED_VT_data{lap_numbers(k),1}(i) = ROTATED_VT_data{lap_numbers(k),1}(i);
                else
                    SNAPPED_VT_data{lap_numbers(k),1}(i) = closestX;
                end
                SN_VT_POS_ON_ARM_PXLS{lap_numbers(k)}(i) = abs((SNAPPED_VT_data{lap_numbers(k),1}(i) - ROTATED_TRACK_COORDS{2,1}(1)));
                
                
            case 2      % arm 2 - vertical - x val same as x of corner 2 & 3
                SNAPPED_VT_data{lap_numbers(k),1}(i) = ROTATED_TRACK_COORDS{2,1}(2);
                if(deltaY(min_ind) == 0)
                    SNAPPED_VT_data{lap_numbers(k),2}(i) = ROTATED_VT_data{lap_numbers(k),2}(i);
                else
                    SNAPPED_VT_data{lap_numbers(k),2}(i) = closestY;
                end
                SN_VT_POS_ON_ARM_PXLS{lap_numbers(k)}(i) = abs( (SNAPPED_VT_data{lap_numbers(k),2}(i) - ROTATED_TRACK_COORDS{2,2}(2)));
                
            case 3      % arm 3 - horizontal - y val same as y of corner 3 & 4
                SNAPPED_VT_data{lap_numbers(k),2}(i) = ROTATED_TRACK_COORDS{2,2}(3);
                if(deltaX(min_ind) == 0)
                    SNAPPED_VT_data{lap_numbers(k),1}(i) = ROTATED_VT_data{lap_numbers(k),1}(i);
                else
                    SNAPPED_VT_data{lap_numbers(k),1}(i) = closestX;
                end
                SN_VT_POS_ON_ARM_PXLS{lap_numbers(k)}(i) = abs( (SNAPPED_VT_data{lap_numbers(k),1}(i) - ROTATED_TRACK_COORDS{2,1}(3)));
                
            case 4      % arm 4 - vertical - x val same as x of corner 4 & 1
                SNAPPED_VT_data{lap_numbers(k),1}(i) = ROTATED_TRACK_COORDS{2,1}(4);
                if(deltaY(min_ind) == 0)
                    SNAPPED_VT_data{lap_numbers(k),2}(i) = ROTATED_VT_data{lap_numbers(k),2}(i);
                else
                    SNAPPED_VT_data{lap_numbers(k),2}(i) = closestY;
                end
                SN_VT_POS_ON_ARM_PXLS{lap_numbers(k)}(i) = abs( (SNAPPED_VT_data{lap_numbers(k),2}(i) - ROTATED_TRACK_COORDS{2,2}(4)));
        end
        SN_ARM_IND{lap_numbers(k)}(i) = min_ind;     
    end    
end
%END function snap_to_track(handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute position on straightened track
% function straighten_cluster(handles)
    
for k=1:total_num_laps    
    for i=1:vt_pts_in_lap(lap_numbers(k))
    
        switch SN_ARM_IND{lap_numbers(k)}(i)
            case 1
                tot_num_of_pxls = SN_VT_POS_ON_ARM_PXLS{lap_numbers(k)}(i);
            case 2
                tot_num_of_pxls = SN_VT_POS_ON_ARM_PXLS{lap_numbers(k)}(i) + TRK_LONG_LEN_PXLS;
            case 3
                tot_num_of_pxls = SN_VT_POS_ON_ARM_PXLS{lap_numbers(k)}(i) + TRK_LONG_LEN_PXLS...
                                  + TRK_SHORT_LEN_PXLS;
            case 4
                tot_num_of_pxls = SN_VT_POS_ON_ARM_PXLS{lap_numbers(k)}(i) + 2* TRK_LONG_LEN_PXLS...
                                  + TRK_SHORT_LEN_PXLS;
        end
        
        SN_VT_TOT_POS{lap_numbers(k)}(i) = tot_num_of_pxls;
        %fprintf('Arm: %d\t|\tPxls: %4.1f\t|\tPercent: %4.1f\n',SN_ARM_IND{ACTIVE_CLUST_NUM}(i),tot_num_of_pxls,100*SN_CLUST_TOT_POS{ACTIVE_CLUST_NUM}(i));
        
    end
end
    
%END function straighten_cluster(handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% determine spike positions on straightened track 

for k=1:total_num_laps
    tt_ctr = 1;
    for i=1:spikes_in_lap(lap_numbers(k))
        while (SMOOTH_VT_data{lap_numbers(k),1}(tt_ctr) < place_cell_firing{lap_numbers(k),1}(i))
            prev_ts = SMOOTH_VT_data{lap_numbers(k),1}(tt_ctr);
            tt_ctr = tt_ctr + 1;
        end
        % if timestamps are equal && SN_VT_TOT_POS is not zero
        if((SMOOTH_VT_data{lap_numbers(k),1}(tt_ctr) == place_cell_firing{lap_numbers(k),1}(i)) && (SN_VT_TOT_POS{lap_numbers(k)}(tt_ctr) ~= 0))
            place_cell_sn_pos_tot{lap_numbers(k)}(i) = SN_VT_TOT_POS{lap_numbers(k)}(tt_ctr);
        else
            prev_ind = tt_ctr-1;
            next_ind = tt_ctr;
            %while previous position is at (0,0) - find previous index
            while((prev_ind > 1) && (SN_VT_TOT_POS{lap_numbers(k)}(prev_ind) == 0) )
                prev_ind = prev_ind - 1;
            end
            
            %while next position is at (0,0) - find next index
            while((next_ind < length(SN_VT_TOT_POS{lap_numbers(k)})) && (SN_VT_TOT_POS{lap_numbers(k)}(next_ind) == 0) )
                next_ind = next_ind + 1;
            end
            
            ts_diff_tot = SMOOTH_VT_data{lap_numbers(k),1}(next_ind) - SMOOTH_VT_data{lap_numbers(k),1}(prev_ind);
            this_ts_diff = place_cell_firing{lap_numbers(k),1}(i) - SMOOTH_VT_data{lap_numbers(k),1}(prev_ind);
            frac_1 = this_ts_diff/ts_diff_tot;
            frac_2 = 1 - frac_1;
            place_cell_sn_pos_tot{lap_numbers(k)}(i) = frac_1 * SN_VT_TOT_POS{lap_numbers(k)}(next_ind) + frac_2 * SN_VT_TOT_POS{lap_numbers(k)}(prev_ind);
            
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine place field - defined by earliest spike and latest spike
% fired by cell
place_field_start_pxl = min(place_cell_sn_pos_tot{lap_numbers(1)});
place_field_end_pxl = max(place_cell_sn_pos_tot{lap_numbers(1)});
for k=1:total_num_laps
    if(min(place_cell_sn_pos_tot{lap_numbers(k)}) < place_field_start_pxl)
        place_field_start_pxl = min(place_cell_sn_pos_tot{lap_numbers(k)});
    end
    if(max(place_cell_sn_pos_tot{lap_numbers(k)}) > place_field_end_pxl)
        place_field_end_pxl = max(place_cell_sn_pos_tot{lap_numbers(k)});
    end
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill up firing_rate_distribution{lap #,1}(1:num_incrs) with mid-points of
% place field bins
% fill up firing_rate_distribution{lap #,7}(1:num_incrs) with start locations of
% each place field bins
% fill up firing_rate_distribution{lap #,8}(1:num_incrs) with end locations of
% each place field bins
% fill up firing_rate_distribution{lap #,4}(1:num_incrs) with 0 - number of
% spikes in bin

place_field_size_pxl = place_field_end_pxl - place_field_start_pxl;

%pxl_incr = 7;
pxl_incr = 5;
place_field_bin_start = place_field_start_pxl-(pxl_incr/2);
num_incrs = ceil((place_field_end_pxl-place_field_start_pxl+pxl_incr)/pxl_incr);

for k=1:total_num_laps
    for i=1:num_incrs
        % midpoint of bin
        firing_rate_distribution{lap_numbers(k),1}(i) = place_field_bin_start + (i-0.5)*pxl_incr;
        % start position of bin
        firing_rate_distribution{lap_numbers(k),7}(i) = place_field_bin_start + (i-1)*pxl_incr;
        % end position of bin
        firing_rate_distribution{lap_numbers(k),8}(i) = place_field_bin_start + i*pxl_incr;
        % initialize number of spikes in bin
        firing_rate_distribution{lap_numbers(k),4}(i) = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine bin for each spike and increment spike counters for each bin - firing_rate_distribution{lap_numbers(k),4}(i)

for k=1:total_num_laps
    for i=1:spikes_in_lap(lap_numbers(k))
        bin_num = ceil((place_cell_sn_pos_tot{lap_numbers(k)}(i) - place_field_bin_start)/pxl_incr);
        firing_rate_distribution{lap_numbers(k),4}(bin_num) = firing_rate_distribution{lap_numbers(k),4}(bin_num) + 1;    
    end
    fprintf('for lap %4.1f, spikes counted %4.1f  actual number of spikes %4.1f\n',lap_numbers(k),sum(firing_rate_distribution{lap_numbers(k),4}),spikes_in_lap(lap_numbers(k)))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine time in each bin - put in firing_rate_distribution{lap #,3}
% interpolate time for starting position and end position of each bin -
% firing_rate_distribution{lap #,5} and firing_rate_distribution{lap #,6}
% Then compute firing rate in each bin - put in firing_rate_distribution{lap #,2}

v_lap_numbers = [];
num_v_laps = 0;
v_ctr = 0;
for k=1:total_num_laps
    tt_ctr = 1;
    % interpolate time for start of first bin
    iinds = find(SN_VT_TOT_POS{lap_numbers(k)} >= firing_rate_distribution{lap_numbers(k),7}(1));
    tt_ctr = iinds(1)-1;
    if(tt_ctr ~= 0) 
        tot_dist_diff = SN_VT_TOT_POS{lap_numbers(k)}(tt_ctr+1) - SN_VT_TOT_POS{lap_numbers(k)}(tt_ctr);
        this_dist_diff = firing_rate_distribution{lap_numbers(k),7}(1) - SN_VT_TOT_POS{lap_numbers(k)}(tt_ctr);
        frac_1 = this_dist_diff/tot_dist_diff;
        frac_2 = 1 - frac_1;
        firing_rate_distribution{lap_numbers(k),5}(1) = frac_1 * SMOOTH_VT_data{lap_numbers(k),1}(tt_ctr+1) + frac_2 * SMOOTH_VT_data{lap_numbers(k),1}(tt_ctr);
        v_ctr = v_ctr + 1;
        v_lap_numbers(v_ctr) = lap_numbers(k);
        num_v_laps = num_v_laps+1;
    else % bad VT data at beginning of place field
        fprintf('for lap %4.1f, bad early VT data\n',lap_numbers(k));
        for i=1:num_incrs
            firing_rate_distribution{lap_numbers(k),5}(i) = 0;
            firing_rate_distribution{lap_numbers(k),6}(i) = 0;
            firing_rate_distribution{lap_numbers(k),3}(i) = 0;
            firing_rate_distribution{lap_numbers(k),2}(i) = 0;
        end
        continue;
    end
        
    for i=1:num_incrs
        % interpolate time for end of each bin
        iinds = find(SN_VT_TOT_POS{lap_numbers(k)} >= firing_rate_distribution{lap_numbers(k),8}(i));

        if(isempty(iinds)) % bad VT data at end of place field
            fprintf('for lap %4.1f, bad late VT data\n',lap_numbers(k));    
            v_ctr = v_ctr - 1;
            num_v_laps = num_v_laps-1;
            for i=1:num_incrs
                firing_rate_distribution{lap_numbers(k),5}(i) = 0;
                firing_rate_distribution{lap_numbers(k),6}(i) = 0;
                firing_rate_distribution{lap_numbers(k),3}(i) = 0;
                firing_rate_distribution{lap_numbers(k),2}(i) = 0;
            end
            break;
        else
            tt_ctr = iinds(1);
            tot_dist_diff = SN_VT_TOT_POS{lap_numbers(k)}(tt_ctr) - SN_VT_TOT_POS{lap_numbers(k)}(tt_ctr-1);
            this_dist_diff = firing_rate_distribution{lap_numbers(k),8}(i) - SN_VT_TOT_POS{lap_numbers(k)}(tt_ctr-1);
            frac_1 = this_dist_diff/tot_dist_diff;
            frac_2 = 1 - frac_1;
            firing_rate_distribution{lap_numbers(k),6}(i) = frac_1 * SMOOTH_VT_data{lap_numbers(k),1}(tt_ctr) + frac_2 * SMOOTH_VT_data{lap_numbers(k),1}(tt_ctr-1);
    
            % enter time in each bin in seconds
            firing_rate_distribution{lap_numbers(k),3}(i) = (firing_rate_distribution{lap_numbers(k),6}(i) - firing_rate_distribution{lap_numbers(k),5}(i))/1000000;

            % compute firing rate in each bin - firing_rate_distribution{lap#,2}
            firing_rate_distribution{lap_numbers(k),2}(i) = firing_rate_distribution{lap_numbers(k),4}(i)/firing_rate_distribution{lap_numbers(k),3}(i);
        
            % enter start time for next bin
            if(i ~= num_incrs)
                firing_rate_distribution{lap_numbers(k),5}(i+1) = firing_rate_distribution{lap_numbers(k),6}(i);
            end
        end    
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute place field size = sum of instantaneous firing rate within each
% pixel = sum(firing_rate_distribution{lap #,2})
% units = cm * Hz

for k=1:total_num_laps
    place_field_size(lap_numbers(k)) = sum(firing_rate_distribution{lap_numbers(k),2});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute cell center of mass
% sum up all time spent in each spatial bin
% sum up all spikes fired in each spatial bin
% get overall firing rate in each bin - cell_firing_rate(1:incr)
% compute center of mass of cell_firing_rate(1:incr) over spatial bins
% then compute variance and skewness relative to cell_center_of_mass

for i=1:num_incrs
    cell_time_in_bin(i) = 0;
    cell_spikes_in_bin(i) = 0;
    cell_firing_rate(i) = 0;
    for k=1:num_v_laps
        cell_time_in_bin(i) = cell_time_in_bin(i) + firing_rate_distribution{v_lap_numbers(k),3}(i);
        cell_spikes_in_bin(i) = cell_spikes_in_bin(i) + firing_rate_distribution{v_lap_numbers(k),4}(i);
    end
    cell_firing_rate(i) = cell_spikes_in_bin(i)/cell_time_in_bin(i);
end
% compute cell center of mass
cell_cmtemp = [];
cell_cmtemp = firing_rate_distribution{v_lap_numbers(1),1}.*cell_firing_rate;
cell_mass = trapz(firing_rate_distribution{v_lap_numbers(1),1},cell_firing_rate);
cell_center_of_mass = trapz(firing_rate_distribution{v_lap_numbers(1),1},cell_cmtemp)/cell_mass;
% compute cell std
cell_cmtemp = firing_rate_distribution{v_lap_numbers(1),1}.*cell_cmtemp;
cell_raw_2nd_moment = trapz(firing_rate_distribution{v_lap_numbers(1),1},cell_cmtemp)/cell_mass;
cell_std = sqrt(cell_raw_2nd_moment - (cell_center_of_mass)^2);
% compute cell skewness
cell_cmtemp = firing_rate_distribution{v_lap_numbers(1),1}.*cell_cmtemp;
cell_raw_3rd_moment = trapz(firing_rate_distribution{v_lap_numbers(1),1},cell_cmtemp)/cell_mass;
cell_skewness = (cell_raw_3rd_moment - 3*cell_raw_2nd_moment*cell_center_of_mass + 2*(cell_center_of_mass)^3)/(cell_std)^3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute lap center of mass of firing rate distribution, 
% std deviation or variance relative to lap center of mass and
% skewness relative to lap center of mass
% and difference between lap center of mass and cell center of mass

for k=1:total_num_laps
    cmtemp = [];
    % do only for laps with more than 2 spikes and where firing rate
    % distribution was computed above - should be only valid laps
    if (sum(firing_rate_distribution{lap_numbers(k),2}) ~= 0 && spikes_in_lap(lap_numbers(k)) > 2)
        mass = trapz(firing_rate_distribution{lap_numbers(k),1},firing_rate_distribution{lap_numbers(k),2});
        cmtemp = firing_rate_distribution{lap_numbers(k),1}.*firing_rate_distribution{lap_numbers(k),2};
        center_of_mass(lap_numbers(k)) = trapz(firing_rate_distribution{lap_numbers(k),1},cmtemp)/mass;
        cmtemp = cmtemp.*firing_rate_distribution{lap_numbers(k),1};
        raw_second_moment = trapz(firing_rate_distribution{lap_numbers(k),1},cmtemp)/mass;
        std_deviation(lap_numbers(k)) = sqrt(raw_second_moment - center_of_mass(lap_numbers(k))^2);
        cmtemp = cmtemp.*firing_rate_distribution{lap_numbers(k),1};
        raw_third_moment = trapz(firing_rate_distribution{lap_numbers(k),1},cmtemp)/mass;
        skewness(lap_numbers(k)) = (raw_third_moment - 3*raw_second_moment*center_of_mass(lap_numbers(k)) + 2*(center_of_mass(lap_numbers(k)))^3)/std_deviation(lap_numbers(k))^3;
    else
        center_of_mass(lap_numbers(k)) = 0;
        std_deviation(lap_numbers(k)) = 0;
        skewness(lap_numbers(k)) = 0;
    end
    centers_wrt_cell_com(lap_numbers(k)) = center_of_mass(lap_numbers(k)) - cell_center_of_mass;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine peak of firing rate distribution

for k=1:total_num_laps
    if (sum(firing_rate_distribution{lap_numbers(k),2}) ~= 0 && spikes_in_lap(lap_numbers(k)) > 2)
        [peak_firing_rate,indx_of_peak] = max(firing_rate_distribution{lap_numbers(k),2});
        distribution_peak(lap_numbers(k)) = firing_rate_distribution{lap_numbers(k),1}(indx_of_peak);
    else
        distribution_peak(lap_numbers(k)) = 0;
    end
end

for k=1:total_num_laps
    dist_peak_wrt_cell_com(lap_numbers(k)) = distribution_peak(lap_numbers(k)) - cell_center_of_mass;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find locations of first and last spikes in each lap

for k=1:total_num_laps
    if (sum(firing_rate_distribution{lap_numbers(k),2}) ~= 0 && spikes_in_lap(lap_numbers(k)) > 2)
        inds_of_nz = find(firing_rate_distribution{lap_numbers(k),4});
        loc_first_spike_in_lap(lap_numbers(k)) = firing_rate_distribution{lap_numbers(k),1}(inds_of_nz(1));
        loc_last_spike_in_lap(lap_numbers(k)) = firing_rate_distribution{lap_numbers(k),1}(inds_of_nz(length(inds_of_nz)));
    else
        loc_first_spike_in_lap(lap_numbers(k)) = 0;
        loc_last_spike_in_lap(lap_numbers(k)) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plot of spike location for each lap

convert_pxls_to_cm = TOTAL_TRK_LEN/TRK_LEN_PXLS;
figure(2);
clf;
hold on;
for k=1:total_num_laps
    tempy = ones(length(place_cell_sn_pos_tot{lap_numbers(k)}))*lap_numbers(k);
    plot(place_cell_sn_pos_tot{lap_numbers(k)}*convert_pxls_to_cm,tempy,'ko');
    if(sum(firing_rate_distribution{lap_numbers(k),2}) ~= 0 && spikes_in_lap(lap_numbers(k)) > 2)
        plot(center_of_mass(lap_numbers(k))*convert_pxls_to_cm,lap_numbers(k),'r*')
    end
    plot([cell_center_of_mass*convert_pxls_to_cm cell_center_of_mass*convert_pxls_to_cm],[lap_numbers(1) lap_numbers(total_num_laps)],'k--');
    xlabel('Position (cm)'), ylabel('Lap number');
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot firing rate distributions
figure(3);
clf;
hold on;
for k=1:total_num_laps
    if(sum(firing_rate_distribution{lap_numbers(k),2}) ~= 0 && spikes_in_lap(lap_numbers(k)) > 2)
        plotnum = lap_numbers(k);
        subplot(5,3,plotnum);
        hold on;
        plot(firing_rate_distribution{lap_numbers(k),1}*convert_pxls_to_cm, firing_rate_distribution{lap_numbers(k),2});
        plot([center_of_mass(lap_numbers(k))*convert_pxls_to_cm center_of_mass(lap_numbers(k))*convert_pxls_to_cm],...
            [0 100]);
        plot([cell_center_of_mass*convert_pxls_to_cm cell_center_of_mass*convert_pxls_to_cm], [0 100], 'k--');
        hold off;
        xlim([firing_rate_distribution{lap_numbers(k),1}(1)*convert_pxls_to_cm ...
                firing_rate_distribution{lap_numbers(k),1}(num_incrs)*convert_pxls_to_cm]);
        title(['Lap ',num2str(lap_numbers(k))]);
        
    end
end
hold off;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot cell firing rate distribution
figure(4);
clf;
hold on;
plot(firing_rate_distribution{lap_numbers(1),1}*convert_pxls_to_cm, cell_firing_rate);
plot([cell_center_of_mass*convert_pxls_to_cm cell_center_of_mass*convert_pxls_to_cm], [0 max(cell_firing_rate)], 'k--');
xlabel('Position (cm)'), ylabel('Firing rate (Hz)');
title('Cell Firing Rate Distribution');
axis tight;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write output

cd(DATADIR);
[filename, pathname] = uiputfile({'*.csv',...
      'Firing Rate Map File (*.csv)'},sprintf('Save Firing Rate Map/Distribution Info'));

        cd(working_dir);
        output_file_name= fullfile(pathname, filename);
        
        OUTPUT_FILE = fopen(output_file_name,'w');
        
        fprintf(OUTPUT_FILE,'This is firing rate map/distribution info.\n');
        fprintf(OUTPUT_FILE,'Conversion from pixels to cm %12.7f \n', convert_pxls_to_cm);
        fprintf(OUTPUT_FILE,'Total number of laps  %6.0f \n',total_num_laps);
        fprintf(OUTPUT_FILE,'place field cell center of mass (in pixels), (in cm)  %12.7f %12.7f \n',cell_center_of_mass, cell_center_of_mass*convert_pxls_to_cm);
        fprintf(OUTPUT_FILE,'place field cell standard deviation (in pixels), (in cm)  %12.7f %12.7f \n',cell_std, cell_std*convert_pxls_to_cm);
        fprintf(OUTPUT_FILE,'place field cell skewness  %12.7f \n',cell_skewness);
        fprintf(OUTPUT_FILE,'lap number, place field size, place field center (pxls), place field center (cm), lap center relative to cell center (pxls), lap center relative to cell center (cm), place field peak (pxls), place field peak (cm), peak relative to cell center (pxls), peak relative to cell center (cm), location of first spike (pxls), location of first spike (cm), location of last spike (pxls), location of last spike (cm), place field width (pxls), place field width (cm), std deviation (pxls), std deviation (cm), skewness \n');
        for k=1:total_num_laps
            place_field_width = loc_last_spike_in_lap(lap_numbers(k))-loc_first_spike_in_lap(lap_numbers(k));
            fprintf(OUTPUT_FILE,'%12.7f, ',lap_numbers(k));
            fprintf(OUTPUT_FILE,'%12.7f, ',place_field_size(lap_numbers(k)));
            fprintf(OUTPUT_FILE,'%12.7f, ',center_of_mass(lap_numbers(k)));
            fprintf(OUTPUT_FILE,'%12.7f, ',center_of_mass(lap_numbers(k))*convert_pxls_to_cm);
            fprintf(OUTPUT_FILE,'%12.7f, ',centers_wrt_cell_com(lap_numbers(k)));
            fprintf(OUTPUT_FILE,'%12.7f, ',centers_wrt_cell_com(lap_numbers(k))*convert_pxls_to_cm);
            fprintf(OUTPUT_FILE,'%12.7f, ',distribution_peak(lap_numbers(k)));
            fprintf(OUTPUT_FILE,'%12.7f, ',distribution_peak(lap_numbers(k))*convert_pxls_to_cm);
            fprintf(OUTPUT_FILE,'%12.7f, ',dist_peak_wrt_cell_com(lap_numbers(k)));
            fprintf(OUTPUT_FILE,'%12.7f, ',dist_peak_wrt_cell_com(lap_numbers(k))*convert_pxls_to_cm);
            fprintf(OUTPUT_FILE,'%12.7f, ',loc_first_spike_in_lap(lap_numbers(k)));
            fprintf(OUTPUT_FILE,'%12.7f, ',loc_first_spike_in_lap(lap_numbers(k))*convert_pxls_to_cm);
            fprintf(OUTPUT_FILE,'%12.7f, ',loc_last_spike_in_lap(lap_numbers(k)));
            fprintf(OUTPUT_FILE,'%12.7f, ',loc_last_spike_in_lap(lap_numbers(k))*convert_pxls_to_cm);
            fprintf(OUTPUT_FILE,'%12.7f, ',place_field_width);
            fprintf(OUTPUT_FILE,'%12.7f, ',place_field_width*convert_pxls_to_cm);
            fprintf(OUTPUT_FILE,'%12.7f, ',std_deviation(lap_numbers(k)));
            fprintf(OUTPUT_FILE,'%12.7f, ',std_deviation(lap_numbers(k))*convert_pxls_to_cm);
            fprintf(OUTPUT_FILE,'%12.7f\n',skewness(lap_numbers(k)));
        end
        
        fprintf(OUTPUT_FILE,'Cell Firing Rate Distribution \n');
        fprintf(OUTPUT_FILE,'location (pxls), location (cm), firing rate in location, time in location, number of spikes in location\n');
        for i=1:length(firing_rate_distribution{lap_numbers(1),1})
            fprintf(OUTPUT_FILE,'%12.7f, ',firing_rate_distribution{lap_numbers(1),1}(i));
            fprintf(OUTPUT_FILE,'%12.7f, ',firing_rate_distribution{lap_numbers(1),1}(i)*convert_pxls_to_cm);
            fprintf(OUTPUT_FILE,'%12.7f, ',cell_firing_rate(i));
            fprintf(OUTPUT_FILE,'%12.7f, ',cell_time_in_bin(i));
            fprintf(OUTPUT_FILE,'%12.7f\n',cell_spikes_in_bin(i));
        end
        
        fprintf(OUTPUT_FILE,'Firing Rate Distribution data for each lap \n');
        fprintf(OUTPUT_FILE,'lap number \n');
        fprintf(OUTPUT_FILE,'location (pxls), location (cm), firing rate in location, time in location (sec), number of spikes in location, start time in location, end time in location, start position of location, end position of location \n');
        
        for k=1:total_num_laps
            fprintf(OUTPUT_FILE,'lap number  %6.0f \n',lap_numbers(k));
            for i=1:length(firing_rate_distribution{lap_numbers(k),1})
               fprintf(OUTPUT_FILE,'%12.7f, ',firing_rate_distribution{lap_numbers(k),1}(i));
               fprintf(OUTPUT_FILE,'%12.7f, ',firing_rate_distribution{lap_numbers(k),1}(i)*convert_pxls_to_cm);
               fprintf(OUTPUT_FILE,'%12.7f, ',firing_rate_distribution{lap_numbers(k),2}(i));
               fprintf(OUTPUT_FILE,'%12.7f, ',firing_rate_distribution{lap_numbers(k),3}(i));
               fprintf(OUTPUT_FILE,'%12.7f, ',firing_rate_distribution{lap_numbers(k),4}(i));
               fprintf(OUTPUT_FILE,'%12.7f, ',firing_rate_distribution{lap_numbers(k),5}(i));
               fprintf(OUTPUT_FILE,'%12.7f, ',firing_rate_distribution{lap_numbers(k),6}(i));
               fprintf(OUTPUT_FILE,'%12.7f, ',firing_rate_distribution{lap_numbers(k),7}(i));
               fprintf(OUTPUT_FILE,'%12.7f\n',firing_rate_distribution{lap_numbers(k),8}(i));

           end
        end
        
        fclose(OUTPUT_FILE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write spike position data

cd(DATADIR);
[filename, pathname] = uiputfile({'*.csv',...
      'Spike Position File (*.csv)'},sprintf('Save Spike Position Info'));

        cd(working_dir);
        output_file_name= fullfile(pathname, filename);
        
        OUTPUT_FILE = fopen(output_file_name,'w');
        
        fprintf(OUTPUT_FILE,'This is spike position info for each lap.\n');
        fprintf(OUTPUT_FILE,'Interpolated positions from interpolated and smooth VT data.\n');
        fprintf(OUTPUT_FILE,'lap number\n');
        fprintf(OUTPUT_FILE,'timestamp, position on straightened track (pxls), position on straightened track (cm) \n');
        for k=1:total_num_laps
            fprintf(OUTPUT_FILE,'lap number  %6.0f \n',lap_numbers(k));
            for i=1:spikes_in_lap(lap_numbers(k))
                fprintf(OUTPUT_FILE,'%12.7f, ',place_cell_firing{lap_numbers(k),1}(i));
                fprintf(OUTPUT_FILE,'%12.7f, ',place_cell_sn_pos_tot{lap_numbers(k)}(i));
                fprintf(OUTPUT_FILE,'%12.7f, ',place_cell_sn_pos_tot{lap_numbers(k)}(i)*convert_pxls_to_cm);
            end
        end
        
fclose(OUTPUT_FILE);  











               


