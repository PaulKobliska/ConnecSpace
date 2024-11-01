function [EEGvalues, EEGts, TTMtx, PosMtx, E_FlagMtx, eeg_chan_to_read,inputfile] = Read_UFF_Slow(inputfile, read_spikes, read_pos, read_eeg,  read_events,eeg_chan_to_read);
%{
Function that reads different types of data from an uff file (spike, 
eeg, position, event flags...). Returns matrices of values.
input:    inputfile:  name (with or without the whole path) shoudl contain .uff extension!
              (if only the name of the file is given, will look for it
              in the current directory)
          read_spikes: 1 if you want to read the spike data, 0 otherwise
          read_pos: 1 if we want to read the spike data, 0 otherwise
          read_eeg: 1 if we want to read the spike data, 0 otherwise
          read_events, 1 if we want to read the spike data, 0 otherwise
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
          eeg_chan_to_read: the channel from which LFP data is extracted. If unspecified, will choose the channel 
              with lowest number available in the data file. if 0, it will
              choose all channels
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
output:   EEGvalues:  vector of values (voltage) of the continuous data (EEG)
          EEGts:  vector of timesteps corresponding to EEGvalues
          TTMtx:  matrix with the time of each spike and its waveform
                col1: ts (in ms); col2: probe; col3: channel
                from col 4 to col 131: waveform amplitude (32 col for each
                channel) organized as follolws: the first 4 (from col 4 to
                col 7) contains values of the first timepoint (out of 32)
                for all channels (col4 ch 1, col5 ch2, col6 ch3, col7 ch4),
                the second 4 (from col 8 to col 11) contains values of the 
                second timepoint (out of 32) for all channels (col5 ch 1, 
                col6 ch2, col7 ch3, col8 ch4), and so on...              
          PosMtx: matrix of timesteps (in sec, col 1) and x (col 2) and y
                   (col 3) corresponding positions
          E_FlagMtx:  matrix of event flags corresponding to each 
              timestep
          eeg_chan_to_read: integer corresponding to the number of the
              actually read channel (generally, 16)


june 2016 modified by Francesca
dec 2013: Modified by L.L. Bologna, E. Duvelle
jan 2013: Modified by E. Duvelle elduvelle@free.fr
modified by PP.

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable to define:
max_chan_num = 32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Length of each block of data in the file (in octets?)
Offset.U = 282; % quatrode Offset (2 more bytes at the end)
Offset.S= 152;  % Stereotrode
Offset.E= 90;   % single Electrode
Offset.P= 12;   % Position
Offset.B = 10;  % Event Flag
Offset.M = 88;  % Message
Offset.R = 1094; % continuous Record
Offset.T = 9008;% Trial descriptor record
Offset.L= 3272; % spike Layout
Offset.Header = 2050;
Offset.DRDB = 2050*9;

% Open the file
fid=fopen(inputfile,'r','l');


%%%%%%%%%%%%%%%%%%%%%
% Data initialization

% Get data dimensions by reading the data file a first time
[eeg_rows_per_channel, Tetrodrows,Tetrodcols,Posrows,Poscols, E_Flagrows, E_Flagcols]= Get_All_MxSize_UFF(fid, Offset,eeg_chan_to_read,max_chan_num);

if read_eeg
    % EEG : Get the channels which contain data
    eeg_channels = find(eeg_rows_per_channel ~= 0);
    % Choose the first channel of the existing channels by default 
    if nargin <6 
        if numel(eeg_channels)~=0
            eeg_chan_to_read = eeg_channels(1);
            txt = ['Setting EEG channel to be read to ' num2str(eeg_chan_to_read)];
            disp(txt)
        end
    end

    if numel(eeg_channels)==0
        disp('Warning! No continuous data found.')
        read_eeg = 0;
        % if no eeg available
        eeg_chan_to_read = -1;
        EEGvalues = [];
        EEGts = [];
    else
        % Initialization of data arrays
%         EEGvalues = zeros(eeg_rows_per_channel(eeg_chan_to_read),max_chan_num);
%         EEGts = zeros(eeg_rows_per_channel(eeg_chan_to_read),max_chan_num); 
          EEGvalues = zeros(eeg_rows_per_channel(1),max_chan_num);
          EEGts = zeros(eeg_rows_per_channel(1),max_chan_num);
    end
    

else
    eeg_chan_to_read = -1;
    EEGvalues = [];
    EEGts = [];
end


if read_spikes
    % Spikes (+ 1 for ts, 1 for probe, 1 for clust = 3)
    %    OLD : FIND THE GAINS!!!!!
    TTMtx = zeros(Tetrodrows,Tetrodcols+3); 
else
    TTMtx = [];
end

if read_pos
    % Position (+1 for the Timestep column)
    PosMtx = zeros(Posrows, Poscols+1); 
else
    PosMtx = [];
end


if read_events
    % Event-flags
    E_FlagMtx = zeros(E_Flagrows, E_Flagcols);
else
    E_FlagMtx = [];
end




%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parsing of the data file


% Skip the file header
fseek (fid, Offset.Header, 'bof');
% Skip the file DRDB (don't really know what it is)
fseek (fid, Offset.DRDB, 'cof');
% create list of increment for each channel (used to add values to  EEGvalues)
EEG_index = ones(max_chan_num,1);
i=1;
j=0;
k=0;
l=0;

while feof(fid) == 0
    
    % Read the character at the beginning of a new block of data.
    % The letter indicates which type of data is stored in this block.
    % If the block doesn't interest us, we just skip to the next block.
    ident2 = fread (fid,1,'uchar');
    valident=char(ident2');

     
    switch valident
        
        case ('U')
            if read_spikes
                % Unit data (spikes, waveforms)

                l=l+1;

                % Probe
                TTMtx(l,2)=fread (fid,1,'uint8');
                % Time : convert it to seconds 
                ts = fread (fid,1,'int32');
                TTMtx(l,1) = ts/10; % time in ms
                % Cluster?
                TTMtx(l,3)=(fread (fid,1,'int16'))';
                fseek (fid, 16, 'cof');
                % OLD    to see them : go 4 by 4
                % TODO: Check this
                TTMtx(l,4:131)=(fread (fid,128,'int16'))';   
                fseek (fid, 2, 'cof');
            else
                % Skip the spikes part
                fseek(fid, Offset.U -1, 'cof');
            end
            
        case ('S')
            fseek (fid, Offset.S -1, 'cof');
        case ('E')
            fseek (fid, Offset.E -1, 'cof');
        case ('P')
            if read_pos
                
                % Position data
                j=j+1;
                fseek (fid,1,'cof');
                % Read the time
                ts = fread (fid,1,'int32');
                % Change the time into mseconds
                PosMtx(j,1) = ts/10;
                % Read the position
                PosMtx(j,2:5)=(fread (fid,4,'uint8'))';
                %fseek (fid, 2,'cof');
                fseek (fid, 2,'cof');
                
            else
                % Skip the position part
                fseek(fid, Offset.P-1, 'cof');
            end
        case ('B')
            if read_events
                % Event-flag data
                % todo get proper characters
                k = k+1;
                fseek (fid,1,'cof');

                % Read the time
                ts = fread(fid,1,'int32');
                % Change the time to mseconds
                E_FlagMtx(k,1) = ts/10;

                % Read the event flags (2 characters)
                eflag_1 = fread(fid, 1, 'uchar');
                eflag_1 = char(eflag_1);
                eflag_2 = fread(fid, 1, 'uchar');
                eflag_2 = char(eflag_2);
                E_FlagMtx(k,2) = eflag_1;
                E_FlagMtx(k,3) = eflag_2;


                fseek (fid, 2,'cof');
            else
                % Skip the event-flags data
                fseek(fid, Offset.B -1, 'cof');
            end
        case ('M')
            fseek (fid, Offset.M -1,'cof');
        case ('R')
            if read_eeg
                % Continuous data
                % read the actual channel
                EEGchan = fread (fid, 1, 'uint8');
                EEGchan = EEGchan+1;
                % ED. Only read if data is from the chosen channel
                 
%                 if EEGchan == eeg_chan_to_read


                    % Read the timestamp (Note: it will be changed into 
                    % seconds afterwards)
                    ts = fread(fid, 1, 'int32');

                    % Read the burst seq buf num (??)
                    brst_freq = fread(fid,1, 'int16');

                    % Read the number of valid samples 
                    % Note: continuous data is stored by chunks of 512 but
                    % one time out of two, only a portion of the data was
                    % stored and the remaining elements up to 512 are 0.
                    % to override this problem we need to take into account
                    % the number of valid samples stored here.
                    num_spl = fread(fid,1, 'int16');
                    %num_spl = 512;

                    % Read the sampling frequency (useful for the timestep 
                    % array)
                    spl_freq = fread(fid,1, 'int32');

                    % Read the volts per AD unit (don't know)
                    volts_per_ad = fread(fid,1, 'int32');

                    % Skip useless (?) data 
                    fseek(fid, 50, 'cof');

                    % Write the EEG values in the EEG values array
                    EEGvalues((EEG_index(EEGchan) : num_spl+(EEG_index(EEGchan))-1),EEGchan) = fread (fid,num_spl,'int16');         

                    % Write the the timestamps in the EEG timestamps array 
                    % Divide the time by 10 to get mseconds. 
                    EEGts((EEG_index(EEGchan):num_spl + (EEG_index(EEGchan)) - 1),EEGchan) = (ts/10 + ((0:num_spl-1)*(1/spl_freq)));

                    % Advance the pointer in case there were zeros in the data
                    fseek(fid, (512 -num_spl)*2, 'cof');

                    % Advance to the next chunk of data
                    fseek (fid, 2, 'cof');

                    % increase the counter for next reading
                    EEG_index(EEGchan) = EEG_index(EEGchan) + num_spl;                
                %end

            else
                % Skip the eeg data
                fseek (fid, Offset.R -1, 'cof');
            end   
            
        case ('T')
            fseek (fid, Offset.T -1, 'cof');
        case ('L')
            fseek (fid, Offset.L -1, 'cof');
            
        otherwise
            if feof(fid) == 1
                bye = 'byebye';
            else
                error ('Probleme d offset majeur in main program')

            end
   
    end

end

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% I added this part to have either 1 channel or all channels

if eeg_chan_to_read
    if read_eeg
        EEGvalues=EEGvalues(:,eeg_chan_to_read+1);
        EEGts=EEGts(:,eeg_chan_to_read+1);
    end
end

% Remove cluster 0 from TTMtx (noise). 
% if read_spikes
%     TTMtx(TTMtx(:,3)==0,:) = [];
% end

