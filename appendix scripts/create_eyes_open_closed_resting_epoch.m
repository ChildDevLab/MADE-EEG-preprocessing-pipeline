%% Often resting state EEG is collected in eyes close and eyes open conditions
% This script inserts event at specific time interval for resting state EEG
% and creates separate events for eyes close and eyes open resting EEG
% data. See MADE paper for details.

% Authors:
% Ranjan Debnath (rdebnatah@umd.edu)
% George A. Buzzell (gbuzzell@umd.edu)
% Santiago Morales  (moraless@umd.edu)
% Stephanie Leach (sleach12@umd.edu)
% Maureen E. Bowers (mbowers1@umd.edu)
% Nathan A. Fox (fox@umd.edu)

% Copyright (C) 2019 Ranjan Debnath, George A. Buzzell, Santiago Morales, Maureen E. Bowers, Nathan A. Fox
% Child Development Lab, University of Maryland, College Park, MD, USA.

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.

% ********************************************************************* %

%% Initialize variables

% 1. Name of eyes close and eyes open markers
rest_event_markers = {'NtCl', 'NTOp'}; % enter eyes close and eyes open markers
% 2. Name of new eyes close and eyes open markers
new_rest_markers={'new_eyes_close_marker', 'new_eyes_open_marker'}; % enter new name
% 3. Does the data have a trial end marker?
trial_end_marker=1; % 0=NO (no trial end marker in data), 1=YES (data have trial end marker)
trial_end_marker_name=('TrEd'); % enter trial end marker name
% 4. Length of epochs in second
rest_epoch_length=2;
% 5. Do you want to create overlapping epoch?
overlap_epoch = 1; % 0 = NO (do not create overlapping epoch), 1 = YES (50% overlapping epoch)

%% Insert markers
if overlap_epoch==1
    time_samples = (rest_epoch_length/2)*EEG.srate; % convert time window into samples or data points
else
    time_samples = rest_epoch_length*EEG.srate;
end
EEG.urevent=EEG.event;

tm=1;
for ue=1:length(EEG.urevent)
    for rm=1:length(rest_event_markers)
        if strcmp(EEG.urevent(ue).type, rest_event_markers{rm})==1
            if trial_end_marker ==1
                for te=ue:length(EEG.urevent)
                    if strcmp(EEG.urevent(te).type, trial_end_marker_name)==1
                        trial_end_latency=EEG.urevent(te).latency;
                        break;
                    end
                end
            else
                trial_end_latency=EEG.urevent(ue+1).latency;
            end
            event_times = EEG.urevent(ue).latency;
            while event_times <trial_end_latency
                EEG.event(tm).type=new_rest_markers{rm};
                EEG.event(tm).latency=event_times;
                event_times = event_times+time_samples;
                tm=tm+1;
            end
        end
    end
end

%% create epoch
EEG = eeg_checkset(EEG);
EEG = pop_epoch(EEG, new_rest_markers, [0 rest_epoch_length], 'epochinfo', 'yes');
EEG = eeg_checkset(EEG);

            