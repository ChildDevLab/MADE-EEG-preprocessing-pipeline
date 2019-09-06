%%%%% Marker editing script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script edits or inserts markers as needed, as well as adds relevant
% information to the markers (e.g. which segments are bad, task information, etc.
%
% This is only a simple example of marker labeling for a specific
% experiment. Typically, additional information would be labeled for each
% marker. It is also important to note that this script only serves as an
% example and will need to be edited to work with the specific experimental
% script that was used to collect your own data. See MADE paper for details.

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

%%
% Add labels to the event structure for: EventType, Congruency, Responded, Accuracy
EEG = pop_editeventfield( EEG, 'indices',  strcat('1:', int2str(length(EEG.event))), 'EventType','NaN','Congruency','NaN','Responded','Accuracy','NaN');
EEG = eeg_checkset( EEG );

% Find all TRSP markers (find the marker numbers that correspond to the TRSP
% markers. The TRSP markers contain information about what happened on a
% given trial, which we will use to label the stimulus and response markers
TrspEvents = find(strcmp({EEG.event.type}, 'TRSP'));

% Loop through all of the TRSP events and label the corresponding stimulus/response markers
% for a given trial according to the TRSP information
for t = TrspEvents %t = EVENT numbers associated with TRSPs (the event number for a TRSP marker is likely not the same as the trial number that the TRSP occurred on)
    
    % Figure out if a response was made on this trial. If so, identify the
    % event # corresponding to the response marker
    if strcmp(EEG.event(t).mffkey_RPTP,'Omission') % if no response
        % Store variable to indicate that there was NO response made on this trial
        Responded = 0;
        % Store variable to indicate there is no response marker (because no response was made)
        RespEventNum = 0;
    else
        % Store vriable to indicate that there was a response made on this trial
        Responded = 1;
        % Identify event # for marker indicating response onset for this trial.
        % When searching for the response marker number, we only search within
        % a range of marker numbers within which the response marker could have
        % occurred. We cannot simply pull the marker that is n markers back,
        % because the number of markers will vary on each trial depending on
        % whether a response was made or not.
        RespEventNum = t - (5 - find(strcmp({EEG.event(t-4:t-1).type}, 'resp')));
        % Assign the value "Resp" to the 'EventType' field
        EEG.event(RespEventNum).EventType = 'Resp';
    end % End conditional to check if response occurred on this trial or not
    
    % Identify event # for marker indicating stimulus onset for this trial.
    % When searching for the stimulus marker number, we only search within
    % a range of marker numbers within which the stimulus marker could have
    % occurred. We cannot simply pull the marker that is n markers back,
    % because the number of markers will vary on each trial depending on
    % whether a response was made or not.
    StimEventNum = t - (5 - find(strcmp({EEG.event(t-4:t-1).type}, 'FLAN')));
    %Assign the value "Stim" to the 'EventType' field
    EEG.event(StimEventNum).EventType = 'Stim';
    
    % Loop to label the stimulus and response markers (if present) for this trial
    for EventNum = [StimEventNum RespEventNum]
        % if trial had no response, RespEventNum will be zero, so we will
        % skip labelling the response marker on this trial (as there is no
        % response marker to label)
        if EventNum ~= 0
            % label congruency of the stimulus shown
            EEG.event(EventNum).Congruency = EEG.event(t).mffkey_CoNo;
            % label whether response occured
            EEG.event(EventNum).Responded = Responded;
            % label accuracy of the response
            EEG.event(EventNum).Accuracy = EEG.event(t).mffkey_eval;
        end % end eventNum=0 conditional
    end % end labelling loop for this trial
end % end loop through TRSP events (all trials)
