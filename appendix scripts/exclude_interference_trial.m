%%
% This script marks interference trials in eeg experiment. See MADE paper 
% for details about interference trials.  

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

%% Read excel file containing trial numbers to be excluded from analysis because of interference during experiment
if subject ==1 % read the file once
    interference_trialfile_location = ''; % location of the file
    interference_trialfile_name = 'interference coding.xlsx'; % excel file name
    interference_trialfile_sheet = 'Sheet1'; % excel sheet to read data from
    interference_trialfile_rwcl = 'A2:C10'; % row and column numbers. starts from A2 to exclude header
    interference_trialfile_header = {'subjectID', 'Condition1', 'Condition2'}; % headers
    
    % Import the data
    [~, ~, raw] = xlsread([interference_trialfile_location, filesep, interference_trialfile_name], interference_trialfile_sheet, interference_trialfile_rwcl);
    interference_data = string(raw);
    interference_data(ismissing(interference_data)) = '';
end

%% Label events
% Add a new field 'TrialNumber' in EEGLAB event list and assign the value 'NaN'
EEG = pop_editeventfield(EEG, 'indices',  strcat('1:', int2str(length(EEG.event))), 'TrialNumber', 'NaN');
EEG = eeg_checkset(EEG);

% Add trial number for each condition
if ~isempty(task_event_markers)==1
    for tm = 1:length(task_event_markers)
        tn=1;
        for em=1:length(EEG.event)
            if strcmp(EEG.event(em).type, task_event_markers{tm})
                EEG.event(em).TrialNumber = tn;
                tn=tn+1;
            end
        end
    end
end

%% Mark interference trials
bad_trials=zeros(length(task_event_markers), length(EEG.event)); 
for subidx=1:length(interference_data)
    if  strcmp(interference_data{subidx}, (datafile_names{subject}(1:end-length(ext))))==1 % check the subjectID from the excle file
        for ii=1:length(task_event_markers)
            bad_trials(ii,(str2num(interference_data{subidx, ii+1})))=1;  % get bad trial list for the subject            
        end
        break;
    end    
end

% Mark trials
for te=1:length(task_event_markers)
    idx_bt=find(bad_trials(te,:)==1);
    for bt=1:length(idx_bt)
        for kk=1:length(EEG.event)
            if strcmp(EEG.event(kk).type, task_event_markers{te})==1 && EEG.event(kk).TrialNumber==idx_bt(bt) %find(bad_trials(te, bt)))
                EEG.event(kk).type = [EEG.event(kk).type '_bad'];
            end
        end
    end                                                                                                                                                                 
end
