function [features, values_m, values_a] = MARA_extract_time_freq_features(EEG)
% function [features, values_m, values_a] = MARA_extract_time_freq_features(EEG)    
%      
%   values_m = avg. resid value btw 20 and 50 Hz
%   values_a = avg. resid value btw 8 and 13 Hz
%   features = list of components with bumps
%                           
data = EEG.data;
fs = EEG.srate;

% transform epoched data into continous data
if length(size(data)) == 3
    s = size(data); 
    data = reshape(data, [EEG.nbchan, prod(s(2:3))]); 
end

%downsample (to 100-200Hz) 
factor = max(floor(EEG.srate/100),1); 
data = data(:, 1:factor:end); 
fs = round(fs/factor); 

%compute icaactivation and standardise variance to 1
icacomps = (EEG.icaweights * EEG.icasphere * data)';
icacomps = icacomps./repmat(std(icacomps,0,1),length(icacomps(:,1)),1);
icacomps = icacomps';

%creating a matrix to help us check for peaks in data
change=[0 linspace(-1,1,50)];
components_with_bumps = [];
components_with_bumps_2 = [];
components_with_plats = [];
components_with_plats_2 = [];
pxx_matrix = zeros(length(icacomps(:,1)),51);
resids_pxx_matrix = zeros(length(icacomps(:,1)),50);

%run through all the components
for ic=1:length(icacomps(:,1)) 
    fprintf('.');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Proc Spectrum for Channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pxx,freq,speccomp,contrib,specstd] = spectopo(icacomps(ic,:), fs, fs);
    pxx = 10*log10(pxx * 100/2);
    pxx_matrix(ic,1:51)=pxx(1:51);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get the residual values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pxx = real(pxx); %need for spectopo
    pxx = pxx*-1; %need for spectopo
    OneOverF = (1./freq);
    lm_pxx = fitlm(OneOverF(2:51),pxx(2:51));
    resids_pxx = lm_pxx.Residuals.Standardized;
    resids_pxx_matrix(ic,1:50)=resids_pxx;
    
    % at one point we were looking at power to try and identify muscle acitivity, but we have not
    % found it as useful as the peak finder
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % The average log band power between 8 and 13 Hz
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     p = 0;
%     for i = 8:13 
%         p = p + resids_pxx(i);
%     end
%     Hz8_13 = p / (13-8+1);
%     values_a(1,ic) = Hz8_13;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % The average residual of the log band power between 20 and 50 Hz
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     p = 0;
%     for fv = 20:49 
%         p = p + resids_pxx(fv);
%     end
%     Hz20_50 = p / (50-20+1);
%     values_m(1,ic) = Hz20_50;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Checking for bumps between 5 and 15 Hz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % look for peaks with the residual values
    [pks, locs, widths, proms] = findpeaks(resids_pxx(4:16),[4:16]); %looks for peaks btw 5 and 15 Hz (1 less and more so any bumps at 5 or 15 are caught)
    for numPeaks = 1:length(pks) %in case we have mult. peaks
        if (proms(numPeaks) > 0.3 && pks(numPeaks)>1 && widths(numPeaks) > 0.9) %conservative
%         if (proms(numPeaks) > 0.5 && pks(numPeaks)>1 && widths(numPeaks) > 0.9) %liberal
            components_with_bumps = [components_with_bumps ic];
        end
    end
    
    %%%%%%%%%%%%%%%%
    % We are not sure if this chunk without the residuals is needed... 
    % the code above just seemed too cautious
    % We are testing it now and might make changes
    
    % look for peaks with the NON-residual values
    % much more liberal threshold
    % above code using the residual values sometimes gets tiny plateaus or noise
    % the below code doesn't catch that stuff
    [pks, locs, widths, proms] = findpeaks(pxx(5:17),[5:17]); %looks for peaks btw 5 and 15 Hz (1 less and more so any bumps at 5 or 15 are caught)
    for numPeaks = 1:length(pks) %in case we have mult. peaks
        if (proms(numPeaks) > 0.15) %conservative
%         if (proms(numPeaks) > 0.35) %liberal
            components_with_bumps_2 = [components_with_bumps_2 ic];
        end
    end
    %%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Checking for plateaus between 5 and 15 Hz for first 7 components
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ic < 8
        [pks, locs, widths, proms] = findpeaks(resids_pxx(4:16),[4:16]);
        for numPeaks = 1:length(pks) %in case we have mult. peaks
            if (proms(numPeaks) > 0.15 && pks(numPeaks)>1) %conservative
%             if (proms(numPeaks) > 0.15 && pks(numPeaks)>1) %liberal
                components_with_plats = [components_with_plats ic];
            end
        end
        
        %%%%%%%%%%%%%%%%
        % We are not sure if this chunk without the residuals is needed... 
        % the code above just seemed too cautious
        % We are testing it now and might make changes
        [pks, locs, widths, proms] = findpeaks(pxx(5:17)+(10*change(5:17)),[5:17]);
        for numPeaks = 1:length(pks) %in case we have mult. peaks
            if (proms(numPeaks) > 0.05) %consevative
%             if (proms(numPeaks) > 0.05) %liberal
                components_with_plats_2 = [components_with_plats_2 ic];
            end
        end
        %%%%%%%%%%%%%%%
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Append Features 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %features(:,ic)= [Hz8_13]; 
end

% intialize these for later
features = []; % this will be our output matrix (our list of components with bumps)
features_2 = [];
% combine the two components with bumps matrices
components_with_bumps = unique(components_with_bumps);
components_with_bumps_2 = unique(components_with_bumps_2);
for check=1:length(components_with_bumps)
    if sum(components_with_bumps(check) == components_with_bumps_2)
        % if they both have a component, add it to the features matrix
        features = [features components_with_bumps(check)];
    end
end
% combine the two components with plateaus matrices
components_with_plats = unique(components_with_plats);
components_with_plats_2 = unique(components_with_plats_2);
for check=1:length(components_with_plats)
    if sum(components_with_plats(check) == components_with_plats_2)
        features_2 = [features_2 components_with_plats(check)];
    end
end
% combine the components with bumps matrix and the plateau matrix
features = unique([features features_2]); %output matrix (list of components with bumps)

% plot componetns with bumps and components w/o bumps
figure; hold on;
bumps=subplot(2,1,1);
for jj=1:length(features)
    plot(bumps,1:50,resids_pxx_matrix(features(jj),:))
    hold on
end
title(bumps,'ICs with Bumps')

good_ics=1:EEG.nbchan;
good_ics(features)=[];
nobumps=subplot(2,1,2);
for gg=1:length(good_ics)
    plot(nobumps,1:50,resids_pxx_matrix(good_ics(gg),:))
    hold on
end
title(nobumps,'ICs without Bumps')

par_id=strsplit(EEG.setname,'_')
saveas(gcf,[par_id{1} '.jpg'])

%open a new figure for the spectopo output
hold off; figure;

% features = unique(components_with_bumps);
disp('.');
end