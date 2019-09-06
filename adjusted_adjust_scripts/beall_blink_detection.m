% Author: Daniel J. Beall
% Email: DJBEALL1101@gmail.com
%
% beall_blink_detection() - Returns list of ICA components marked as a blink 
%
% Usage:
%   >> [blinks]=beall_blink_detection(cwb,topog,chanlocs,n);
%
% Inputs:
%   cwb        - list of components with bumps
%   topog      - topographies vector
%   chanlocs   - EEG.chanlocs struct
%   n          - number of ICs
%   nchannels  - number of channels
% Outputs:
%   blinks     - list of detected blinks
%

function [blinks] = beall_blink_detection(cwb,topog,chanlocs,n)

blinks = [];
nchannels=length(chanlocs);

%% Get Z-Scores Matrix
zmatrix = zeros(n,n);
for i=1:n % for each topography 
    curr_vect = zeros(n);
    %deep copy ICA matrix
    for h=1:n
        zmatrix(i,h) = topog(i,h);
    end
    for h=1:n
        curr_vect(h,1) = topog(i,h);
    end
    z_vect = zscore(curr_vect);
    for h=1:n
        zmatrix(i,h) = z_vect(h,1);
    end    
end
%            _                       _
% zmatrix = | z11    z12    z13   ... | row1 = IC1
%           | z21    z22    z23   ... | row2 = IC2
%           | z31    z32    z33   ... | ...
%           |_...    ...    ...   zij_|
%          col1=E1  col2=E2 ...             nxn matrix (n=#components)

%% Define scalp zones

% Find electrodes in Frontal Area (FA)
dimlefteyes=0; %number of electrodes
index1=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
   %if (chanlocs(1,k).theta > -28 && chanlocs(1,k).theta < 0) && (chanlocs(1,k).radius>0.45) && (chanlocs(1,k).radius<0.58) %electrodes are in FA
   if (chanlocs(1,k).theta > -60 && chanlocs(1,k).theta < 0) && (chanlocs(1,k).radius>0.45) && (chanlocs(1,k).radius<0.60)
        dimlefteyes=dimlefteyes+1; %count electrodes
        index1(1,dimlefteyes)=k;
    end
end

dimrighteyes=0; %number of electrodes
index2=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
   %if (chanlocs(1,k).theta < 28 && chanlocs(1,k).theta > 0) && (chanlocs(1,k).radius>0.45) && (chanlocs(1,k).radius<0.58) %electrodes are in FA
   if (chanlocs(1,k).theta < 60 && chanlocs(1,k).theta > 0) && (chanlocs(1,k).radius>0.45) && (chanlocs(1,k).radius<0.60)
        dimrighteyes=dimrighteyes+1; %count electrodes
        index2(1,dimrighteyes)=k;
    end
end

dimcentereyes=0; %number of electrodes
index3=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
    if(abs(chanlocs(1,k).theta) < 20) && (chanlocs(1,k).radius>0.45) %electrodes are in FA
        dimcentereyes=dimcentereyes+1; %count electrodes
        index3(1,dimcentereyes)=k;
    end
end
%add in code in case E5 and E10 are deleted in a 64-chan net
if dimcentereyes==0
    for k=1:nchannels
        if(abs(chanlocs(1,k).theta) == 0) && (chanlocs(1,k).radius>0.39)
            dimcentereyes=dimcentereyes+1; %count electrodes
            index3(1,dimcentereyes)=k;
        end
    end
end

% Find electrodes in everywhere else
dimcenter=0; %number of electrodes
indexc=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
    if( (abs(chanlocs(1,k).theta)) > 35 && abs(chanlocs(1,k).theta) < 109 && (chanlocs(1,k).radius<0.45) )
        dimcenter=dimcenter+1; %count electrodes
        indexc(1,dimcenter)=k;
    end
end

dimbackleft=0; %number of electrodes
indexbl=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
   if (chanlocs(1,k).theta <= -109) && (chanlocs(1,k).radius<0.55) %electrodes are in FA
        dimbackleft=dimbackleft+1; %count electrodes
        indexbl(1,dimbackleft)=k;
    end
end

dimbackright=0; %number of electrodes
indexbr=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
   if (chanlocs(1,k).theta >= 109) && (chanlocs(1,k).radius<0.55) %electrodes are in FA
        dimbackright=dimbackright+1; %count electrodes
        indexbr(1,dimbackright)=k;
    end
end

%% z-scores of the front and back

% NOTE: 
% ica struct is used for debugging purposed only
% and has no impact on the actual script run,
% just ignore it.

blinksCount = 1;
for i=1:n % for each topography
    
    % create FA electrodes vector
    zfrontlefteyes=zeros(1,dimlefteyes);
    for h=1:dimlefteyes
        zfrontlefteyes(1,h)=zmatrix(i,index1(1,h));
    end
    ica(i).lefteyeArr = zfrontlefteyes;
    
    zfrontrighteyes=zeros(1,dimrighteyes);
    for h=1:dimrighteyes
        zfrontrighteyes(1,h)=zmatrix(i,index2(1,h));
    end
    ica(i).righteyeArr = zfrontrighteyes;
    
    zcentereyes=zeros(1,dimcentereyes);
    for h=1:dimcentereyes
        zcentereyes(1,h)=zmatrix(i,index3(1,h));
    end
    ica(i).centereyeArr = zcentereyes;
    
    % create other electrodes vectors
    zcenter=zeros(1,dimcenter);
    for h=1:dimcenter
        zcenter(1,h)=zmatrix(i,indexc(1,h));
    end
    ica(i).chArr = zcenter;
    
    zbackl=zeros(1,dimbackleft);
    for h=1:dimbackleft
        zbackl(1,h)=zmatrix(i,indexbl(1,h));
    end
    ica(i).blArr = zbackl;

    zbackr=zeros(1,dimbackright);
    for h=1:dimbackright
        zbackr(1,h)=zmatrix(i,indexbr(1,h));
    end
    ica(i).brArr = zbackr;
   
    % Calculates the mean of each z-score vector
    % Sums the absolute value of each element in the vector,
    % then divides that by the number of elements.
    le = 0; 
    for j=1:dimlefteyes
        le = le + abs(zfrontlefteyes(1,j)); 
    end
    le = le/dimlefteyes;
    ica(i).lemean = le;
    
    re = 0; 
    for j=1:dimrighteyes
        re = re + abs(zfrontrighteyes(1,j)); 
    end
    re = re/dimrighteyes;
    ica(i).remean = re;
    
    ce = 0; 
    for j=1:dimcentereyes
        ce = ce + zcentereyes(1,j); %removed abs() around zcentereyes b/c they should be the same sign (1/29/19)
    end
    ce = ce/dimcentereyes;
    ica(i).cemean = ce;

    ch = 0; 
    for j=1:dimcenter
        ch = ch + abs(zcenter(1,j)); 
    end
    ch = ch/dimcenter;
    ica(i).chmean = ch;

    bl = 0; 
    for j=1:dimbackleft
        bl = bl + abs(zbackl(1,j)); 
    end
    bl = bl/dimbackleft;
    ica(i).blmean = bl;

    br = 0; 
    for j=1:dimbackright
        br = br + abs(zbackr(1,j)); 
    end
    br = br/dimbackright;
    ica(i).brmean = br;

    % if one of the means of the eye vectors pass a certain threshold and 
    % all of the back are below a certain threshold, then we reject
    if (abs(le) > 2 || abs(ce) > 2.5 || abs(re) > 2) %if high activity at front
        if (mean([abs(bl) abs(br) abs(ch)]) < 1) % if low activity in back
            if ( (var([zbackl zbackr]') < 0.15 || var([abs(zbackl) abs(zbackr)]')<0.075) && var([zfrontlefteyes zfrontrighteyes]')>var([zbackl zbackr]') ) %if low variance in back and greater variance in front
                blinks(blinksCount) = i;
                blinksCount = blinksCount + 1;
            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEPH EDIT TO CODE 1/18/19                                           %
%find all ICs with bumps btw 5 and 15 Hz                              %
%cwb=MARAcode_alphapeak(EEG);                                         %
%loop through artifacted ICs and remove any ICs that have bumps       %
remove_indices=[];                                                    %
for ii = 1:length(blinks)                                             %
    for jj = 1:length(cwb)%components_with_bumps)                     %
        if blinks(ii) == cwb(jj)%components_with_bumps(jj)            %
            remove_indices = [remove_indices ii];                     %
        end                                                           %
    end                                                               %
end                                                                   %
blinks(remove_indices)=[];                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEPH EDIT TO ADJUST CODE 2/8/19                              %
if length(blinks)>0                                            %
    % eye_electrodes is entered into the func below to check   %
    % that the eye artifacts doesn't go too far into the head  %
    eye_electrodes=[];                                         %
    for e=1:n                                                  %
        if(chanlocs(1,e).radius>0.45 && chanlocs(1,e).radius<0.54) 
            if(abs(chanlocs(1,e).theta)<17)                   %
                eye_electrodes(e) = e;                         %
            end                                                %
        end                                                    %
    end                                                        %
    eye_electrodes = nonzeros(eye_electrodes);                 %
    sec = Spatial_Info_eyes(blinks,eye_electrodes,topog,chanlocs,n); 
    remove_indices=[];                                         %
    for ii = 1:length(blinks)                                  %
        for jj = 1:length(sec)%blinks_to_keep                  %
            if blinks(ii) == sec(jj)                           %
                remove_indices = [remove_indices ii];          %
            end                                                %
        end                                                    %
    end                                                        %
    blinks(remove_indices)=[];                                 %
end                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('End Blink Detection');
