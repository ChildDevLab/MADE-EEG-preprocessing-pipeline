function ICs_to_keep = Spatial_Info_eyes(blinks,eye_electrodes,topog,chanlocs,n)
% function discs_to_keep = Spatial_Info_discs(EEG)
%
% figures out if eye stuff rejected by adjust should be kept
%    - keeps them if there's too much spread

%define variables
ICs_to_keep = [];
arb_dist = 5.5;
arb_zval = 2;

%Code from Beall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Z-Scores Matrix
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%list of channels on the outside of the head
% Find electrodes in outer ring
dimouterring=0; %number of electrodes
indexor=zeros(1,n); %indices of electrodes
for k=1:n
   if (chanlocs(1,k).radius>0.51)
        dimouterring=dimouterring+1; %count electrodes
        indexor(1,dimouterring)=k;
   end
end
outer_ring = nonzeros(indexor);
% 64-chan net: outer_ring = [1 17 18 23 24 29 32 43 47 52 55 58];
% 128-chan net: outer_ring = [48 43 49 56 63 68 73 81 88 94 99 107 113 120 119 17 128 125 127 126];
outer_ring = unique([outer_ring' eye_electrodes']);

dimhigherthresh=0;
indexht=zeros(1,n);
for k=1:n
    if (chanlocs(1,k).radius>0.35 && chanlocs(1,k).radius<0.5 && abs(chanlocs(1,k).theta)<20)
        dimhigherthresh=dimhigherthresh+1;
        indexht(1,dimhigherthresh)=k;
    end
end
higherthresh = nonzeros(indexht)';

keep_IC = []; % will be the list of ICs to keep
for ic=1:n
    if sum(ic==blinks)>0 % if this IC matches any numbers in the eye artifact list
        % find any channels/electrodes that have high activity (z>2)
        chans=find(abs(zmatrix(ic,:))>arb_zval); %grabs any chans/electrodes that are above thresh
        
        % delete the outer ring(s) of electrodes
        % We only care about the spread over more central sites for eye stuff
        chans_to_delete=[];
        for e=1:length(chans)
            if find(outer_ring == chans(e)) > 0
                % delete this chan/electrode from the list of chans (we don't care about it)
                chans_to_delete=[chans_to_delete e];
            end
            if find(higherthresh == chans(e)) > 0 
                if zmatrix(ic,chans(e))<2.5
                    % this channel isn't above the higher threshold set for this ring on the net
                    chans_to_delete=[chans_to_delete e];
                end
            end
        end
        chans(unique(chans_to_delete))=[];
        
        % calculate the distance between the points
        pts_matrix_3D=[];
        dist_matrix=[];
        for c=1:length(chans) % loop through list of chans above the threshold value and get a list of points
            pts_matrix_3D=[pts_matrix_3D; chanlocs(1,chans(c)).X chanlocs(1,chans(c)).Y chanlocs(1,chans(c)).Z];
        end
        dist_matrix = pdist(pts_matrix_3D,'euclidean'); % calculate distances btw all pts. above the threshold
        
        % figure out if component should be kept
        corresponding_pts = []; % gives us the indices of the points
        for jj=1:length(chans)-1
            corresponding_pts = [corresponding_pts; (jj*ones(1,length(chans)-jj))' ((jj+1):length(chans))'];
        end
        dists_below_thresh = find(dist_matrix<arb_dist); % tells us distances below the threshold
        pts_within_dist = corresponding_pts(dists_below_thresh,:); % tells us which point pairs are within the distnace threshold
        within_diff = [];
        if length(pts_within_dist)>0 % if any pts are within the arbitrary distance
            most_com = mode(pts_within_dist(:)); % lets us know if eye stuff spreads too far into the head
            if length(find(pts_within_dist(:) == most_com))>1
                keep_IC = [keep_IC ic];
            end
        end
        ICs_to_keep = unique(keep_IC);
    else
        % skip (we don't need to waste time looking at it if it won't be rejected)
    end
end

end % end of function