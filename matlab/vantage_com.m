function [com2bbMean, com3d, com2bb] = vantage_com(varargin)
%VANTAGE_COM Calculates a cyclist's center of mass position using
%reflective marker locations from a Retul Vantage motion capture system.
%
% This script uses an output file (.pose) from the Retul Vantage Motion
% Capture system (left or right side) and generates a prediction and
% animation of the rider's Center of Mass (CoM) position and displacement.
% This animation is saved as a video file (.mp4). The mean position of the
% rider's CoM relative to the bottom bracket over the entire capture period
% is also calculated.
%
%   Examples:
%   *Note: inputs need to be left blank or in exact order. I.e. next input
%   in sequence cannot be accessed if previous input is not present.
%   Default values can be entered if necessary.
%   Order sequence: file, anim, sex, mass, sacralAngle, sh_width, hp_width
%       - com3d = vantage_com 
%       - [com3d, com2bb, com2bbMean] = vantage_com 
%       - vantage_com('foo.pose')
%       - vantage_com('foo.pose','off')
%       - vantage_com('foo.pose','off','female')
%       - vantage_com('foo.pose','off','female', 60)
%       - vantage_com('foo.pose','off','female', 60, 63)
%       - vantage_com('foo.pose','off','female', 60, 63, 343)
%       - vantage_com('foo.pose','off','female', 60, 63, 343, 288)
%
%   Inputs:
%        file - string containing the pose file name that you wish to load
%        (leave blank to choose from dialog box). File needs to exist in
%        current working directory.
%
%        anim - string containing whether to produce an animation of the
%        rider 'on' or 'off' (default - 'on')
%        
%        sex - string containing the sex of rider 'male' or 'female'
%        (default - 'male')
%        
%        mass - mass of rider in kilograms (default - 78.05 kg for
%        male, 64.26 kg for female)
%
%        sacralAngle - angle of lower back counter-clockwise from
%        horizontal (default - 54 degrees.)
%
%        sh_width - width of rider's shoulders in millimeters (default -
%        411 for male, 367 for female)
%
%        hp_width - width of rider's hips in millimeters (default - 296 for
%        male, 291 for female)  
%
%        shoeMass - mass of each shoe in kilograms (default - 1 kg)
% 
%     Outputs: 
%         com3d - mx3 array of frame-by-frame three-dimensional CoM
%         position within the vantage-coordinate system, m = number of
%         frames. Data is upsampled to 200 Hz.
%
%         com2bb - mx1 vector of fore-aft CoM position relative to the
%         estimated bottom bracket position. Data is upsampled to 200 Hz.
%
%         com2bbMean - mean value of com2bb across all frames
%
% Copyright (C) Ross Wilkinson 2020- vantage_com.m by Ross Duncan Wilkinson
% is licensed under a Creative Commons
% Attribution-NonCommercial-NoDerivatives 4.0 International License.
% ======================================================================= %
% Background: Determining the human body's center of mass is an important
% tool for analysing the biomechanics and energetics of human motion. In a
% biomechanical research setting, the most accurate estimates of each
% segment's CoM position requires placing and tracking the
% three-dimensional position of more than 38 markers (Tisserand et al.,
% 2016). This method is expensive and time consuming, which is impractical
% for certain applications like bike fitting. Therefore another approach is
% to use a reduced number of markers to estimate whole body CoM position
% (Dumas & Wojtusch, 2017). In either case, the technique involves
% determining the end points of each segment and estimates of body segment
% inertial parameters (BSIPs). BSIPs include BSIPs can be obtained in
% different ways including direct measurements on cadavers or
% photogrammetry and medical imaging on living humans, but they are more
% generally estimated by regression equations (based on those
% measurements).
% 
% The following approach uses BSIPs based on the regression equations of De
% Leva (1996) adjusted from the data of Zatsiorsky et al. (1990) in
% combination with Retul Vantage data (8 markers) to estimate the
% whole-body CoM position of a 16-segment rigid body biomechanical model
% (head with neck, upper trunk, middle trunk, lower trunk, upper arm (x2),
% forearm (x2), hand (x2), thigh (x2), shank (x2), and foot (x2).
% 
% Beyond the limitations inherent to estimating BSIPs, the main assumptions
% for this approach are:
% 
% # Retul Vantage marker placements correspond to segment end-points
% # Motion of the right and left limbs are symmetrical
% # The length of each subject's "head with neck" segment is the same
% within each sex
% # The "head with neck" segment is roughly aligned with the upper trunk
% (~45 degrees)
% # The length of each hand is 0
% # The length of each foot is from the calcaneus to the MTP joint
% 
% *References*
% 
% # de Leva, P. (1996). Adjustments to Zatsiorsky-Seluyanov?s segment
% inertia parameters. _Journal of Biomechanics_, _29_(9), 1223?1230.
% <https://doi.org/10.1016/0021-9290(95)00178-6
% https://doi.org/10.1016/0021-9290(95)00178-6>
% # Tisserand, R., Robert, T., Dumas, R., & Chèze, L. (2016). A simplified
% marker set to define the center of mass for stability analysis in dynamic
% situations. _Gait and Posture_, _48_, 64?67.
% <https://doi.org/10.1016/j.gaitpost.2016.04.032
% https://doi.org/10.1016/j.gaitpost.2016.04.032>
% # Dumas, R., & Wojtusch, J. (2017). Estimation of the Body Segment
% Inertial Parameters for the Rigid Body Biomechanical Models Used in
% Motion Analysis. In _Handbook of Human Motion_.
% <https://doi.org/10.1007/978-3-319-30808-1
% https://doi.org/10.1007/978-3-319-30808-1>
% ======================================================================= %

%% Set default values
animate = 'on'; sex = 'male'; mass = 78.05; sacralAngle = 54; sh_width = 411;
hp_width = 296; shoeMass = 1;
if contains(computer('arch') ,'win')
    path = [pwd '\'];
else
    path = [pwd '/'];
end

%% Evaluate inputs
switch nargin
    case 0
        [file, path] = uigetfile('*.*');
    case 1
        file = varargin{1};
    case 2
        file = varargin{1}; animate = varargin{2};
    case 3
        file = varargin{1}; animate = varargin{2}; sex = varargin{3};
        if strcmp(sex,'female')
            mass = 64.26; sh_width = 367; hp_width = 291;
        end  
    case 4
        file = varargin{1}; animate = varargin{2}; sex = varargin{3};
        mass = varargin{4};
        if strcmp(sex,'female')
            sh_width = 367; hp_width = 291;
        end 
    case 5
        file = varargin{1}; animate = varargin{2}; sex = varargin{3};
        mass = varargin{4}; sacralAngle = varargin{5};
        if strcmp(sex,'female')
            sh_width = 367; hp_width = 291;
        end
    case 6
        file = varargin{1}; animate = varargin{2}; sex = varargin{3};
        mass = varargin{4}; sacralAngle = varargin{5}; sh_width = varargin{6};
        if strcmp(sex,'female')
             hp_width = 291;
        end
    case 7
        file = varargin{1}; animate = varargin{2}; sex = varargin{3};
        mass = varargin{4}; sacralAngle = varargin{5}; sh_width = varargin{6};    
        hp_width = varargin{7};
    otherwise
        file = varargin{1}; animate = varargin{2}; sex = varargin{3};
        mass = varargin{4}; sacralAngle = varargin{5}; sh_width = varargin{6};    
        hp_width = varargin{7}; shoeMass = varargin{8};       
end

%% Save file in .xml format
% The file is then saved in .xml format within the same directory as the
% original .pose file.
cd(path)
filename = strrep(file,'.pose','.xml');
copyfile(file,filename)

%% Load Vantage data into structure
% The Vantage Data is loaded into a MATLAB structure. 
tree = xml_read([path filename]);

%% Place marker data into new structure 'S'
% First, the 'COMMENT' field is removed to leave only the marker names as
% fields. The time and position coordinates for each marker at each frame
% are converted from a field to an element within each array. This just
% makes it alot easier to perform calculations on.
tree.stxyz = rmfield(tree.stxyz, 'COMMENT');
markerList = fieldnames(tree.stxyz(1));
nMarkers = length(fieldnames(tree.stxyz(1))); 
nFrames = length(tree.stxyz);

for i = 1:nMarkers
    for j = 1:nFrames
        str = tree.stxyz(j).(markerList{i});
        S.(markerList{i})(j,:) = str2double(split(str));
    end
end

%% Set Body Segment Inertial Parameters (BSIPs)
% Briefly, the segment mass is computed as a percentage of the body mass.
% The position of the segment's CoM is computed as a percentage of the
% segment length, defined as the distance between the segment end-points.
nSegments = 16;

if strcmp(sex, 'male') % Male
    % Mass (% of whole body)
    M = [0.0694;... % head
        0.1596;...  % upper trunk
        0.1633;...  % middle trunk
        0.1117;...  % lower trunk
        0.0271;...  % upper arm right
        0.0162;...  % forearm right
        0.0061;...  % hand right
        0.1416;...  % thigh right
        0.0433;...  % shank right
        0.0137;...  % foot right
        0.0271;...  % upper arm left
        0.0162;...  % forearm left
        0.0061;...  % hand left
        0.1416;...  % thigh left
        0.0433;...  % shank left
        0.0137];    % foot left
    % Length to segment CoM (% from proximal end-point)
    L = [0.5002;... % head
        0.2999;...  % upper trunk
        0.4502;...  % middle trunk
        0.6115;...  % lower trunk
        0.5772;...  % upper arm right
        0.4608;...  % forearm right
        0.7948;...  % hand right
        0.4095;...  % thigh right
        0.4459;...  % shank right
        0.4415;...  % foot right
        0.5772;...  % upper arm left
        0.4608;...  % forearm left
        0.7948;...  % hand left
        0.4095;...  % thigh left
        0.4459;...  % shank left
        0.4415];    % foot left     
else % Female
    % Mass (% of whole body)
    M = [0.0668;... % head
        0.1545;...  % upper trunk
        0.1465;...  % middle trunk
        0.1247;...  % lower trunk
        0.0255;...  % upper arm right
        0.0138;...  % forearm right
        0.0056;...  % hand right
        0.1478;...  % thigh right
        0.0481;...  % shank right
        0.0129;...  % foot right
        0.0255;...  % upper arm left
        0.0138;...  % forearm left
        0.0056;...  % hand left
        0.1478;...  % thigh left
        0.0481;...  % shank left
        0.0129];    % foot left
    % Length to segment CoM (% from proximal to end-point)
    L = [0.4841;... % head
        0.2077;...  % upper trunk
        0.4512;...  % middle trunk
        0.4920;...  % lower trunk
        0.5754;...  % upper arm right
        0.4592;...  % forearm right
        0.7534;...  % hand right
        0.3612;...  % thigh right
        0.4416;...  % shank right
        0.4014;...  % foot right
        0.5754;...  % upper arm left
        0.4559;...  % forearm left
        0.7474;...  % hand left
        0.3612;...  % thigh left
        0.4416;...  % shank left
        0.4014];    % foot left
end

%% add shoe mass to each foot
relativeShoeMass = shoeMass / mass;
M = M - (relativeShoeMass/nSegments);
M(10) = M(10) + (relativeShoeMass/2);
M(16) = M(10);

%% Upsample time series data
% Vantage markers are initialized asynchronously. Upsampling the data from
% 18 Hz to 200 Hz decreases the absolute time delay in initialization
% between markers, which reduces the synchronization error.
for i = 1:nMarkers
    start_time = S.(markerList{i})(1,2);
    end_time = S.(markerList{i})(end,2);
    freq = 1 / ((end_time - start_time) / nFrames);
    new_freq = 200;
    A.(markerList{i})(:,1) = interp(S.(markerList{i})(:,2),round(new_freq/freq));
    A.(markerList{i})(:,2) = interp(S.(markerList{i})(:,3),round(new_freq/freq));
    A.(markerList{i})(:,3) = interp(S.(markerList{i})(:,4),round(new_freq/freq));
    A.(markerList{i})(:,4) = interp(S.(markerList{i})(:,5),round(new_freq/freq));
end

%% Synchronize time series data
% The wrist marker is initialized last. Therefore, clip all other marker
% signals collected prior to the first frame of the wrist marker.
for i = 1:nMarkers
    t1 = find(A.(markerList{i})(:,1) >= A.wr(1,1),1,'first');
    trim = ceil(t1/30)*30;
    C.(markerList{i})(:,1) = A.(markerList{i})(t1:t1+(end-trim),1);
    C.(markerList{i})(:,2) = A.(markerList{i})(t1:t1+(end-trim),2);
    C.(markerList{i})(:,3) = A.(markerList{i})(t1:t1+(end-trim),3);
    C.(markerList{i})(:,4) = A.(markerList{i})(t1:t1+(end-trim),4);
end

%% Create symmetrical contralateral markers
% Phase-shift ipsilateral markers by finding the half cycle period of the
% x-coordinate of the meta-tarsal marker ("ft") signal. This should be a
% good estimate of the time taken to complete a half crank cycle.
[~, locs] = findpeaks(C.ft(:,2),'minPeakDistance',60);
wavelength = mean(diff(locs));
delay = round(wavelength/2);

% At the same time, adjust marker names.
for i = 1:nMarkers
    signal = C.(markerList{i});
    if strcmp(tree.viewedSide,'R')
        C.([markerList{i} '_r']) = signal(1:end-(delay-1),:);
        C.(markerList{i}) = signal(delay:end,:);
    else
        C.([markerList{i} '_r']) = signal(delay:end,:);
        C.(markerList{i}) = signal(1:end-(delay-1),:);
    end
end

%% Create estimate of lower trunk length
% Use BSIPs from De Leva (1996) for the distribution of trunk length.
% Use estimated sacral angle of 54 deg. or use measured angle from vantage
% data.
for i =  1:length(C.sh)
    trunkLength(i) = norm(mean([C.sh(i,[2 3]); C.sh_r(i,[2 3])])...
        - mean([C.hp(i,[2 3]) C.hp_r(i,[2 3])]));
end

if strcmp(sex, 'male') % Male
    % use male full trunk length average from DeLeva to get initial estimate of lower trunk length
    lowertrunkLength = trunkLength * (145.7/531.9); 
else % Female
    % use female full trunk length average from DeLeva to get initial estimate of lower trunk length
    lowertrunkLength = trunkLength * (181.5/529.3);
end

%% Create virtual marker at proximal endpoint of lower trunk
lt_x = cos(deg2rad(sacralAngle)) * lowertrunkLength;
lt_y = sin(deg2rad(sacralAngle)) * lowertrunkLength;

if strcmp(tree.viewedSide,'R')
    C.lt = C.hp_r;
    C.lt(:,2) = C.hp_r(:,2) - lt_x';
    C.lt(:,3) = C.hp_r(:,3) - lt_y';
else
    C.lt(:,2) = C.hp(:,2) + lt_x'; % +ve for left
    C.lt(:,3) = C.hp(:,3) - lt_y';
end

%% Create estimate of head and middle trunk length
% Use BSIPs from De Leva (1996) for head length. Estimate whole trunk
% length using shoulder and hip markers. Use lower trunk marker to
% determine residual trunk length to shoulder.
for i = 1:length(C.sh)
    residualtrunkLength(i) = ...
        norm(C.lt(i,[2 3]) - mean([C.sh(i,[2 3]); C.sh_r(i,[2 3])]));
    residualtrunkAngle(i) = ...
        angle2Points(mean([C.sh(i,[2 3]); C.sh_r(i,[2 3])]),C.lt(i,[2 3]));
end

if strcmp(sex, 'male') % Male
    headLength = 242.9;
    middletrunkLength = residualtrunkLength * (215.5/(170.7+215.5));
else % Female
    headLength = 243.7;
    middletrunkLength = residualtrunkLength * (205.3/(142.5+205.3));
end

%% Create virtual markers at the proximal enpoint of the head and middle trunk
% Estimate the position of vertex of head by adding head length to shoulder
% position and using residual trunk angle. Estimate the position of the
% proximal end of the middle trunk by adding middle trunk length to lower
% trunk marker at residual trunk angle.
mt_x = cos(residualtrunkAngle) .* middletrunkLength;
mt_y = sin(residualtrunkAngle) .* middletrunkLength;
hd_x = cos(residualtrunkAngle) .* headLength; 
hd_y = sin(residualtrunkAngle) .* headLength;
    
if strcmp(tree.viewedSide,'R')
    C.hd = C.sh_r;
    C.hd(:,2) = C.sh_r(:,2) - hd_x';
    C.hd(:,3) = C.sh_r(:,3) - hd_y';
    C.mt = C.lt;
    C.mt(:,2) = C.lt(:,2) - mt_x';
    C.mt(:,3) = C.lt(:,3) - mt_y';
else
    C.hd(:,2) = C.sh(:,2) - hd_x'; %-ve for left as well
    C.hd(:,3) = C.sh(:,3) - hd_y'; %-ve for left as well
    C.mt = C.lt; % add this for left side files as well
    C.mt(:,2) = C.lt(:,2) - mt_x'; % -ve for left as well
    C.mt(:,3) = C.lt(:,3) - mt_y'; % C.lt not C.hp
end

%% Adjust z-coordinates of contralateral markers
% Use CDC data (<https://www.cdc.gov/nchs/data/series/sr_11/sr11_249.pdf
% https://www.cdc.gov/nchs/data/series/sr_11/sr11_249.pdf>) for average
% shoulder and hip breadth for males and females. Then adjust z-coordinates
% of contralateral head, shoulder, elbow, wrist, hip, knee, ankle, heel,
% and foot markers.
if strcmp(tree.viewedSide,'R')
    C.hd(:,4) = C.sh_r(:,4) - (sh_width/2);
    C.mt(:,4) = C.sh_r(:,4) - (sh_width/2);
    C.lt(:,4) = C.sh_r(:,4) - (sh_width/2);
    C.sh(:,4) = C.sh_r(:,4) - sh_width;
    C.el(:,4) = C.sh_r(:,4) - sh_width - (C.el_r(:,4)-C.sh_r(:,4));
    C.wr(:,4) = C.sh_r(:,4) - sh_width - (C.wr_r(:,4)-C.sh_r(:,4));
    C.hp(:,4) = C.hp_r(:,4) - hp_width;
    C.kn(:,4) = C.hp_r(:,4) - hp_width - (C.kn_r(:,4)-C.hp_r(:,4));
    C.an(:,4) = C.hp_r(:,4) - hp_width - (C.an_r(:,4)-C.hp_r(:,4));
    C.he(:,4) = C.hp_r(:,4) - hp_width - (C.he_r(:,4)-C.hp_r(:,4));
    C.ft(:,4) = C.hp_r(:,4) - hp_width - (C.ft_r(:,4)-C.hp_r(:,4));
else
    C.hd(:,4) = C.sh(:,4) - (sh_width/2);
    C.mt(:,4) = C.sh(:,4) - (sh_width/2);
    C.lt(:,4) = C.sh(:,4) - (sh_width/2);
    C.sh_r(:,4) = C.sh(:,4) - sh_width;
    C.el_r(:,4) = C.sh(:,4) - sh_width - (C.el(:,4)-C.sh(:,4));
    C.wr_r(:,4) = C.sh(:,4) - sh_width - (C.wr(:,4)-C.sh(:,4));
    C.hp_r(:,4) = C.hp(:,4) - hp_width;
    C.kn_r(:,4) = C.hp(:,4) - hp_width - (C.kn(:,4)-C.hp(:,4));
    C.an_r(:,4) = C.hp(:,4) - hp_width - (C.an(:,4)-C.hp(:,4));
    C.he_r(:,4) = C.hp(:,4) - hp_width - (C.he(:,4)-C.hp(:,4));
    C.ft_r(:,4) = C.hp(:,4) - hp_width - (C.ft(:,4)-C.hp(:,4));
end

%% Convert data to left side Vantage coordinate system
% Left side view XYZ = right side view (-X,Y,-Z).
for i = 1:nMarkers
    C.(markerList{i})(:,[2 4]) = C.(markerList{i})(:,[2 4]) * tree.forward(1);
    C.([markerList{i} '_r'])(:,[2 4]) = C.([markerList{i} '_r'])(:,[2 4]) * tree.forward(1);
end
C.hd(:,[2 4]) = C.hd(:,[2 4]) * tree.forward(1);
C.mt(:,[2 4]) = C.mt(:,[2 4]) * tree.forward(1);
C.lt(:,[2 4]) = C.lt(:,[2 4]) * tree.forward(1);

%% Estimate segment end points
% Use marker positions to estimate proximal and distal end-points of each segment.
nSamples = length(C.sh);
[prox, dist] = deal(zeros(nSegments,3,nSamples));

for i = 1:nSamples
    % Proximal Marker                           % Segment           % Proximal Endpoint    
    % ===============                           % =======           % =================
    prox(:,:,i) = [...
        C.hd(i,2:4);...                         % head              vertex
        mean([C.sh(1,2:4);C.sh_r(1,2:4)]);...   % upper trunk       jugular notch
        C.mt(i,2:4);...                         % middle trunk      xyphion
        C.lt(i,2:4);...                         % lower trunk       omphalion
        C.sh_r(i,2:4);...                       % upper arm right   shoulder JC right
        C.el_r(i,2:4);...                       % forearm right     elbow JC right
        C.wr_r(i,2:4);...                       % hand right        stylion right
        C.hp_r(i,2:4);...                       % thigh right       hip JC right      
        C.kn_r(i,2:4);...                       % shank right       knee JC right
        C.he_r(i,2:4);...                       % foot right        calcaneus right
        C.sh(i,2:4);...                         % upper arm left    shoulder JC left
        C.el(i,2:4);...                         % forearm left      elbow JC left
        C.wr(i,2:4);...                         % hand left         stylion left         
        C.hp(i,2:4);...                         % thigh left        hip JC left
        C.kn(i,2:4);...                         % shank left        knee JC left
        C.he(i,2:4)];                           % foot left         calcaneus left
    
    % Distal Marker                             % Segment           % Distal Endpoint    
    % =============                             % =======           % ===============
    dist(:,:,i) = [...
        mean([C.sh(1,2:4);C.sh_r(1,2:4)]);...   % head              mid cervicale
        C.mt(i,2:4);...                         % upper trunk       xyphion
        C.lt(i,2:4);...                         % middle trunk      omphalion
        mean([C.hp_r(i,2:4);C.hp(i,2:4)]);...   % lower trunk       mid hip JC
        C.el_r(i,2:4);...                       % upper arm right   elbow JC
        C.wr_r(i,2:4);...                       % forearm right     stylion
        C.wr_r(i,2:4);...                       % hand right        3rd metacarpale
        C.kn_r(i,2:4);...                       % thigh right       knee JC
        C.an_r(i,2:4);...                       % shank right       lateral malleolus    
        C.ft_r(i,2:4);...                       % foot right        toe tip
        C.el(i,2:4);...                         % upper arm left    elbow JC
        C.wr(i,2:4);...                         % forearm left      stylion
        C.wr(i,2:4);...                         % hand left         3rd metacarpale
        C.kn(i,2:4);...                         % thigh left        knee JC
        C.an(i,2:4);...                         % shank left        lateral maleolus
        C.ft(i,2:4)];                           % foot left         toe tip
    % create segments for animation 
    head(:,:,i) = [C.hd(i,2:4); mean([C.sh_r(i,2:4);C.sh(i,2:4)])];
    arm_r(:,:,i) = [C.sh_r(i,2:4); C.el_r(i,2:4); C.wr_r(i,2:4)];
    arm_l(:,:,i) = [C.sh(i,2:4); C.el(i,2:4); C.wr(i,2:4)];
    trunk(:,:,i) = [mean([C.sh(1,2:4); C.sh_r(1,2:4)]); C.mt(i,2:4);...
        C.lt(i,2:4); mean([C.hp_r(i,2:4);C.hp(i,2:4)])];
    leg_r(:,:,i) = [C.hp_r(i,2:4); C.kn_r(i,2:4); C.an_r(i,2:4); C.he_r(i,2:4);...
        C.ft_r(i,2:4); C.an_r(i,2:4)];    
    leg_l(:,:,i) = [C.hp(i,2:4); C.kn(i,2:4); C.an(i,2:4); C.he(i,2:4);...
        C.ft(i,2:4); C.an(i,2:4)];    
end

%% Determine segment CoM coordinates and segment torque about origin
% The center of mass is an ideal point about which the torques due to body
% segment weights is zero. Segment CoM coorindates will be equal to the
% proximal coordinates plus the relative length to the CoM towards the
% distal coordinates.
segCOM = prox + L .* (dist-prox);

%% Determine segment torque about origin
% Segment torque around the origin will be equal to the product of the CoM
% coordinates and the relative mass of the segment.
segTRQ = segCOM .* M;

%% Determine whole body CoM coordinates
% Sum the torques about the origin.
wbCOM = sum(segTRQ);

%% Determine CoM relative to BB (anterior-posterior)
% Estimate bottom bracket position using left and right meta-tarsal
% markers. Then estimate the anterior-posterior position of the CoM
% relative to the bottom bracket.
bb = mean([dist(10,:,:);dist(16,:,:)]);
bb = reshape(bb,[3,length(bb)])';
com2bb = bb(:,1) - reshape(wbCOM(:,1,:),[length(wbCOM),1]);
com2bbMean = mean(com2bb);
com3d = wbCOM;

%% Save data
save(strrep(filename,'.xml','.mat'))

%% Calculate right crank angle
% Use the meta-tarsal marker of the right foot and the bottom bracket position 
% to estimate crank angle.
P1 = C.ft_r(:,2:3);
P2 = bb(:,1:2);
P3 = [bb(:,1) bb(:,2)+1];
angle = (2*pi - angle3Points(P1,P2,P3)) * (180/pi);

figure('Name','Crank Angle')
plot(angle)
xlabel('Samples (200 Hz)')
ylabel('Crank angle (degrees)')

[pks, locs] = findpeaks(angle);
n = 0;
for i = 1:length(locs)-1
    n = n + pks(i);
    angle(locs(i)+1:locs(i+1)) = angle(locs(i)+1:locs(i+1))+n;
    if i == length(locs)-1
        n = n + pks(end);
        angle(locs(end)+1:end) = angle(locs(end)+1:end)+n;
    end
end

%% Plot initial rider and CoM position
% Coordinate system changes depending on which side the Vantage data is
% filmed from.
% 
% For example:
% * Vantage system (Left) (X,Y,Z) = Vantage system (Right) (-X,Y,-Z)
% * Vantage system (Left) (X,Y,Z) = MATLAB coordinate system (-Y,-Z,X)
% * Vantage system (Right) (X,Y,Z) = MATLAB coordinate system (Y,-Z,-X)

%% Subplot 1: 3D rider
fig = figure;
color1 = 'k';
color2 = 'r';
set(fig,'name','Rider Animation','color','w','position',[30 80 1220 620])
ax1 = subplot(3,3,[1 4 7]);
p1 = plot3(head(:,3,1),-head(:,1,1),-head(:,2,1),'-o',...
    'color',color1,'markerFaceColor',color1,'markerSize',4); % head
hold on
p2 = plot3(arm_r(:,3,1),-arm_r(:,1,1),-arm_r(:,2,1),'r-o',...
    'color',color2,'markerFaceColor',color2,'markerSize',4); % right arm red
p3 = plot3(arm_l(:,3,1),-arm_l(:,1,1),-arm_l(:,2,1),'-o',...
    'color',color1,'markerFaceColor',color1,'markerSize',4); % left arm
p4 = plot3(trunk(:,3,1),-trunk(:,1,1),-trunk(:,2,1),'-o',...
    'color',color1,'markerFaceColor',color1,'markerSize',4); % torso
p5 = plot3(leg_r(:,3,1),-leg_r(:,1,1),-leg_r(:,2,1),'-o',...
    'color',color2,'markerFaceColor',color2,'markerSize',4); % right leg red
p6 = plot3(leg_l(:,3,1),-leg_l(:,1,1),-leg_l(:,2,1),'-o',...
    'color',color1,'markerFaceColor',color1,'markerSize',4); % left leg
p7 = scatter3(wbCOM(:,3,1),-wbCOM(:,1,1),-wbCOM(:,2,1),50,'g','filled'); % CoM
p8 = line([dist(10,3,1) dist(10,3,1)],[-dist(10,1,1) -bb(1,1)],...
    [-dist(10,2,1) -bb(1,2)],'Color','k','LineStyle','-','LineWidth',10); % cr_r
p9 = line([dist(16,3,1) dist(16,3,1)],[-dist(16,1,1) -bb(1,1)],...
    [-dist(16,2,1) -bb(1,2)],'Color','k','LineStyle','-','LineWidth',10); % cr_l
p10 = line([dist(10,3,1) dist(16,3,1)],[-bb(1,1) -bb(1,1)],[-bb(1,2) -bb(1,2)],...
    'Color','k','LineStyle','-','LineWidth',5); % spindle
% p11 = line([bb(3) bb(3)],[-bb(1) -bb(1)],[-600 1000],'Color','k',...
% 'LineStyle','--', 'LineWidth',1); % BB line
p12 = line([C.sh(1,4) C.sh_r(1,4)],[-C.sh(1,2) -C.sh_r(1,2)],[-C.sh(1,3)...
    -C.sh_r(1,3)],'Color',[1 1 1]*0.5,'LineStyle','-','LineWidth',1); % cervicale
p13 = line([C.hp(1,4) C.hp_r(1,4)],[-C.hp(1,2) -C.hp_r(1,2)],[-C.hp(1,3)...
    -C.hp_r(1,3)],'Color',[1 1 1]*0.5,'LineStyle','-','LineWidth',1); % pelvis

% Set axis limits and labels
ax1.View = [-90 0];
midline = mean(C.hd(:,4));
axis equal
axis([midline-500 midline+500 -bb(1,1)-800 -bb(1,1)+600 -bb(1,2)-400 ...
    -bb(1,2)+1600])
grid on
box on
title({'Rider position (3D View)','\rmred = right side'})
xlabel('z (mm)')
ylabel('x (mm)')
zlabel('y (mm)')

%% Subplot 2: 3D CoM
subplot(3,3,[2 5 8])
xDis = wbCOM(:,3,:) - mean(wbCOM(:,3,:));
yDis = -wbCOM(:,1,:) - mean(-wbCOM(:,1,:));
zDis = -wbCOM(:,2,:) - mean(-wbCOM(:,2,:));

lim = ceil(max(abs([xDis(:); yDis(:); zDis(:)])*1.25));
if lim < 30
    lim = 30;
end

p14 = scatter3(xDis(1),yDis(1),zDis(1),300,'g','filled'); % CoM disp
hold on
p15 = line([0 0],[-lim lim],[0 0],'Color','k','LineStyle','-','LineWidth',1);
p16 = line([-lim lim],[0 0],[0 0],'Color','k','LineStyle','-','LineWidth',1);
p17 = line([0 0],[0 0],[-lim lim],'Color','k','LineStyle','-','LineWidth',1);
box on
grid on
axis equal
axis([-lim lim -lim lim -lim lim])
title({'CoM Displacement (3D View)','\rm[0,0,0] = mean CoM position'})
xlabel('\Deltaz (mm)')
ylabel('\Deltax (mm)')
zlabel('\Deltay (mm)')

%% Plot CoM projection on each plane [x-y, y-z, x-z]
p21 = animatedline;
p22 = animatedline;
p23 = animatedline;

%% Subplot 3: CoM displacement (x-axis) with respect to crank angle
% Movement in the x-axis represents anterior-posterior displacement.
subplot(3,3,3)
p24 = animatedline('LineStyle','--');
hold on
p18 = animatedline;
axis([0 angle(end) -lim lim])
ax = gca;
ax.XTick = 0:180:360*(length(locs)+1);
ax.XTickLabels = repmat({'','TDC'},1,length(locs)+1);
str = num2str(round(com2bbMean));
title({'CoM Displacement (fwd-bck)',['\rmBB to CoM = ' str 'mm']})
ylabel({'\Deltax (mm)','fwd \leftrightarrow bck'})
grid on
box on

%% Subplot 4: CoM displacement (y-axis) with respect to crank angle
% Movement in the y-axis represents inferrior-superior displacement.
subplot(3,3,6)
p19 = animatedline;
axis([0 angle(end) -lim lim])
ax = gca;
ax.XTick = 0:180:360*length(locs)+1;
ax.XTickLabels = repmat({'','TDC'},1,length(locs)+1);
title('CoM Displacement (vertical)')
ylabel({'\Deltay (mm)','down \leftrightarrow up'})
grid on
box on

%% Subplot 5: CoM displacement (z-axis) with respect to crank angle
% Movement in the z-axis represents medial-lateral displacement.
subplot(3,3,9)
p20 = animatedline;
axis([0 angle(end) -lim lim])
ax = gca;
ax.XTick = 0:180:360*length(locs)+1;
ax.XTickLabels = repmat({'','TDC'},1,length(locs)+1);
title('CoM Displacement (right-left)')
xlabel('Crank angle (Right TDC-TDC)')
ylabel({'\Deltaz (mm)','right \leftrightarrow left'})
grid on
box on

%% Choose whether to visualise data
if strcmp(animate, 'on')   
    % Select directory to save video
    filename = strrep(filename,'.xml','.mp4');
    v = VideoWriter([path '/' filename],'MPEG-4');

    %% Animate all subplots
    % Loop through data from each sample point. Do not close the figure!
    open(v)
    for j = 1:3:length(prox)  
        set(p1,'xdata',head(:,3,j),'ydata',-head(:,1,j),'zdata',-head(:,2,j));
        set(p2,'xdata',arm_r(:,3,j),'ydata',-arm_r(:,1,j),'zdata',-arm_r(:,2,j));
        set(p3,'xdata',arm_l(:,3,j),'ydata',-arm_l(:,1,j),'zdata',-arm_l(:,2,j));
        set(p4,'xdata',trunk(:,3,j),'ydata',-trunk(:,1,j),'zdata',-trunk(:,2,j));
        set(p5,'xdata',leg_r(:,3,j),'ydata',-leg_r(:,1,j),'zdata',-leg_r(:,2,j));
        set(p6,'xdata',leg_l(:,3,j),'ydata',-leg_l(:,1,j),'zdata',-leg_l(:,2,j));
        set(p7,'xdata',wbCOM(:,3,j),'ydata',-wbCOM(:,1,j),'zdata',-wbCOM(:,2,j));
        set(p8,'xdata',[dist(10,3,j) dist(10,3,j)],'ydata',[-dist(10,1,j) -bb(j,1)],...
            'zdata',[-dist(10,2,j) -bb(j,2)]); % crank right
        set(p9,'xdata',[dist(16,3,j) dist(16,3,j)],'ydata',[-dist(16,1,j) -bb(j,1)],...
            'zdata',[-dist(16,2,j) -bb(j,2)]); % crank left
        set(p12,'xdata',[C.sh(j,4) C.sh_r(j,4)],'ydata',[-C.sh(j,2) -C.sh_r(j,2)],...
            'zdata',[-C.sh(j,3) -C.sh_r(j,3)]); % cervicale
        set(p13,'xdata',[C.hp(j,4) C.hp_r(j,4)],'ydata',[-C.hp(j,2) -C.hp_r(j,2)],...
            'zdata',[-C.hp(j,3) -C.hp_r(j,3)]); % pelvis

        set(p14,'xdata',xDis(j),'ydata',yDis(j),'zdata',zDis(j)); % CoM
        addpoints(p18,angle(j),yDis(j))
        addpoints(p19,angle(j),zDis(j))
        addpoints(p20,angle(j),xDis(j))
        addpoints(p21,xDis(j),lim,zDis(j)) % coronal plane
        addpoints(p22,lim,yDis(j),zDis(j)) % sagittal plane
        addpoints(p23,xDis(j),yDis(j),-lim) % transverse plane
        %addpoints(p24,[angle(1) angle(end)],[-com2bb(j,1) -com2bb(j,1)]) % com2bb
        %pause(0.01)
        drawnow
        frame = getframe(fig);
        writeVideo(v,frame)
    end
    close(v)
else
end
%% Check z coords
figure('Name','Check for Z-coord. dropout','color','w')
for iM = 1:numel(markerList)
    x = seconds(1/new_freq:1/new_freq:length(C.hd)/new_freq);
    plot(x,C.(markerList{iM})(:,4))
    hold on
    xlabel('Time')
    ylabel('Z coordinates (mm)')
    title('PASSED. Range of marker Z coordinates looks normal.')
    if range(ylim) >1000
        title('WARNING. Range of marker Z coordinates is higher than expected.')
        warning('Possible errors in CoM results due to marker dropout. Please inspect the plot of Z-coorindates.'); 
    end
    box off
end

% Use a standard media player to play the video. The figure should still be
% open. You can analyse the data in the figure or close it. All variables
% should still exist in the workspace.

end

