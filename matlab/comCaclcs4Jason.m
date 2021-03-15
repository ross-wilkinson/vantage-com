%% Running vantage_com function for Jason and Todd

% Ross change this
mydir = '/Users/rosswilkinson/Google Drive/projects/vantage-com/data/COM Pose files/Jason Pose files';
resDir = '/Users/rosswilkinson/Google Drive/projects/vantage-com/results';
cd(mydir)
myDir = dir('*.pose');
% myDir(1) = []; myDir(1) = [];
mass = 91.6;
myComList = struct();
Hip = 370;
Shoulder = 450;

for i=1:length(myDir)
%     myPose = fullfile(myDir(i).folder,myDir(i).name);
    cd(mydir)
    myPose = myDir(i).name;
    [com3d, com2bb, com2bbMean] = vantage_com(myPose,'off','male', mass, 54, Shoulder, Hip);
    myComList(i).rider = 'Jason';
    myComList(i).poseFile = myDir(i).name;
    myComList(i).com2bbMean = com2bbMean;
end

myTable = struct2table(myComList);
cd(resDir)
writetable(myTable,'com2bbJason.csv');

%% FOR TODD
%Ross change this
mydir = '/Users/rosswilkinson/Google Drive/projects/vantage-com/data/COM Pose files/Todd Pose files';
resDir = '/Users/rosswilkinson/Google Drive/projects/vantage-com/results';
cd(mydir)
myDir = dir('*.pose');
% myDir(1) = []; myDir(1) = [];
myComList = struct();
mass = 78.8;
Hip = 296;
Shoulder = 411;

for i=1:length(myDir)
%     myPose = fullfile(myDir(i).folder,myDir(i).name);
    cd(mydir)
    myPose = myDir(i).name;
    [com3d, com2bb, com2bbMean] = vantage_com(myPose,'off','male', mass, 54, Shoulder, Hip);
    myComList(i).rider = 'Todd';
    myComList(i).poseFile = myDir(i).name;
    myComList(i).com2bbMean = com2bbMean;
end

myTable = struct2table(myComList);
cd(resDir)
writetable(myTable,'com2bbTodd.csv');

