%% Specify the folder where the files live.
myFolder = '‪C:\Users\gaopa\Documents\Ongoing';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end

%% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*_MembraneTest.txt'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
% If want to get a list of all files in the folder, and its subfolders, with the desired file name pattern.
% filePattern = fullfile(myFolder, '**/*.txt'); % Change to whatever pattern you need.

% preallocate a matrix
neurons = 100;
collection = zeros(600,neurons);

for theFiles_i = 1 : length(theFiles)
    baseFileName = theFiles(theFiles_i).name;
    fullFileName = fullfile(theFiles(theFiles_i).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    % such as reading it in as an array with readtable and table2array
    % sgl_table = readtable(fullFileName,'PreserveVariableNames',true);
    sgl_table = readtable(fullFileName);% substitute this sentence above used for old version of matlab
    collection(:,theFiles_i) = table2array(sgl_table);
end

%% Intrinsic properties
sample_rate_kHz = 20;
T = 1/sample_rate_kHz;

% trim preallocated matrix 
collection = collection(:,1:theFiles_i); 

% loop the analysis for different cells
[~,collect_Size] = size(collection); % column size
for collect_i = 1 : collect_Size
    baseFileName = theFiles(collect_i).name;
    [Cm_fast,Cm_slow,Rm,Rs,Rsquare,Rsquare_adj,RMSE,tau_fast,tau_slow] = MembraneTest_fun(collection(:,collect_i),baseFileName,T);
    theFiles(collect_i).Cm_fast = Cm_fast;
    theFiles(collect_i).Cm_slow = Cm_slow;
    theFiles(collect_i).Rm = Rm;
    theFiles(collect_i).Rs = Rs;
    theFiles(collect_i).Rsquare = Rsquare;
    theFiles(collect_i).Rsquare_adj = Rsquare_adj;
    theFiles(collect_i).RMSE = RMSE;
    theFiles(collect_i).tau_fast = tau_fast;
    theFiles(collect_i).tau_slow = tau_slow;
end