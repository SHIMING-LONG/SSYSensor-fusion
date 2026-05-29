close all; clear; clc;
addpath(genpath('./functions'));

m = mobiledev('iPhone - smeroree🐟');

[calAcc, calGyr, calMag] = calibrateiOS(m);

analyzeData(calAcc, calGyr, calMag);

save('calibration.mat', 'calAcc', 'calGyr', 'calMag');
