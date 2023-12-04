%% UV-VIS SPECTRUM ANALYSIS V2

%% DESCRIPTION

% This is a script for analyzing absorbance spectrum data where you have a photocaged compound
% that is uncaged over time following UV irradiation.

% UV-Vis data is saved as an .txt file. Save files in a way that Absorbance T0 is
% alphabetically first in folder relative to the Absorbance T series. 
% Script will plot three figures; an Absorbance vs. Wavelength (nm) comparing the 
% No LED T0 scan and the T0 scan of the time series; an Absorbance vs Wavelength (nm)
% with 2 subplots, showing absorbance scan at different irradiation time points; 
% and Absorbance vs Time (min) with Tau plotted.

% written by Janelle Tobias June 2022.
% updated May 2023.

%% Matlab setup

clear all % this deletes all the variables in the workspace
close all % closes all old figures

%% EXPERIMENT INPUT PARAMETERS - ADJUST HERE

% Figure parameters
titlename = 'OCT-Cap'; % File name for saving figures
LED = '365 nm'; % LED used for uncaging
Lwidth = 1.6; % sets line thickness for plotting

% Row of wavelength value for plotting absorbanc vs time
Wave1 = 480; % row 586 = 305.7770 nm
Wave2 = 848; % row 854 = 356.5610 nm

% Time between scans
scantime = 5; % time between scans, in seconds

% Location of data in text file
RWave = 14; % row where wavelength data begins in .txt file1 where file1 is a single scan
RAbsT0 = 15; % row where absorbance data begins in .txt file1 where file1 is a single scan
RAbsSeries = 2; % row where absorbance data begins in table of file2 where file2 is a series data
CData = 3; % column where data begins in all files

% Establish estimates for exponetial curve fitting
% f(x)= a*exp(b*x) + c*exp(d*x) where a and c=coefficients, b=decay
% constant due to photoswitching, and d=decay constant not due to photoswitching
aEst = 1; % can estimate a to be the abosorbance at t = 0 min
cEst = 0; % can estimate c to be the absorbance at t = end 
bEst = -0.5*(aEst+cEst); % can estimate b to be the slope of the expontial change due to photoswitching, if decay b should be negative, if growth b should be positive
dEst = 0.01; % can estimate d to be the slope of expontial change unrelated to photoswitching

%% DATA PROCESSING - Import data and extract absorbance values

% Finds current directory folder.
CurrentFolder = dir('*.txt');  % Imports all .txt files as a structure array.

% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 3650); % Set to 3650 because that is the amount of wavelength data collected on the spectrometer

% Specify range and delimiter for Import Options
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Dataset1 - T0, No LED single Absorbance scan
file1 = readtable(CurrentFolder(1).name, opts); % creates a table including all contents of .txt file
wavelength = transpose(str2double(file1{RWave,CData:end})); % extract wavelength data
AbsDataT0 = transpose(str2double(file1{RAbsT0,CData:end})); % extract T0 absorbance scan

% Dataset2 - Tseries, multiple Absorbance scans
file2 = readtable(CurrentFolder(2).name); % creates table including all contents of text file
AbsDataTseries = transpose(file2{RAbsSeries:end,CData:end}); % extract series absorbance scans

%% DATA PROCESSING - Timestamp

dt = scantime/60; % time between scans, in minutes
[~,nscans] = size(AbsDataTseries); % # of rows in absorbacne data = # of scans 
Time = transpose(dt*((1:nscans))-dt); % creates matrix with time in minutes

%% DATA PROCESSING - Normalize data + calculate kinetics 

% Extract absorbance values at wavelength of interest
waveOI1 = wavelength(Wave1); % wavelength of interest, in nm
NormAbsW1 = transpose(AbsDataTseries(Wave1,:)/AbsDataTseries(Wave1,1)); % normalize data against T0 at Wave1

% Extract absorbance values at wavelength of interest
waveOI2 = wavelength(Wave2); % wavelength of interest, in nm
NormAbsW2 = transpose(AbsDataTseries(Wave2,:)/AbsDataTseries(Wave2,1)); % normalize data against T0 at Wave2



% Set up fittype and options.
ft = fittype('exp2');
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
% opts.StartPoint = [0.563205581488667 -1.865774081796 0.202630029403522 0.0104519550026605];
opts.StartPoint = [aEst bEst cEst dEst];

% Fit model to data.
[fitresult, gof] = fit(Time, NormAbsW2, ft, opts );
fitVariables = coeffvalues (fitresult);
Tau = abs(1/(fitVariables(2)));
TauAbs = fitVariables(1)*exp(fitVariables(2)*Tau)+fitVariables(3)*exp(fitVariables(4)*Tau);

%% Plotting Spectrum - Absorbance vs Wavelength - comparing No LED T0 vs T0 in time series

figure
hold on
title('Comparing T0 Scans');
ylabel('Absorbance (A.U.)'); % y axis label
xlabel('Wavelength (nm)'); % x axis label
hold on
plot(wavelength, AbsDataT0, 'k', 'LineWidth',Lwidth, 'DisplayName', 'No LED'); 
hold on
plot(wavelength, AbsDataTseries (:,1), 'b', 'LineWidth',Lwidth, 'DisplayName', 'T0 of time series');
axis([265 500 -0.02 0.4]); % sets axis limits
legend ('Location', 'Best'); % adds legend and puts in best location on graph
hold off
saveas(gcf,[num2str(titlename),'_',num2str(LED),'T0_AbsVsWave'],'png');

%% Plotting Spectrum - Absorbance vs Wavelength
#F44336
figure ('Position', [696,422,544,556]) % Sets size of figure
p1 = subplot (2,1,1);
hold on
set (p1, 'Position', [0.13, 0.58, 0.78, 0.34]);
ylabel('Absorbance (A.U.)'); % y axis label
xlabel('Wavelength (nm)'); % x axis label
title(['',num2str(titlename),' ',num2str(LED),' Uncaging Over 1 min']); % Title of subplot 1

hold on
plot(wavelength, AbsDataTseries (:,1),'b', 'LineWidth', Lwidth, 'DisplayName', 't=0s');
plot(wavelength, AbsDataTseries (:,2), 'LineWidth', Lwidth, 'DisplayName', 't=5s');
plot(wavelength, AbsDataTseries (:,7), 'm', 'LineWidth', Lwidth, 'DisplayName', 't=30s');
plot(wavelength, AbsDataTseries (:,10), 'LineWidth', Lwidth, 'DisplayName', 't=45s');
plot(wavelength, AbsDataTseries (:,13), 'k', 'LineWidth', Lwidth, 'DisplayName', 't=1min');
axis([265 500 -0.02 0.5]); % sets axis limits
legend ('Location', 'Best'); % adds legend and puts in best location on graph

p2 = subplot (2,1,2);
hold on
set (p2, 'Position', [0.13, 0.11, 0.78, 0.34]);
ylabel('Absorbance (A.U.)'); % y axis label
xlabel('Wavelength (nm)'); % x axis label
title(['',num2str(titlename),' ',num2str(LED),' Uncaging Over 5 min']); % Title of subplot 2

hold on
plot(wavelength, AbsDataTseries (:,1),'b', 'LineWidth', Lwidth, 'DisplayName', 't=0s');
plot(wavelength, AbsDataTseries (:,7), 'LineWidth', Lwidth, 'DisplayName', 't=30s');
plot(wavelength, AbsDataTseries (:,19), 'k',  'LineWidth', Lwidth, 'DisplayName', 't=1.5min');
plot(wavelength, AbsDataTseries (:,37),  'LineWidth', Lwidth, 'DisplayName', 't=3min');
plot(wavelength, AbsDataTseries (:,61),'r',  'LineWidth', Lwidth, 'DisplayName', 't=5min');
axis([265 500 -0.02 0.5]); % sets axis limits
legend ('Location', 'Best'); % adds legend and puts in best location on graph
hold off
saveas(gcf,[num2str(titlename),'_',num2str(LED),'TSeries_AbsVsWave'],'png');

%% Plotting Spectrum - Absorbance vs Time

figure
hold on
title(['',num2str(titlename),' Uncaging Kinetics']);
hold on
plot(Time, NormAbsW1, 'k', 'LineWidth',Lwidth, 'DisplayName', ['Normalized at ',num2str(waveOI1),' nm']); 
hold on
plot(Time, NormAbsW2/3, 'r', 'LineWidth',Lwidth, 'DisplayName', ['Normalized at ',num2str(waveOI2),' nm']); 
% hold on
% plot(fitresult); 
% set (legend, 'visible', 'off');
% hold on
% T = plot ([Tau Tau], [0 2]); % sets first compound addition line position
% T.LineWidth = avewidth; % sets first addition line width
% T.Color = 'r'; % sets first addition line color to red
% T.LineStyle = ':'; % sets first addition line style, dotted
hold on
axis([0 5 0.2 4.5]);
legend ('Location', 'Best'); % adds legend and puts in best location on graph
ylabel('Normalized Absorbance'); % y axis label
xlabel('Time (min)'); % x axis label
hold off
saveas(gcf,[num2str(titlename),'_',num2str(LED),'_Kinetics'],'png');



% carmel figures

figure ('Position', [718,137,423,359]) % Sets size of figure

plot(wavelength, AbsDataTseries (:,1), 'Color','#000000','LineWidth', Lwidth, 'DisplayName', 't=0s');
hold on
plot(wavelength, AbsDataTseries (:,2), 'Color','#611a15', 'LineWidth', Lwidth, 'DisplayName', 't=5s');
plot(wavelength, AbsDataTseries (:,7), 'Color','#c3352b', 'LineWidth', Lwidth, 'DisplayName', 't=30s');
plot(wavelength, AbsDataTseries (:,13), 'Color','#F44336','LineWidth', Lwidth, 'DisplayName', 't=1min');
plot(wavelength, AbsDataTseries (:,19), 'Color','#f9a19a', 'LineWidth', Lwidth, 'DisplayName', 't=1.5min');
plot(wavelength, AbsDataTseries (:,61), '--', 'Color','#fabdb8', 'LineWidth', Lwidth, 'DisplayName', 't=5min');

ylabel('Absorbance (A.U.)'); % y axis label
xlabel('Wavelength (nm)'); % x axis label
axis([265 400 -0.02 0.8]); % sets axis limits
box off
legend ('Location', 'Best'); % adds legend and puts in best location on graph
saveas(gcf,[num2str(titlename),'_',num2str(LED),'TSeries_AbsVsWave_final'],'emf');


a_normalized = normalize(NormAbsW2, 'range');

figure ('Position', [718,137,423,359]) % Sets size of figure
plot(Time, NormAbsW1, 'k', 'LineWidth',Lwidth, 'DisplayName', ['Normalized at ',num2str(waveOI1),' nm']); 
hold on
plot(Time, NormAbsW2, 'r', 'LineWidth',Lwidth, 'DisplayName', ['Normalized at ',num2str(waveOI2),' nm']); 

ylabel('Normalised Absorbance'); % y axis label
xlabel('Time (min)'); % x axis label
%axis([0 2 0 1]); % sets axis limits
axis square
box off
legend ('Location', 'Best'); % adds legend and puts in best location on graph
saveas(gcf,[num2str(titlename),'_',num2str(LED),'_Kinetics_final'],'png');




%% Final figure

Ncells_365=4;
averagetrace_365 = mean(alldata_365')'; % calculates average intensity at each time point of all cells from the normalized data
SEMAverage_365 = std(alldata_365')'/sqrt(Ncells_365); % calculates SEM from normalized data

Ncells_415=4;
averagetrace_415 = mean(alldata_415')'; % calculates average intensity at each time point of all cells from the normalized data
SEMAverage_415 = std(alldata_415')'/sqrt(Ncells_415); % calculates SEM from normalized data

Ncells_470=2;
averagetrace_470 = mean(alldata_470')'; % calculates average intensity at each time point of all cells from the normalized data
SEMAverage_470 = std(alldata_470')'/sqrt(Ncells_470); % calculates SEM from normalized data

Ncells_565=4;
averagetrace_565 = mean(alldata_565')'; % calculates average intensity at each time point of all cells from the normalized data
SEMAverage_565 = std(alldata_565')'/sqrt(Ncells_565); % calculates SEM from normalized data

figure
hold on
ylabel('Normalized Absorbance') % y axis label
xlabel('Time (s)') % x axis label
hold on
fill([Time;flipud(Time)],[averagetrace_365-SEMAverage_365;flipud(averagetrace_365+SEMAverage_365)],[0.8 0.7 0.8],'linestyle','none'); % adds shaded SEM around average trace
alpha(.5)
plot(Time,averagetrace_365,'Color','#611a15', 'LineWidth',2,'DisplayName', '365 nm'); % plots average trace

fill([Time;flipud(Time)],[averagetrace_415-SEMAverage_415;flipud(averagetrace_415+SEMAverage_415)],[137/255 174/255 161/255],'linestyle','none'); % adds shaded SEM around average trace
alpha(.6)
plot(Time,averagetrace_415,'Color','#53796C', 'LineWidth',2,'DisplayName', '415 nm'); % plots average trace

fill([Time;flipud(Time)],[averagetrace_470-SEMAverage_470;flipud(averagetrace_470+SEMAverage_470)],[139/255 137/255 174/255],'linestyle','none'); % adds shaded SEM around average trace
alpha(.6)
plot(Time,averagetrace_470,'Color','#4F4E70', 'LineWidth',2,'DisplayName', '470 nm'); % plots average trace

fill([Time;flipud(Time)],[averagetrace_565-SEMAverage_565;flipud(averagetrace_565+SEMAverage_565)],'k','linestyle','none'); % adds shaded SEM around average trace
alpha(.4)
plot(Time,averagetrace_565,'Color','#000000', 'LineWidth',2,'DisplayName', '565 nm'); % plots average trace
box off
saveas(gcf,'OCTCAP_kineticsUncaging_allWavelengths','png')
saveas(gcf,'OCTCAP_kineticsUncaging_allWavelengths','emf')
hold off