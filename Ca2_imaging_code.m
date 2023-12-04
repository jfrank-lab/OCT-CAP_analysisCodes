%% Ca2+ ANALYSIS - AVE / INDIV TRACES; HEAT MAP; OSCILLATION FREQUENCY; AUC; V3.0

%% DESCRIPTION

% This is a script for analyzing Ca2+ data for a photocaged compound where you have a genetically
% encoded sensor, two irradiation, and positive control addition.

% Cell data is a .csv file

% Script will plot 4 figures. First, a figure with the averaged cell Ca2+ data normalized to the average 
% of the first five frames, with error bars plotted as SEM. Second, five representative cell traces.
% Third, a heat map of the individual cell traces, normalized to the positive control. Indivial cells correspond 
% to the cell number in the heatmap data. Traces for represented cells can be picked based on heat map data.
% Fourth, a box chart plotting the number of peaks counted pre UV and post UV. 

% written by Janelle Tobias July 2023.

%% Matlab setup

clear all % this deletes all the variables in the workspace
close all % closes all old figures

%% EXPERIMENT INPUT PARAMETERS - ADJUST HERE

filename1 = 'OCT-CAP-5uM_TRPV1_SNAPmutant_3repeats_29cells_data.csv'; % Input filename for dataset 1

% File name for saving figures
titlename = 'example'; 

UV1f = 30; % frame# where first irradiation began
UVdur = 15; % number of frames for irradiation
PosContrf = 97;  % frame# when positive control was added
frametime = 4.019; % frame time in seconds

UV1f = 70; % frame# where first irradiation began
UVdur = 1; % number of frames for irradiation
PosContrf = 140;  % frame# when positive control was added
frametime = 4.019; % frame time in seconds

cell1 = 1; % extracts cell to plot for individual trace, corresponds to cell number in heat map
cell2 = 3; % extracts cell to plot for individual trace, corresponds to cell number in heat map
cell3 = 6; % extracts cell to plot for individual trace, corresponds to cell number in heat map
cell4 = 14; % extracts cell to plot for individual trace, corresponds to cell number in heat map
cell5 = 17; % extracts cell to plot for individual trace, corresponds to cell number in heat map

yPeakThreshold = 0.5; % sets threshold for finding peaks, peak must be at least this value greater than signal before / after
xPeakThreshold = 2; % finds tallest peak and eliminates peak in x direction on either side of peak

preUVf = UV1f-1; % frames of baseline
postUVf = PosContrf-UV1f; % calculated number of frames between first UV and positive control

%% APPEARANCE PARAMETERS - ADJUST HERE

avewidth = 1.6; % sets line thickness for Fig 1 trace
indwidth = 1.2;  % sets line thickness for Fig 2 traces
stimwidth = 1.5; % sets stimulation bar thickness for Fig 1/2
boxcolor = [0.65 0 1]; % sets color of light stimulation boxes: 375 nm = [0.65 0 1], Fig 1/2/4/5
MkS = 20; % defines marker size in figure 4

%% Unit Conversion - frames to time

UV1 = UV1f*frametime-frametime; % time of UV1 in seconds
UVdur = UVdur*frametime; % duration of irradiation in seconds
PosContr = PosContrf*frametime-frametime; % time of positive control addition in seconds

preUV1 = (preUVf*frametime-frametime)/60; % duration of baseline in min
postUV1 = (postUVf*frametime)/60; % converts frame to min

%% Import and Data Processing - Averaged Trace with Standard Error

alldata = csvread(filename1);

time=alldata(:,1); % reading first column for the timestamp

alldata(:,1) = []; % remove timestamp column to isolate the fluorescence intensities
[Nframe,Ncells] = size(alldata); % calculates the size of the new matrix; #rows = number of frames; #columns = number of cells

% Normalize against first 5 frames
datanorm = alldata./mean(alldata(1:5,:)); % normalizes the data against the first 5 frames
averagetrace = mean(datanorm')'; % calculates average intensity at each time point of all cells from the normalized data
totaltime=Nframe*frametime-frametime; % calculates the total length of the recording
SEMAverage = std(datanorm')'/sqrt(Ncells); % calculates SEM from normalized data

% %Normalize against KCl
% FMaxf = alldata(PosContrf-1:end,:); % extracts frames after positive control from each column for Fmax calculation
% FMax = max(FMaxf); % calculates maximum value of normalized data after positive control for each cell
% %  
% FMinf = alldata(1:UV1f-1,:); % extracts baseline frames from each column for Fmin calculation
% FMin = min(FMinf); % calculates minimum value of normalized baseline for each cell
% % 
% range = FMax-FMin;
% datanorm = alldata./range;
% averagetrace = mean(datanorm')'; % calculates average of all cells from the normalized data
% [framenumber,numberofcells] = size(datanorm); % calculates number of frames and number of cells from the normalized dataset
% totaltime=framenumber*frametime-frametime; % calculates the total length of the recording
% stderror = std(datanorm')'/sqrt(numberofcells); % calculates std error from normalized data

upperlimit = max(averagetrace); % finds maximum of averaged trace
upperlimit = max(upperlimit); % finds the largest of maximums
ymax = upperlimit*1.2; % calculates y-axis maximum for Figure 1
lowerlimit = min(averagetrace); % finds minimum of averaged trace
lowerlimit = min(lowerlimit); % finds the lowest of minimus
ymin = lowerlimit*0.8; % calculates y-axis minimum for Figure 1

%% Import and Data Processing - Individual Cell Trace

Cell1 = datanorm (:,cell1); % extracts column from normalized data
Cell2 = datanorm (:,cell2); % extracts column from normalized data
Cell3 = datanorm (:,cell3); % extracts column from normalized data
Cell4 = datanorm (:,cell4); % extracts column from normalized data
Cell5 = datanorm (:,cell5); % extracts column from normalized data

indivcells = [Cell1 Cell2 Cell3 Cell4 Cell5]; % combines individual traces into a single array
indupperlimit = max(indivcells); % finds maximum of each individual trace
indupperlimit = max(indupperlimit); % finds the largest of the maximums
ymaxind = indupperlimit*1.05; % sets y-axis maximum for Figure 2
indlowerlimit = min(indivcells); % finds minimum of each individual traces
indlowerlimit = min(indlowerlimit); % finds the lowest of the minimums
yminind = indlowerlimit*0.95; % sets y-axis minimum for Figure 2

%% Heat Map - Normalize Against Positive Control

FnormMaxf = alldata(PosContrf-1:end,:); % extracts frames after positive control from each column for Fmax calculation
FnormMax = max(FnormMaxf); % calculates maximum value of normalized data after positive control for each cell
 
FnormMinf = alldata(1:preUVf,:); % extracts baseline frames from each column for Fmin calculation
FnormMin = min(FnormMinf); % calculates minimum value of normalized baseline for each cell

range = FnormMax-FnormMin; % calculates the range in intensity values relative to the max value after the positive control and min value from the baseline
datanormPosContr = alldata./range; % calculates each intensity value as a fraction of the range
transdatanorm = datanormPosContr.'; % Transform data to switch axes for heat map, rows are individual cells and columns are normalized intensities over time

%% PLOTTING FIGURE 1 - averaged trace with error bars

figure
hold on
PosContr1 = plot ([PosContr PosContr], [ymin ymax]); % sets first compound addition line position
PosContr1.LineWidth = stimwidth; % sets first addition line width
PosContr1.Color = 'r'; % sets first addition line color to red
PosContr1.LineStyle = ':'; % sets first addition line style, dotted
hold on
scaleBarX = plot([25 85], [1.5 1.5]);
scaleBarX.LineWidth = 2; % sets first addition line width
scaleBarX.Color = 'k'; % sets first addition line color to black
scaleBarY = plot([25 25], [1.5 1]);
scaleBarY.LineWidth = 2; % sets first addition line width
scaleBarY.Color = 'k'; % sets first addition line color to black
ylabel('F/Fmin') % y axis label
xlabel('Time (s)') % x axis label
title(['Average (N = ',num2str(Ncells),' cells)'])
rectangle('Position',[UV1 ymin UVdur ymax],'FaceColor',boxcolor,'EdgeColor',boxcolor) % first UV stimulation marker
hold on
fill([time;flipud(time)],[averagetrace-SEMAverage;flipud(averagetrace+SEMAverage)],'k','linestyle','none'); % adds shaded SEM around average trace
alpha(.4)
plot(time,averagetrace,'k', 'LineWidth',avewidth); % plots average trace
axis([0 totaltime ymin ymax]) % sets axes limits
box off
saveas(gcf,[num2str(titlename),'_AverageTrace_SB_60sec_0-5FFmin'],'png')
saveas(gcf,[num2str(titlename),'_AverageTrace_SB_60sec_0-5FFmin'],'emf')
hold off


%% PLOTTING FIGURE 2 - individual cell traces for 5 cells

figure
subplot(1,3,2) %sets up plot space as an mxn grid, placing graphs in each position (m,n,position)
hold 
g1=subplot(5,1,1); % sets up the first plot of a 5x1 plot
hold on
PosContr21 = plot ([PosContr PosContr], [yminind ymaxind]); % sets first compound addition line position
PosContr21.LineWidth = stimwidth; % sets first addition line width
PosContr21.Color = 'r'; % sets first addition line color to red
PosContr21.LineStyle = ':'; % sets first addition line style, dotted
hold on
rectangle('Position',[UV1 yminind UVdur ymaxind],'FaceColor',boxcolor,'EdgeColor',boxcolor) % first UV stimulation marker
plot(time,Cell1,'k','LineWidth',indwidth); % plots trace
axis([0 totaltime yminind ymaxind]) % sets axis
ylabel('F/Fmin') % y label
%xlabel('Time (s)') % x label
title('Cell #1') % title
box off
hold off

hold
g2=subplot(5,1,2); % sets up the second plot of a 4x1 plot
hold on
PosContr22 = plot ([PosContr PosContr], [yminind ymaxind]); % sets first compound addition line position
PosContr22.LineWidth = stimwidth; % sets first addition line width
PosContr22.Color = 'r'; % sets first addition line color to red
PosContr22.LineStyle = ':'; % sets first addition line style, dotted
hold on
rectangle('Position',[UV1 yminind UVdur ymaxind],'FaceColor',boxcolor,'EdgeColor',boxcolor) % first UV stimulation marker
plot(time,Cell2,'k','LineWidth',indwidth); %plots trace
axis([0 totaltime yminind ymaxind]) % sets axis
ylabel('F/Fmin') % y label
%xlabel('Time (s)') % x label
title('Cell #2') % title
box off
hold off

hold
g3=subplot(5,1,3); % sets up the third plot of a 4x1 plot
hold on
PosContr23 = plot ([PosContr PosContr], [yminind ymaxind]); % sets first compound addition line position
PosContr23.LineWidth = stimwidth; % sets first addition line width
PosContr23.Color = 'r'; % sets first addition line color to red
PosContr23.LineStyle = ':'; % sets first addition line style, dotted
hold on
rectangle('Position',[UV1 yminind UVdur ymaxind],'FaceColor',boxcolor,'EdgeColor',boxcolor) % first UV stimulation marker
plot(time,Cell3,'k','LineWidth',indwidth); % plots trace
axis([0 totaltime yminind ymaxind]) % sets axis
ylabel('F/Fmin') % y label
%xlabel('Time (s)') % x label
title('Cell #3') % title
box off
hold off

hold
g4=subplot(5,1,4); % sets up the fourth plot of a 4x1 plot
hold on
PosContr24 = plot ([PosContr PosContr], [yminind ymaxind]); % sets first compound addition line position
PosContr24.LineWidth = stimwidth; % sets first addition line width
PosContr24.Color = 'r'; % sets first addition line color to red
PosContr24.LineStyle = ':'; % sets first addition line style, dotted
hold on
rectangle('Position',[UV1 yminind UVdur ymaxind],'FaceColor',boxcolor,'EdgeColor',boxcolor) % first UV stimulation marker
plot(time,Cell4,'k','LineWidth',indwidth); % plots trace
axis([0 totaltime yminind ymaxind]) % sets axis
ylabel('F/Fmin') % y label
%xlabel('Time (s)') % x label
title('Cell #4') % title
box off
hold off

hold
g5=subplot(5,1,5); % sets up the fourth plot of a 4x1 plot
hold on
PosContr25 = plot ([PosContr PosContr], [yminind ymaxind]); % sets first compound addition line position
PosContr25.LineWidth = stimwidth; % sets first addition line width
PosContr25.Color = 'r'; % sets first addition line color to red
PosContr25.LineStyle = ':'; % sets first addition line style, dotted
hold on
rectangle('Position',[UV1 yminind UVdur ymaxind],'FaceColor',boxcolor,'EdgeColor',boxcolor) % first UV stimulation marker
plot(time,Cell5,'k','LineWidth',indwidth); % plots trace
axis([0 totaltime yminind ymaxind]) % sets axis
ylabel('F/Fmin') % y label
xlabel('Time (s)') % x label
title('Cell #5') % title
box off
saveas(gcf,[num2str(titlename),'_IndivTrace'],'png')
hold off

%% PLOTTING FIGURE 3 - Heat Map

figure
H = heatmap (transdatanorm,'Colormap',parula); % generates a heat map of transformed, normalized data
H.Title = (['Normalized Against KCl (N = ',num2str(Ncells),' Cells)']);
H.GridVisible = 'off';
saveas(gcf,[num2str(titlename),'_HeatMap'],'png')


