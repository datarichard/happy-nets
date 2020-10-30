load Data_USEconModel
DataTable.RGDP = DataTable.GDP./DataTable.GDPDEF*100;

%% Plot the data to look for trends

figure
subplot(3,1,1)
plot(DataTable.Time,DataTable.RGDP,'r');
title('Real GDP')
grid on
subplot(3,1,2);
plot(DataTable.Time,DataTable.M1SL,'b');
title('M1')
grid on
subplot(3,1,3);
plot(DataTable.Time,DataTable.TB3MS,'k')
title('3-mo T-bill')
grid on

%% Detrend the data
% To counter the trends in real GDP and M1, take a difference of the
% logarithms of the data. 

rgdpg = price2ret(DataTable.RGDP);
m1slg = price2ret(DataTable.M1SL);

% Also, stabilize the T-bill series by taking the 
% first difference. 

dtb3ms = diff(DataTable.TB3MS);

% Synchronize the date series so that the data has the 
% same number of rows for each column.

Data = array2timetable([rgdpg m1slg dtb3ms],...
    'RowTimes',DataTable.Time(2:end),...
    'VariableNames',{'RGDP' 'M1SL' 'TB3MS'});

figure
subplot(3,1,1)
plot(Data.Time,Data.RGDP,'r');
title('Real GDP')
grid on
subplot(3,1,2);
plot(Data.Time,Data.M1SL,'b');
title('M1')
grid on
subplot(3,1,3);
plot(Data.Time,Data.TB3MS,'k'),
title('3-mo T-bill')
grid on

%% Preprocess the timeseries
% The scale of the first two columns is about 100 times smaller than the
% third. Multiply the first two columns by 100 so that the time series are
% all roughly on the same scale. This scaling makes it easy to plot all the
% series on the same plot. More importantly, this type of scaling makes
% optimizations more numerically stable (for example, maximizing
% loglikelihoods).

Data{:,1:2} = 100*Data{:,1:2};
figure
plot(Data.Time,Data.RGDP,'r');
hold on
plot(Data.Time,Data.M1SL,'b');
datetick('x')
grid on
plot(Data.Time,Data.TB3MS,'k');
legend('Real GDP','M1','3-mo T-bill');
hold off

%% Clean the data
% Remove all missing values from the beginning of the series.
idx = all(~ismissing(Data),2);
Data = Data(idx,:);

%% Create four VAR models

numseries = 3;
dnan = diag(nan(numseries,1));
seriesnames = {'Real GDP','M1','3-mo T-bill'};
VAR2diag = varm('AR',{dnan dnan},'SeriesNames',seriesnames);
VAR2full = varm(numseries,2);
VAR2full.SeriesNames = seriesnames;
VAR4diag = varm('AR',{dnan dnan dnan dnan},'SeriesNames',seriesnames);
VAR4full = varm(numseries,4);
VAR4full.SeriesNames = seriesnames;

