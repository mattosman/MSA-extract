function [MSA] = processMSAarray(maxYr,minYr,maxStdYr,minStdYr,smoother,group,plotOpt)
%% Inputs and Outputs
% 
% Written by M. Osman (osmanm@mit.edu)
% 
% Inputs:
%   maxYr = maximum (most recent) year to plot
%   minYr = minimum (oldest) year to plot
%   maxStdYr = standardization range; youngest (most recent) year
%   minStdYr = standardization range; oldest year
%   smoother = low pass filter data to include frequencies of 1/smoother and smaller (where smoother is in units of years)
%   group = choose either 'Greenland'', ''NorthernGreenland'', or ''SouthernGreenland''
%   plotOpt = if == 1 plots, if = otherwise, don't plot 
% Outputs: 
%   year = year range to plot against, chosen as range between max_year and min_year inputs; of size M, where M is the number of years
%   data_std = data incorporated into plots - the data has been standardized over the year range pertaining to std_max_year and std_min_year; of size M x N, where N is the number of sites
%   data_lp = as in data_std, but low pass filtered to include frequencies of 1/smoother and smaller, of size M x N
%   colNames = names of the sites heading the data_std and data_lp columns, of size N
% 
% To run, requires the following to be in the Current Folder:
%   1. lowpass.m
%   2. GrIS_MSA_recs.xlsx

%%
if nargin<7; plotOpt = []; end;
    if isempty(plotOpt)
        plotOpt = 0;
    end
%%

% must have 'MSA_greenland_array_repository.xlsx' in current folder!
[data, ~] = xlsread('GrIS_MSA_recs.xlsx','Sheet1');

years = data(:,1);
    Iplot = years <= maxYr & years >= minYr; 
    years = years(Iplot);
    Istd = years <= maxStdYr & years >= minStdYr; 
    
% standardize_data
% Greenland sites
Summit2010_std = data(Iplot,3);  Summit2010_std = (Summit2010_std - nanmean(Summit2010_std(Istd)))/nanstd(Summit2010_std(Istd));
D4_std         = data(Iplot,4);  D4_std         = (D4_std - nanmean(D4_std(Istd)))/nanstd(D4_std(Istd));
GC_std         = data(Iplot,5);  GC_std         = (GC_std - nanmean(GC_std(Istd)))/nanstd(GC_std(Istd));
D20_std        = data(Iplot,6);  D20_std        = (D20_std - nanmean(D20_std(Istd)))/nanstd(D20_std(Istd));
GRIP93a_std    = data(Iplot,7);  GRIP93a_std    = (GRIP93a_std - nanmean(GRIP93a_std(Istd)))/nanstd(GRIP93a_std(Istd));
NGRIP_std      = data(Iplot,8);  NGRIP_std      = (NGRIP_std - nanmean(NGRIP_std(Istd)))/nanstd(NGRIP_std(Istd));
B16_std        = data(Iplot,9);  B16_std        = (B16_std - nanmean(B16_std(Istd)))/nanstd(B16_std(Istd));
TUNU_std       = data(Iplot,2);  TUNU_std       = (TUNU_std - nanmean(TUNU_std(Istd)))/nanstd(TUNU_std(Istd));
B18_std        = data(Iplot,10); B18_std        = (B18_std - nanmean(B18_std(Istd)))/nanstd(B18_std(Istd));
B20_std        = data(Iplot,11); B20_std        = (B20_std - nanmean(B20_std(Istd)))/nanstd(B20_std(Istd));
B21_std        = data(Iplot,12); B21_std        = (B21_std - nanmean(B21_std(Istd)))/nanstd(B21_std(Istd));
B26_std        = data(Iplot,13); B26_std        = (B26_std - nanmean(B26_std(Istd)))/nanstd(B26_std(Istd));

TF_G =         strcmp(group,'Greenland');
TF_NG = strcmp(group,'NorthernGreenland');
TF_SG = strcmp(group,'SouthernGreenland');

    if TF_G == 1
        zDat =  horzcat(Summit2010_std,D4_std,GC_std,D20_std,GRIP93a_std,NGRIP_std,B16_std,TUNU_std,B18_std,B20_std,B21_std,B26_std);
        colNames = {'Summit2010','D4','GC','20D','GRIP93a','NGRIP','NGT-B16','TUNU','NGT-B18','NGT-B20','NGT-B21','NGT-B26'};
    elseif TF_NG == 1
        zDat =  horzcat(TUNU_std,B18_std,B20_std,B21_std, B26_std); % B21_std, B26_std
        colNames = {'TUNU','NGT-B18','NGT-B20','NGT-B21','NGT-B26'}; % ,'NGT-B21','NGT-B26'}
    elseif TF_SG == 1
        zDat =  horzcat(Summit2010_std,D4_std,GC_std,D20_std,GRIP93a_std,NGRIP_std,B16_std);
        colNames = {'Summit2010','D4','GC','20D','GRIP93a','NGRIP','NGT-B16'};
    else
        disp(['Input a group to plot! Choose either ''Greenland'', ''NorthernGreenland'', or ''SouthernGreenland'''])
        return
    end

% how many records vs. time? 
for k = 1:length(zDat(:,1))
    numRecs(k,1) = sum(~isnan(zDat(k,:)));
end
Ione = numRecs == 1;

% take stack mean
for k = 1:length(zDat(:,1))
    muDat(k,1) = nanmean(zDat(k,:));
    sigDat(k,1) = nanstd(zDat(k,:));
end
muDat(Ione) = NaN;
sigDat(Ione) = NaN;

if plotOpt == 1
    
    % plot the group
    fig1 = figure; 

    subplot(3,6,[2:6, 8:12])
    for i = 1:length(zDat(1,:))
        hold on
        plot(years,zDat(:,i),'-');% ,'color',[0.65 0.65 0.65],'linewidth',1);
        ylabel('MS^{-} (z-score)')
        xlim([min(years) max(years)])
        ylim([-4 4])
        legendInfo{i} = [colNames{i}]; % or whatever is appropriate
        set(gca,'Fontsize',14,'Linewidth',2)
        hold off
    end
    % hold on; plot(years,mean_data_std,'-','linewidth',2.5,'Color',[0 0 0]); hold off;
    box on; grid on; legend(legendInfo)

    subplot(3,6,[14:18])
        plot(years,numRecs,'-k','linewidth',3,'Color',[0.6 0.6 0.6]); 
        ylabel('Number of Records')
        xlim([min(years) max(years)])
        ylim([0 max(numRecs)+1])
        xlabel('Year (A.D.)')
        set(gca,'Fontsize',14,'Linewidth',2)
        box on; grid on;

    set(fig1,'PaperPositionMode','auto');         
    set(fig1,'PaperOrientation','landscape');
    set(fig1,'Position',[50 50 1000 500]); 

end

    %% low pass filter
    
if abs(maxYr - minYr) >= 30
    
    % Greenland sites
    ind_Summit2010 = ~isnan(Summit2010_std);
    ind_D4         = ~isnan(D4_std);
    ind_GC         = ~isnan(GC_std);
    ind_D20        = ~isnan(D20_std);
    ind_GRIP93a    = ~isnan(GRIP93a_std);
    ind_NGRIP      = ~isnan(NGRIP_std);
    ind_B16        = ~isnan(B16_std);
    ind_TUNU       = ~isnan(TUNU_std);
    ind_B18        = ~isnan(B18_std);
    ind_B20        = ~isnan(B20_std);
    ind_B21        = ~isnan(B21_std);
    ind_B26        = ~isnan(B26_std);
    % Stack mean
    ind_nums_stack_mean = ~isnan(muDat);
    ind_nums_stack_std  = ~isnan(sigDat);

        % this is just in case only on record exists in a region
        if sum(~isnan(muDat)) < 1 % if no non-NaN values exist...
            for k = 1:length(zDat(:,1))
                muDat(k,1) = nanmean(zDat(k,:)); % recalculate the mean_data_std variable
                sigDat(k,1) = nanstd(zDat(k,:)); % recalculate the mean_data_std variable
            end
            ind_nums_stack_mean   = ~isnan(muDat);
            ind_nums_stack_std    = ~isnan(sigDat);
        end


    % Greenland sites
    [Summit2010_lp,~] = lowpass(Summit2010_std(~isnan(Summit2010_std)),1/smoother,1,1);
    [D4_lp,~]         = lowpass(D4_std(~isnan(D4_std)),1/smoother,1,1);
    [GC_lp,~]         = lowpass(GC_std(~isnan(GC_std)),1/smoother,1,1);
    [D20_lp,~]        = lowpass(D20_std(~isnan(D20_std)),1/smoother,1,1);
    [GRIP93a_lp,~]    = lowpass(GRIP93a_std(~isnan(GRIP93a_std)),1/smoother,1,1);
    [NGRIP_lp,~]      = lowpass(NGRIP_std(~isnan(NGRIP_std)),1/smoother,1,1);
    [B16_lp,~]        = lowpass(B16_std(~isnan(B16_std)),1/smoother,1,1);
    [TUNU_lp,~]       = lowpass(TUNU_std(~isnan(TUNU_std)),1/smoother,1,1);
    [B18_lp,~]        = lowpass(B18_std(~isnan(B18_std)),1/smoother,1,1);
    [B20_lp,~]        = lowpass(B20_std(~isnan(B20_std)),1/smoother,1,1);
    [B21_lp,~]        = lowpass(B21_std(~isnan(B21_std)),1/smoother,1,1);
    [B26_lp,~]        = lowpass(B26_std(~isnan(B26_std)),1/smoother,1,1);
    % low pass the stack mean
    [muStack_lp,~] = lowpass(muDat(~isnan(muDat)),1/smoother,1,1);
    [sigStack_lp,~]  = lowpass(sigDat(~isnan(sigDat)),1/smoother,1,1);

    % Greenland sites
    Summit2010_LP = nan(size(years));    Summit2010_LP(ind_Summit2010) = Summit2010_lp;
    D4_LP         = nan(size(years));    D4_LP(ind_D4)                 = D4_lp;
    GC_LP         = nan(size(years));    GC_LP(ind_GC)                 = GC_lp;
    D20_LP        = nan(size(years));    D20_LP(ind_D20)               = D20_lp;
    GRIP93a_LP    = nan(size(years));    GRIP93a_LP(ind_GRIP93a)       = GRIP93a_lp;
    NGRIP_LP      = nan(size(years));    NGRIP_LP(ind_NGRIP)           = NGRIP_lp;
    B16_LP        = nan(size(years));    B16_LP(ind_B16)               = B16_lp;
    TUNU_LP       = nan(size(years));    TUNU_LP(ind_TUNU)             = TUNU_lp;
    B18_LP        = nan(size(years));    B18_LP(ind_B18)               = B18_lp;
    B20_LP        = nan(size(years));    B20_LP(ind_B20)               = B20_lp;
    B21_LP        = nan(size(years));    B21_LP(ind_B21)               = B21_lp;
    B26_LP        = nan(size(years));    B26_LP(ind_B26)               = B26_lp;
    % Stack mean
    muStack_LP = nan(size(years));    muStack_LP(ind_nums_stack_mean) = muStack_lp;
    sigStack_LP  = nan(size(years));    sigStack_LP(ind_nums_stack_std)   = sigStack_lp;

    % concatenate all datasets
        if TF_G == 1
            data_lp =  horzcat(Summit2010_LP,D4_LP,GC_LP,D20_LP,GRIP93a_LP,NGRIP_LP,B16_LP,TUNU_LP,B18_LP,B20_LP,B21_LP,B26_LP);
            colNames = {'Summit2010','D4','GC','20D','GRIP93a','NGRIP','NGT-B16','TUNU','NGT-B18','NGT-B20','NGT_B21','NGT_B26'};
        elseif TF_NG == 1
            data_lp = horzcat(TUNU_LP,B18_LP,B20_LP,B21_LP,B26_LP); % B21_LP,B26_LP
            colNames = {'TUNU','NGT-B18','NGT-B20','NGT-B21','NGT-B26'}; % ,'NGT-B21','NGT-B26'
        elseif TF_SG == 1
            data_lp = horzcat(Summit2010_LP,D4_LP,GC_LP,D20_LP,GRIP93a_LP,NGRIP_LP,B16_LP);
            colNames = {'Summit2010','D4','GC','20D','GRIP93a','NGRIP','NGT-B16'};
        else
            disp('Input a group to plot! Input either ''All'', ''NorthernGreenland'', ''SouthernGreenland''.')
            return
        end

    if plotOpt == 1
        % Plot
        fig2 = figure; 

        subplot(3,6,[2:6, 8:12])
        for i = 1:length(data_lp(1,:))
            hold on
%             if i < 8
                plot(years,(data_lp(:,i)),'-');% ,'Color',[0.65 0.65 0.65],'linewidth',1.5); % svalbard 0.6 0.2 0.8 % canadian arctic 0.2 0 1 % northern greenland 1 0.6 0 % southern grenland 1 0.4 0.2
            ylabel('MS^{-} (z-score)')
            xlim([min(years) max(years)])
            ylim([-3 3])
            legendInfo{i} = [colNames{i}]; % or whatever is appropriate
            xlabel('Year (A.D.)')
            set(gca,'Fontsize',14,'Linewidth',2)
            hold off
        end
        % hold on; plot(years,Stack_LP_mean,'-','linewidth',2.5,'Color',[0 0 0]); hold off; % svalbard 0.4 0 0.6 % northern greenland 1 0.4 0
        box on; grid on; legend(legendInfo)

        subplot(3,6,[14:18])
            plot(years,numRecs,'-k','linewidth',3,'Color',[0.6 0.6 0.6]); 
            ylabel('Number of Records')
            ylim([0 max(numRecs)+1])
            xlim([min(years) max(years)])
            xlabel('Year (A.D.)')
            set(gca,'Fontsize',14,'Linewidth',2)
            box on; grid on;

        set(fig2,'PaperPositionMode','auto');         
        set(fig2,'PaperOrientation','landscape');
        set(fig2,'Position',[50 50 1000 500]); 

        % plot records vs. time
    end
    
    else
        warning('Warning! The assigned time-range must represent >= 30 yrs of data to lowpass filter the data.  Skipping lowpass calculation.')
        data_lp = [];
        muStack_LP  = []; 
        sigStack_LP   = []; 
end

% compute correlation matrix
for i = 1:size(zDat,2)
    for j = i:size(zDat,2)
        ind_overlapping_vals = ~isnan(zDat(:,i).*zDat(:,j));
        [r, p] = corrcoef(zDat(ind_overlapping_vals,i),zDat(ind_overlapping_vals,j));
        corrTable.correlation(i,j) = r(2,1);  % CorrTable.correlation(j,i) = CorrTable.correlation(i,j);
        corrTable.significance(i,j) = p(2,1); % CorrTable.significance(j,i) = CorrTable.significance(i,j);
        corrTable.deg_freedom(i,j) = sum(ind_overlapping_vals) - 1;% CorrTable.deg_freedom(j,i) = CorrTable.deg_freedom(i,j);
    end
end

MSA.year              = years;
MSA.data              = zDat; 
MSA.stack_data        = muDat; 
MSA.data_lp           = data_lp; 
MSA.stack_data_lp     = muStack_LP;
MSA.stack_data_lp_std = sigStack_LP;
MSA.number_of_records = numRecs; 
MSA.colNames          = colNames;
MSA.stdev_data_std    = sigDat;
MSA.CorrTable         = corrTable;

end