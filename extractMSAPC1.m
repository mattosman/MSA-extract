%% Data infilling using EOF's;
% Script written by M. Osman (MIT/WHOI; osmanm@mit.edu); Aug. 2017
% Extracts MSA-PC1 from the GrIS MSA ice core array, as described in Osman et al., 2019 (Nature)
% To run, requires the following to be in the Current Folder:
%   1. processMSAarray.m
%   2. eof_infill.m
%   3. lowpass.m
%   4. GrIS_MSA_recs.xlsx

clear all
PC_mode = 1; % retain only the first mode of variability

% first , infill the end of the MSA record, prior to 1985 (the youngest 20D year)
[msa] = processMSAarray(2013,1767,1985,1821,10,'Greenland'); 
age = msa.year;

% define missing data and the data variance
missing_data  = isnan(msa.data);
data_variance = nanvar(msa.data);

% number of cross-validation procedures;
n_xval = 100; % note: Osman et al. (2019) run 1000 realizations
% preallocate output
    data = cell(n_xval,1); 
    rmse = cell(n_xval,1);  min_rmse = zeros(size(rmse,1),1);
    nEOF = cell(n_xval,1); max_nEOF = zeros(size(rmse,1),1);
    r2_cv = cell(n_xval,1); max_r2_cv = zeros(size(rmse,1),1);
    noiseDat = cell(n_xval,1); 
    pc_orig = cell(n_xval,1);  % save as cell to store the entire PC_orig matrix for subs. bs adjustments
    pve_orig = zeros(n_xval,1);
    
for i = 1:n_xval
    
    while max_nEOF(i,1) <= 2
        [data{i}, rmse{i}, nEOF{i}, ~, noiseDat{i}, r2_cv{i}] = eof_infill(msa.data,1e-5);
        min_rmse(i,1)   = min(rmse{i});    
        max_nEOF(i,1) = max(nEOF{i});
        max_r2_cv(i,1) = abs(max(r2_cv{i}));
    end
    
    % now, add remaining nvariance back into the imputed values as autocorrelated noise
    autocorr_noise = zeros(size(noiseDat{i}));
    for j = 1:size(noiseDat{i},2)
        [autoCov, lags] = xcov(noiseDat{i}(:,j),'coef');
         autoCov = autoCov(lags >= 0);
         autoCov = toeplitz(autoCov);
         [~,pos_def] = chol(autoCov);
            while pos_def > 0 % if near-singular performs slightly crude "just in case" regularization
            	autoCov = autoCov + 1e-9.*(eye(size(autoCov)));
                [~,pos_def] = chol(autoCov);
            end
            autocorr_noise(:,j) = (autoCov')^0.5 * randn(size(noiseDat{i}(:,1),1),1);
            % scale by the remaining variance of the noise
            autocorr_noise(:,j) = var(noiseDat{i}(:,j)) .* autocorr_noise(:,j);
    end
    % add the noise back into the imputed data:
    data{i}(missing_data) = data{i}(missing_data) + autocorr_noise(missing_data);
        
    % compute the Principal components for each, and retain only the denoted "mode"
    [U, S, V] = svd(data{i});
    pc_orig{i} = (V' * data{i}')' ; 
    pve = diag(S).^2/sum(diag(S).^2);
        pve_orig(i) = max(pve);
    
    % require the orientation to be consistent between each PC; 
    if i == 1 % don't assign during the first iteration; this will be the PC1 reference orientation
    else
        check_orientation = corrcoef(pc_orig{1}(:,PC_mode),pc_orig{i}(:,PC_mode));
        if check_orientation(2,1) < 0
            pc_orig{i}(:,PC_mode) = -1.*pc_orig{i}(:,PC_mode);
        end
    end
    
end
clearvars U S V PC PVE missing_data W n D autocorr autocorr_noise pos_def data_variance check_orientation

% now run the bootstrap N times for each n_xval iteration.
if n_xval <= 100
    N = n_xval; % 40*size(MSA.data,2); % suggested number of bootstrap iterations
else
    N = n_xval/10;
end

pc_boot = cell(n_xval,1);
pve_boot = zeros(n_xval,N);
boot_index = cell(n_xval,1); % store which timeseries are called just for reference
boot_data = zeros(size(msa.data)); % hold the random bootstrap series for each iteration
% bootstrap loop
for j = 1:n_xval
    disp(['Bootstrap iteration: ',num2str(j)]);
    for i = 1:N
        
        boot_index{j}(i,:) = randi([1, size(msa.data,2)], 1, size(msa.data,2));
        for g = 1:size(msa.data,2)
            boot_data(:,g) = data{j}(:,boot_index{j}(i,g)); 
        end
    
        [U, S, V] = svd(boot_data);
        pc = (V' * boot_data')' ; 
        pve = diag(S).^2/sum(diag(S).^2);
            pve_boot(j,i) = max(pve);
            
    % rotate + save PC output
    [~, Z] = procrustes(pc_orig{j},pc);
    pc_boot{j}(:,i) =  Z(:,PC_mode);
    end
end
clearvars U S V PC PVE boot_data Z

% re-center and scale PC's to full period;
for j = 1:n_xval
    for i = 1:N
        pc_boot{j}(:,i) = zscore(pc_boot{j}(:,i)); 
    end
end

% To check, compute the normal non-bias corrected 95% confidence limits from bootstrap tests
% resize the PC_boot cell array into a matrix
for i = 1:length(pc_boot)
    if i == 1
        pc_bootstrap = pc_boot{i};
    else
        pc_bootstrap = horzcat(pc_bootstrap,pc_boot{i});
    end
end

% Finally! Extract MSA-PC1
low = zeros(size(msa.data,1),1);
med = zeros(size(msa.data,1),1);
upp = zeros(size(msa.data,1),1);
sig = zeros(size(msa.data,1),1);
for i = 1:size(msa.data,1)
    low(i,1) = prctile(pc_bootstrap(i,:),2.5);
    med(i,1) = prctile(pc_bootstrap(i,:),50);
    upp(i,1) = prctile(pc_bootstrap(i,:),97.5);
    sig(i,1) = nanstd(pc_bootstrap(i,:));
end

%% Plot it up

h = figure; hold on;
	set(h,'units','centimeters','position',[1,1,25,10]);
	ax = gca; ax.Visible = 'off';
	set(h,'PaperPositionMode' ,'auto');         
	set(h,'PaperOrientation','landscape');
	set(h,'Color',[1 1 1]); 

% [MSA]-PC1
    ax1 = axes('Position',[0.20 0.25 0.43 0.60]); hold on; box on; grid off;
    h1 = fill([age(1);  age;  flipud([(age);  age(end)])],...
    [upp(1);  smooth(low,1);  flipud([smooth(upp,1);  low(end)])],[0.5 0.5 0.5]);
    set(h1,'edgecolor','none','facealpha',0.6);
    p1 = plot(age, med,'-','color',[0.1 0.1 0.2],'linewidth',0.7); p1.Color(4) = 0.75;
    p2 = plot(age, lowpass(med,1/10,1,1),'-','color',[1 0.1 0],'linewidth',1.8); p2.Color(4) = 0.95;
    set(gca,'xlim',[min(age) max(age)],'ylim',[-2.3 2.3])
    set(gca,'Color','none','LineWidth',1.5,'Fontsize',12)
    xlabel('Year (A.D.)'); ylabel('[MSA]-PC1 ({\itz})'); 
    
% Percent variance explained
    % calculate prob density function
    pve = 100.*pve_boot(:);
    xbin = 0:2:100;
    [n, ~] = histc(pve, xbin);
    n_norm = n./trapz(xbin,n);
    
    ax2 = axes('Position',[0.65 0.25 0.15 0.60]); hold on; box on; grid off;
    b = bar(xbin, n_norm, 'hist');
    set(b,'FaceColor',[0.7 .7 .7],'Edgecolor',[0.3 0.3 0.3],'LineWidth',1.5); hold on;   
    set(gca,'xlim',[20 70],'ylim',[0 0.07],'xtick',[30 50 70])
    set(gca,'YAxisLocation','right','Color','none','LineWidth',1.5,'Fontsize',12) 
    xlabel('Var. Expl. (%)'); ylabel('PDF (%^{-1})'); 
