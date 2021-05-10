function [Xa, rmse, neof, Xa_stack, Xa_noise, r2_cv] = eof_infill(Xo,delrms,verbal)

%% Data infilling using EOF's;
% Script written by M. Osman (MIT/WHOI; osmanm@mit.edu) Aug. 2017 as described in : 
% Beckers and Rixen, "EOF Calculations and Data Filling from  Incomplete Oceanographic Datasets",
% J. of Atm. and Oc. Technology, 20, 1839-1856, 2003., and as motivated by discussion 
% in Kinnard et al., 2011 (Nature) SM
% 
% Inputs:
% Xo is an MxN data field with missing values set as NaN
% delrms is the convergence threshhold; suggested as 1e-5 (rms = root mean square difference)
% verbal is an option to output iterations, must put either "true" or "false"; defaults to false

% Outputs:
% Xa = the infilled data product
% rmse = the root mean squared error of the infilled data product (a vector
%   denoting number of iterations required to reach convergence)
% neof = the number of EOFs required to reach convergence (a vector
%   denoting number of iterations required to reach convergence)
% Xa_stack = retention of Xa at each iteration 
% Xa_noise = the noise components of the final Xa matrix; calculated
%   using the (total(EOF) - max(nEOF)) components
% r2_cv = cross-validation squared pearson correlation coefficient

if nargin < 3; verbal = 'false';   end
verbal_true = strcmp(verbal,'true');
verbal_false = strcmp(verbal,'false');
if ischar(verbal) && (verbal_true ~= 1) && (verbal_false ~= 1)
    warning('Please input either ''true'' or ''false'' for the ''verbal'' input and try again.'); return;
elseif ~ischar(verbal) 
    warning('Please input either ''true'' or ''false'' for the ''verbal'' input and try again.'); return; end
    
%% Initialization

Im = isnan(Xo); % locate the location of the missing values in the data field, Xo
    if sum(Im(:)) == 0 % quit the code right here if there are no missing values in Xo
        if verbal_true == 1
        disp(' '); disp('Caution! There are no missing values (NaNs) in the dataset!')
        disp('Infilling all zero-values with NaNs, but double checking the input Xo is recommended'); disp(' ')
        end
        pause(1)
        Xo(Xo == 0) = NaN;
        Im = isnan(Xo);
    end
    
    Ir = find(~isnan(Xo)); % takes all the real number values in Xo, and stacks them as a vector
    % now randomly reorder the real elements in "all_real_values":
    ix = randsample(Ir,length(Ir));% ix contains the positions of the real numbers in Xo in random order

    % use the randomly ordered real_values of Xo to fill into a vector, initialize_missing_values_for_Xo, 
    % which will be used as a worthy a priori "guess" for the root-mean-square error field
    if floor(0.05*length(Ir(:))) > 30 % if 5% of the number of real values in Xo is greater than 30
        ref_pos = zeros(floor(0.05*length(Ir(:))),1); % count the number of missing values
    else % if 30 is greater than 5% of the real values in Xo, set as 30
        ref_pos = zeros(30,1);
    end
        % sample the first i indices of the real values created in ix
        for i = 1:length(ref_pos)
            ref_pos(i) = ix(i);
        end
        
    % initialize and define Xa as Xo with 0 input into the missing values(Pg. 184) 
    Xa = Xo; 
        Xa(Im) = 0; % input zeros into the Xa matrix, from the locations of missing data in Xo
        Xa(ref_pos(:)) = 0; % input zeros into elements of Xo that were originally numbers
        
        rmse_prev = inf;
        rmse_curr = sqrt(mean((Xa(ref_pos) - Xo(ref_pos)).^2)); 

    % initialize some variables
    neof_curr = 1; % k = mode of variability included; initialize as 1 (the first mode of variability);
    rmse = rmse_curr; % keeps track of the RMS at each iteration
    neof = neof_curr; % keeps track of the EOF at each iteration
    Xa_stack = Xa; % keeps track of the updated Xa at each iteration
    Xa_best = Xa; % records the optimal RMS
    neof_best = neof; % record the optimal number of EOFs
    r2_cv = 0;
    
    while (rmse_prev - rmse_curr) > delrms   &&    neof_curr < length(Xo(1,:)) % loop for increasing number of eofs
        
        while (rmse_prev - rmse_curr) > delrms  % loop for increasing number of eofs
            rmse_prev = rmse_curr; % update the rms value for this iteration to the prior iteration
            [U, S, V] = svd(Xa,0); % SVD decomposition - Xa is an updated version of itself at the missing + randomly selected data points each time
            rec_X = U(:,1:neof_curr)*S(1:neof_curr,1:neof_curr)*V(:,1:neof_curr)'; % recontruct the X data field using n_eof fields
            Xa(Im) = rec_X(Im); % input revised values for locations originally missing data back into the Xa matrix
            Xa(ref_pos(:)) = rec_X(ref_pos(:)); % input the recovered EOFs values at the locations of the randomly removed values back into Xa 
            rmse_curr = sqrt(mean((Xa(ref_pos) - Xo(ref_pos)).^2)); 
            
            % R2 = corrcoef(Xa(ref_pos),Xo(ref_pos)); R2 = R2(2,1)^2;
            r2 = 1 - (mean((Xa(ref_pos) - Xo(ref_pos)).^2))./mean(Xo(ref_pos).^2);  % from Michaelsen et al., 1987
               
            if verbal_true == 1
            disp(['EOF-',num2str(neof_curr),'; RMSE = ',num2str(rmse_curr),'; R2_cv = ',num2str(r2),])
            end
            
            rmse =      vertcat(rmse,rmse_curr); % track the RMSE!
            neof =    vertcat(neof,neof_curr);
            Xa_stack =  cat(3,Xa_stack,Xa);
            r2_cv = vertcat(r2_cv,r2);
            if rmse_curr == min(rmse)
            	Xa_best = Xa;
            	neof_best = neof_curr;
            end
        end
       
        neof_curr = neof_curr + 1;
        rmse_prev = rmse_curr;
       
        [U, S, V] = svd(Xa,0);
        rec_X = U(:,1:neof_curr)*S(1:neof_curr,1:neof_curr)*V(:,1:neof_curr)';
       
        Xa(Im) = rec_X(Im); % input revised values from rec_X into the Xa matrix
        Xa(ref_pos(:)) = rec_X(ref_pos(:)); % input revised values from rec_X
        
        rmse_curr = sqrt(mean((Xa(ref_pos) - Xo(ref_pos)).^2)); 
       
        % R2 = corrcoef(Xa(ref_pos),Xo(ref_pos)); R2 = R2(2,1)^2;
        r2 = 1 - (mean((Xa(ref_pos) - Xo(ref_pos)).^2))./mean(Xo(ref_pos).^2);  % from Michaelsen et al., 1987
               
        if verbal_true == 1
            disp(['EOF-',num2str(neof_curr),'; RMSE = ',num2str(rmse_curr),'; R2_cv = ',num2str(r2),])
        end
       
        rmse         = vertcat(rmse,rmse_curr); % track the RMSE!
        neof       = vertcat(neof,neof_curr);
        Xa_stack     = cat(3,Xa_stack,Xa);
        r2_cv = vertcat(r2_cv,r2);
        if rmse_curr == min(rmse)
            Xa_best = Xa;
            neof_best = neof_curr;
        end
                
    end
    
    Xa          = Xa_best;
    neof_curr   = neof_best;
    Xa(ref_pos) = Xo(ref_pos); % add the withheld points back in to the infilled matrix
   
    % option to add noise back in
    [U, ~, ~] = svd(Xa,0);
    remainder = neof_best+1 : size(U,2);
        % Xa_signal = U(:,1:nEOF_best)*U(:,1:nEOF_best)'*Xa;
        Xa_noise = U(:,remainder)*U(:,remainder)'*Xa;
       
    if verbal_true == 1; disp(' '); end
    r2_cv(isnan(r2_cv)) = [];
   
end
