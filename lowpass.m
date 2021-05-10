function [smoothed,mse] = lowpass(indata,frequency,iconb,icone)
% Adapted by M. Osman from Kinnard et al.(2011), Nature
% Aug. 2017;
%
% [smoothed,mse] = lowpass(indata,frequency,iconb,icone)
% lowpass filter.
% assumes unit time spaced series 
% 
% lowpass cutoff at f="frequency" (in cycles/time unit)
%
% icon is any one of several possible boundary constraints on the smooth at 
% the begin/end x limits
% (0) minimum norm,
% (1) minimum slope
% (2) minimum roughness 
%
% see:
% Park, J., Envelope estimation for quasi-periodic geophysical 
%  signals in noise: A multitaper approach, in Statistics in the Environmental 
%  and Earth Sciences, edited by A. T. Walden and P. Guttorp, pp. 189-219, 
%  Edward Arnold, London, 1992; 
% Ghil, M., Allen, M.R., Dettinger, M.D., Ide, K., Kondrashov, D., Mann, M.E., 
%  Robertson, A.W., Tian, Y., Varadi, F., Yiou, P., 
%  Advanced Spectral Methods for Climatic Time Series, Reviews of Geophysics, 
%  40 (1), 1003, doi: 10.1029/2000RG000092, 2002.]
%
% The lowpass routine employs a 10 point butterworth filter. 

% Rather than implementing contraints (0)-(2) in the frequency domain (as in Park, 
% Ghil et al) we use the following approximate implementations of the boundary 
% constraints
% 
% (0) pad series with mean value beyond x boundary
% (1) pad series with values over last 1/2 filter width reflected w.r.t. x
%     [imposes a local maximum/minimum beyond x boundary]
% (2) pad series with values over last 1/2 filter width reflected w.r.t. x and y (w.r.t. final value)
%     [imposes a point of inflection beyond x boundary]
%
ipts=10;
fn=frequency*2;
npad=1/fn;
nn=length(indata);
padded(npad+1:npad+nn)=indata;
padded(npad+nn+1:nn+2*npad)=indata(nn);
padded(1:npad)=indata(1);
for j=nn+npad+1:nn+2*npad
   ipad=j-nn-npad;
   if (icone==0) 
       apad=mean(indata);
   else if (icone==1)
           apad=indata(nn-ipad);
        else
           apad=2*indata(nn)-indata(nn-ipad);
        end
   end
   padded(j:j)=apad;
end
for j=1:npad
   ipad=j;
   if (iconb==0) 
       apad=mean(indata);
   else if (iconb==1)
           apad=indata(npad-ipad+1);
        else
           apad=2*indata(1)-indata(npad-ipad+1);
        end
   end
   padded(j:j)=apad;
end
%
% smoothed=padded;
[b,a]=butter(ipts,fn,'low');
smoothed0=filtfilt(b,a,padded);
smoothed=smoothed0(npad+1:nn+npad)';
%
% determine mse relative to raw data
%
asum = 0.0;
resid=smoothed-indata;
mse=var(resid)/var(indata);





