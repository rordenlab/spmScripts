function filt_data = fsl_temporal_filt(orig_data,HP_sigma,LP_sigma,flag)

% Perform bandpass filtering à la FSL (non-linear highpass and Gaussian 
% linear lowpass). Filters can either be entered in terms of the cutoff (in
% volumes not seconds) or the filter itself. Set either sigma<0 to skip
% that filter. 
% 
% There are several ways to calculate sigma from traditional Hz cutoff
% value (2nd method used in GTG):
%      Typical FSL method:  sigma = ((1./Hz_cutoff)./2)./TR;
%      Slightly more exact: sigma = ((1./Hz_cutoff)./sqrt(8*log(2)))./TR;
%      Alternative method:  sigma = (1/TR)./(2*pi*Hz_cutoff);
% 
% 
% 
% Usage:  filt_data = fsl_temporal_filt(orig_data,HP_sigma,LP_sigma,flag);
% 
% Inputs: orig_data = input data
%         HP_sigma  = high-pass cutoff sigma (in volumes) or filter
%         LP_sigma  = low-pass cutoff sigma (in volumes) or filter
%         flag      = method of dealing with values outside bounds; 1 
%                     (default): filter is truncated (FSL method), 2: 
%                     values outside bounds are computed via reflection
% 
% Output: filt_data = filtered data
% 
%
%
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 02.09.16
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

orig_data = double(squeeze(orig_data(:)));

if ~exist('HP_sigma','var')
    HP_sigma = -1;
elseif ~exist('LP_sigma','var')
    LP_sigma = -1;
end

if ~exist('flag','var')
    flag = 1;
end

if length(HP_sigma)==1
    if HP_sigma>0
        HP_filt_size = round(HP_sigma*8);
        HP_lin = linspace((-HP_filt_size/2),(HP_filt_size/2),HP_filt_size);
        HP_gfilt = exp(-(HP_lin.^2)/(2*(HP_sigma^2)));
        HP_gfilt = HP_gfilt/sum(HP_gfilt);
    end
    if LP_sigma>0
        LP_filt_size = round(LP_sigma*8);
        LP_lin = linspace((-LP_filt_size/2),(LP_filt_size/2),LP_filt_size);
        LP_gfilt = exp(-(LP_lin.^2)/(2*(LP_sigma^2)));
        LP_gfilt = LP_gfilt/sum(LP_gfilt);
    end
else
    HP_gfilt = HP_sigma(:)';
    LP_gfilt = LP_sigma(:)';
end

if HP_sigma>0
    if flag==1
        filt_data = zeros(size(orig_data));
        back  = floor((HP_filt_size-1)/2);
        front = ceil((HP_filt_size-1)/2);
        for t = 1:length(orig_data)
            if (t-back<1) && (t+front>length(orig_data))
                trunc_HP_gfilt = HP_gfilt((back-t+2):(end-(t+front-length(orig_data))));
                trunc_HP_gfilt = trunc_HP_gfilt/sum(trunc_HP_gfilt);
                filt_data(t) = trunc_HP_gfilt*orig_data;
            elseif t-back<1
                trunc_HP_gfilt = HP_gfilt((back-t+2):end);
                trunc_HP_gfilt = trunc_HP_gfilt/sum(trunc_HP_gfilt);
                filt_data(t) = trunc_HP_gfilt*orig_data(1:t+front);
            elseif t+front>length(orig_data)
                trunc_HP_gfilt = HP_gfilt(1:(end-(t+front-length(orig_data))));
                trunc_HP_gfilt = trunc_HP_gfilt/sum(trunc_HP_gfilt);
                filt_data(t) = trunc_HP_gfilt*orig_data(t-back:end);
            else
                filt_data(t) = HP_gfilt*orig_data(t-back:t+front);
            end
        end
        filt_data = orig_data-filt_data;
    elseif flag==2
        filt_data = orig_data-imfilter(orig_data,HP_gfilt','symmetric');
    end
else
    filt_data = orig_data;
end

if LP_sigma>0
    if flag==1
        filt_data_orig = filt_data;
        back  = floor((LP_filt_size-1)/2);
        front = ceil((LP_filt_size-1)/2);
        for t = 1:length(filt_data)
            if (t-back<1) && (t+front>length(filt_data_orig))
                trunc_LP_gfilt = LP_gfilt((back-t+2):(end-(t+front-length(filt_data_orig))));
                trunc_LP_gfilt = trunc_LP_gfilt/sum(trunc_LP_gfilt);
                filt_data(t) = trunc_LP_gfilt*filt_data_orig(t-back:t+front);
            elseif t-back<1
                trunc_LP_gfilt = LP_gfilt((back-t+2):end);
                trunc_LP_gfilt = trunc_LP_gfilt/sum(trunc_LP_gfilt);
                filt_data(t) = trunc_LP_gfilt*filt_data_orig(1:t+front);
            elseif t+front>length(filt_data_orig)
                trunc_LP_gfilt = LP_gfilt(1:(end-(t+front-length(orig_data))));
                trunc_LP_gfilt = trunc_LP_gfilt/sum(trunc_LP_gfilt);
                filt_data(t) = trunc_LP_gfilt*filt_data_orig(t-back:end);
            else
                filt_data(t) = LP_gfilt*filt_data_orig(t-back:t+front);
            end
        end
    elseif flag==2
        filt_data = imfilter(filt_data,LP_gfilt','symmetric');
    end
end