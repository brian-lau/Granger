% Format data for fitting autoregressive model using ARFIT or ARFIT2
%
%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     Released under the BSD license. The license and most recent version
%     of the code can be found on GitHub:
%     https://github.com/brian-lau/Granger
%
%   datChunk = format_arfit_input(dat.data,winStart(i),winEnd(i),flag,demean);

function [datChunk,nContr] = format_arfit_input(dat,winStart,winEnd,flag,demean)

ind = winStart:winEnd;

% For ARFIT2
minTrialDur = ceil(0.10*length(ind)); % each trial must be at least this long
minTrials = 3; % This many trials must contribute to the average for demeaning
nTrials = length(dat);

switch flag
   case 'arfit'
      for j = 1:nTrials
         datChunk(:,:,j) = dat(j).v(ind,:);
      end
      if demean
         temp = mean(datChunk,3);
         datChunk = datChunk - repmat(temp,[1 1 nTrials]);
      end
   case 'arfit2'
      % Put each trial into a 3-d array, filled with NaNs where there is no data
      n = size(dat(1).v,2);
      datChunk = NaN*zeros(length(ind),n,nTrials);
      for j = 1:nTrials
         % Trial length
         n2 = length(dat(j).v);
         % If trial ends before start of requested window
         if n2 > (winStart+minTrialDur)
            ind2 = winStart:min(winEnd,n2);
            datChunk(1:length(ind2),:,j) = dat(j).v(ind2,:);
         end
      end

      if demean
         % count # of contributing trials
         nContr = sum(~isnan(datChunk),3);
         % mean across trials
         temp = nanmean(datChunk,3);
         % NaN-out those points that don't have enough trials contribtuting to an average
         temp2 = temp;
         temp2(nContr<minTrials) = NaN;
         datChunk = datChunk - repmat(temp2,[1 1 nTrials]);
      end
      
      if sum(isnan(temp2)) ==  length(ind)
         datChunk = [];
      else
         % Format for ARFIT2
         v = [];
         for j = 1:nTrials
            ind = find(isnan(datChunk(:,1,j)));
            if isempty(ind)
               % data fills window
               v = [v ; NaN*ones(1,n) ; datChunk(:,:,j)];
            else
               % data shorter than window
               ind = ind(1);
               v = [v ; NaN*ones(1,n) ; datChunk(1:ind,:,j)];
            end
         end
         datChunk = v;
      end
   otherwise
      error('Bad FLAG to FORMAT_ARFIT_INPUT');
end

