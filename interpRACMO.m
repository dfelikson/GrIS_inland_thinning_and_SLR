function [dataout] = interpRACMO(X,Y,string,options_values)
%interpRACMO - interpolate RACMO data
%
%   Available data:
%      1.  smb
%      2.  smb_downscaled
%      3.  runoff
%
%   Usage:
%      [dataout] = interpRACMO(X,Y,string)
%      [dataout] = interpRACMO(X,Y,string,options_values)

if ~exist('options_values','var')
   options_values = {};
end

dataout = [];

% Glacier
idx = find(strcmp(options_values,'glacier'));
if ~isempty(idx)
   glacier = options_values{idx+1};
else
   glacier = evalin('base','glacier');
end

%Read data{{{
switch string
   case 'smb_downscaled'
      
      if ~isempty(find(strcmp(options_values,'stable')))
         % This is the mean SMB from 01Jan1971-31Dec1988 (so, a "stable ice sheet" mean climate):
         racmodatasets = {['/Users/dfelikso/Research/Data/RACMO/RACMO2.3/RACMO2.3p1anomaly/' glacier '_smb_1971-1988_mean.tif']};
         
      elseif ~isempty(find(strcmp(options_values,'present-day')))
         % This is the mean SMB from 01Aug2000-31Jul2015 (so, the present-day climate):
         racmodatasets={'/Users/dfelikso/Research/Data/RACMO/RACMO2.3/RACMO2.ep1anomaly/GrIS_smb_downscaled_mean_01Aug2000-31Jul2015_fromMean_01Aug2000-31Jul2015_dh.tif'};

      else
         idx = find(strcmp(options_values,'startyear'));
         if ~isempty(idx)
            % Parse start/end years
            idx = find(strcmp(options_values,'startyear'));
            startyear = options_values{idx+1};
            idx = find(strcmp(options_values,'endyear'));
            endyear = options_values{idx+1};
            
            years = startyear:endyear;
         else
            % Find requested years
            idx = find(strcmp(options_values,'years'));
            years = options_values{idx+1};
         end

         disp(['     -- RACMO: loading time-varying SMB from ' sprintf('%4d',years(1)) ' to ' sprintf('%4d',years(end))]);
         % Look for files on local machine
         if exist(['/Users/dfelikso/Research/Data/RACMO/RACMO2.3/RACMO2.3p1anomaly/' glacier '_smb_' sprintf('%4d',years(1)) '.tif'],'file') & ...
            exist(['/Users/dfelikso/Research/Data/RACMO/RACMO2.3/RACMO2.3p1anomaly/' glacier '_smb_' sprintf('%4d',years(end)) '.tif'],'file')

            racmodatasets = {};
            for year = years
               racmodatasets{end+1} = ['/Users/dfelikso/Research/Data/RACMO/RACMO2.3/RACMO2.3p1anomaly/' glacier '_smb_' sprintf('%4d',year) '.tif'];
            end

         else
            % Bounding box
            xmin = strtrim(sprintf('%16.6f', min(X)));
            xmax = strtrim(sprintf('%16.6f', max(X)));
            ymin = strtrim(sprintf('%16.6f', min(Y)));
            ymax = strtrim(sprintf('%16.6f', max(Y)));
            boundingbox = sprintf('''%s %s %s %s''',xmin,xmax,ymin,ymax);

            fprintf(['\n\033[' '103;30' 'm   WARNING: RACMO SMB tiffs not set up for ' glacier '.' '\033[0m']);
            fprintf(['\n\033[' '103;30' 'm            Use /home/student/denis/ModeledInlandThinning/Analysis/climate/clipRACMO_avgYear.py on melt.' '\033[0m']);
            fprintf(['\n\033[' '103;30' 'm            Bounding box: ' boundingbox '\033[0m \n']);
            return

            % % Call a script on melt to process
            % cmdstr = ['ssh -t denis@melt.ig.utexas.edu /home/student/denis/ModeledInlandThinning/Analysis/climate/clip_RACMO_smb.sh ' ...
            % glacier ' ' startdate ' ' enddate ' ' boundingbox];
            % system(cmdstr);
            % return
         end
      end

   otherwise
      error(['variable ' string ' not supported'])
end
%}}}

dataout = nan * ones(length(X),numel(racmodatasets));
for i = 1:numel(racmodatasets)
   racmodataset = racmodatasets{i};
   [filepath, filename, fileext] = fileparts(racmodataset);

   %if strcmp(fileext,'.nc')
   %%{{{
   %   if strcmp(opt,'years')
   %      startyear = optvalue1;
   %      endyear   = optvalue2;
   %      time = [];
   %      
   %      disp('     -- RACMO: loading coordinates');
   %      racmofile = strrep(racmodataset,'????',num2str(startyear));
   %      if ~exist(racmofile)
   %         disp(['missing file: ' racmofile])
   %         return
   %      end
   %      
   %      xdata = double(ncread(racmofile,'x'));
   %      ydata = double(ncread(racmofile,'y'));
   %      dataout = [];
   %
   %      for year = startyear:endyear
   %         disp(['     -- RACMO: reading year: ' num2str(year)]);
   %         time = [time; ncread(racmofile,'time')];
   %         for month = 1:12
   %            datayear = double(ncread(racmofile, string, [1 1 month], [length(xdata) length(ydata) 1], [1 1 1]))';
   %            dataout = cat(2,dataout,InterpFromGrid(xdata,ydata,datayear,X,Y));
   %         end
   %      end
   %   else
   %      disp('must specify date range')
   %      return
   %   end
   %end
   %%}}}

   if strcmp(fileext,'.tif')
      disp(['     -- RACMO: loading ' string ' from ' filename]);
      [data,R] = geotiffread(racmodataset);
      pos=find(data==0); data(pos)=NaN;
      data=double(flipud(data));
      xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
      xdata =(xdata(1:end-1)+xdata(2:end))/2;
      ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
      ydata =(ydata(1:end-1)+ydata(2:end))/2;

      disp(['     -- RACMO: interpolating ' string]);
      if strcmpi(string,'LandMask');
         smb = InterpFromGrid(xdata,ydata,data,X,Y,'nearest');
      else
         smb = InterpFromGrid(xdata,ydata,data,X,Y);
      end
      smb(isnan(smb)) = 0;
      smb(smb==-9999) = 0;
      dataout(:,i) = reshape(smb,size(X,1),size(X,2));
   end
end

