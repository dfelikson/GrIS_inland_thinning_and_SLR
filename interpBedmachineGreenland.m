function output = interpBedmachineGreenland(X,Y,string,ncdate),

if nargin<3, string = 'bed'; end
if nargin<4,
	%ncdate='2013-05-21';
	%ncdate='2013-06-27';
	%ncdate='2013-07-18';
	%ncdate='2013-11-15';
	%ncdate='2013-12-03';
	%ncdate='2014-02-26';
	%ncdate='2014-03-24';
	%ncdate='2014-07-31';
	%ncdate='2014-11-14';
	%ncdate='2015-03-03';
	%ncdate='2015-03-10';
	%ncdate='2015-03-26';
	%ncdate='2015-03-30';
	%ncdate='2015-04-27'; %BedMachine v2
	%ncdate='2015-07-30';
	%ncdate='2015-10-02';
	%ncdate='2016-03-21';
	%ncdate='2016-05-12';
	ncdate='2016-07-06';
	ncdate='2016-08-04';
	ncdate='2016-10-26';
	ncdate='2016-11-23';
	ncdate='2016-12-21';
	ncdate='2017-01-19';
	ncdate='2017-03-28';
	ncdate='2017-05-10';
	ncdate='2017-07-21';
	ncdate='2017-09-25'; %BedMachine v3
	ncdate='2017-09-20'; %BedMachine v3ish
   ncdate='v5';
end

%if exist('datetime','file') 
%	date1 = sscanf(ncdate,'%d-%d-%d');
%	date2 = datetime(date1(1),date1(2),date1(3));
%	if date2<datetime(2016,10,24),
%		basename = 'MCdataset'; 
%	else
%		basename = 'BedMachineGreenland';
%	end
%else
%  basename = 'BedMachineGreenland';
%end

switch strtok(oshostname(),'.'),
   case {'gs615serac','gsslw17081229156068103180'}
      morlighem2013nc=['/Users/dfelikso/Research/Data/GreenlandBed/MCbed/BedMachineGreenland-' ncdate '/BedMachineGreenland-' ncdate '.nc'];
	case {'murdo','thwaites','astrid'}
		morlighem2013nc=['/u/astrid-r1b/ModelData/ModelData/MCdataset-' ncdate '.nc']';
	case {'ronne'}
		morlighem2013nc=['/home/ModelData/Greenland/BedMachine/' basename '-' ncdate '.nc'];
	otherwise
      morlighem2013nc=['/Users/dfelikso/Research/Data/GreenlandBed/MCbed/BedMachineGreenland-' ncdate '/BedMachineGreenland-' ncdate '.nc'];
      if ~exist(morlighem2013nc, 'file')
         error('machine not supported yet');
      end
end

disp(['   -- BedMachine Greenland version: ' ncdate]);
xdata = double(ncread(morlighem2013nc,'x'));
ydata = double(ncread(morlighem2013nc,'y'));

offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(xdata<=xmax);
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
posy=find(ydata>=ymin);
id1y=max(1,find(ydata<=ymax,1)-offset);
id2y=min(numel(ydata),posy(end)+offset);

disp(['   -- BedMachine Greenland: loading ' string]);
data  = double(ncread(morlighem2013nc,string,[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
xdata=xdata(id1x:id2x);
ydata=ydata(id1y:id2y);
data(find(data==-9999))=NaN;

disp(['   -- BedMachine Greenland: interpolating ' string]);
if strcmp(string,'mask') | strcmp(string,'source'),
	%Need nearest neighbor to avoid interpolation between 0 and 2
	output = InterpFromGrid(xdata,ydata,data,double(X),double(Y),'nearest');
else
	output = InterpFromGrid(xdata,ydata,data,double(X),double(Y));
end

%TEST https://www.mathworks.com/matlabcentral/fileexchange/10772-fast-2-dimensional-interpolation
