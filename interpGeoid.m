function output = interpGeoid(X,Y),

prefix = 'BedMachineGreenland';
ncdate='2017-09-20';
bednc = sprintf('~/Research/Data/GreenlandBed/MCbed/%s-%s/%s-%s.nc', prefix, ncdate, prefix, ncdate);
%switch oshostname(),
%	case {'murdo','thwaites','astrid'}
%		bednc=['/u/astrid-r1b/morlighe/issmjpl/proj-morlighem/DatasetGreenland/Output/MCdataset-' ncdate '.nc']';
%	case {'ronne'}
%		bednc=['/home/ModelData/Greenland/BedMachine/MCdataset-' ncdate '.nc'];
%	otherwise
%		error('machine not supported yet');
%end

disp(['   -- Bedmachine Greenland: loading geoid']);
xdata = double(ncread(bednc,'x'));
ydata = double(ncread(bednc,'y'));

offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(xdata<=xmax);
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
posy=find(ydata>=ymin);
id1y=max(1,find(ydata<=ymax,1)-offset);
id2y=min(numel(ydata),posy(end)+offset);

data  = double(ncread(bednc,'geoid',[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
xdata=xdata(id1x:id2x);
ydata=ydata(id1y:id2y);
data(find(data==-9999))=NaN;

output = InterpFromGrid(xdata,ydata,data,double(X),double(Y));
