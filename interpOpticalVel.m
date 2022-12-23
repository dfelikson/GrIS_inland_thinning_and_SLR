function [vxout vyout] = interpOpticalVel(X,Y,Location,Date),
%  [vxout vyout] = interpOpticalVel(X,Y,Location,Date)
%     Location must be a string with format: (e.g.) Wcoast-69.95N
%     Date must be a string with format: yyyy-mm

% EQI: Wcoast-69.95N

switch oshostname(),
	case {'nilas.local'}
      rootname = ['/Users/denisfelikson/GoogleDrive/Research/Data/Velocity/MEaSUREs/nsidc-0646'];
	otherwise
		error('machine not supported yet');
end
verbose = 1;

if ~exist(rootname,'dir'),
	error(['file ' rootname ' not found']);
end

rootname = [rootname '/' Location '/' 'OPT_' strrep(Location,'coast-','') '_' Date '/' 'OPT_' strrep(Location,'coast-','') '_' Date];

if verbose, disp('   -- Optical: loading vx'); end
[data,R] = geotiffread([rootname '.vx.tif']);
pos=find(data==-99999.000); data(pos)=NaN;
data=double(flipud(data));
xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
xdata =(xdata(1:end-1)+xdata(2:end))/2;
ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
ydata =(ydata(1:end-1)+ydata(2:end))/2;
if verbose, disp('   -- Optical: interpolating vx'); end
vxout = InterpFromGrid(xdata,ydata,data,X,Y); %,'nearest');
vxout = reshape(vxout,size(X,1),size(X,2));

if verbose, disp('   -- Optical: loading vy'); end
[data,R] = geotiffread([rootname '.vy.tif']);
pos=find(data==-99999.000); data(pos)=NaN;
data=double(flipud(data));
xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
xdata =(xdata(1:end-1)+xdata(2:end))/2;
ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
ydata =(ydata(1:end-1)+ydata(2:end))/2;
if verbose, disp('   -- Optical: interpolating vy'); end
vyout = InterpFromGrid(xdata,ydata,data,X,Y); %,'nearest');
vyout = reshape(vyout,size(X,1),size(X,2));

%v = sqrt(vxout.^2 + vyout.^2);
%pos = find(v<=0);
%vxout(pos) = NaN; vyout(pos) = NaN;
return

% There's more in interpJoughin.m ...

