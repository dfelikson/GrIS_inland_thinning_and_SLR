function sout = interpWVdem(glacier,X,Y),

% Translate glacier name
glacier_translated = glacier_abbrs(glacier, 'new_abbr', 'centerline_name');

DEMdir='/Users/denisfelikson/Research/Projects/CentralWestGrISGlaciers/Data/DEMs/WVfilled-WVmerged';

% Find the DEM corresponding to the glacier
d = dir([DEMdir filesep glacier_translated '*-20m.tif']);
if numel(d) == 0
   disp(['WARNING: No WorldView DEMs found for glacier ' glacier '.'])
   sout = [];
   return
elseif numel(d) > 1
   disp(['WARNING: Multiple WorldView DEMs found for glacier ' glacier '. Must resolve manually.'])
   sout = [];
   return
else
   DEMpath = [DEMdir filesep d(1).name];
end

usemap = 1;
if license('test','map_toolbox')==0,
	disp('WARNING: map toolbox not installed, trying house code');
	usemap = 0;
elseif license('checkout','map_toolbox')==0
	disp('WARNING: map toolbox not available (checkout failed), trying house code');
	usemap = 0;
end

if usemap,
	[data,R] = geotiffread(DEMpath);
	data=double(flipud(data));
	xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
	xdata =(xdata(1:end-1)+xdata(2:end))/2;
	ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
	ydata =(ydata(1:end-1)+ydata(2:end))/2;
else

	%Get image info
	Tinfo = imfinfo(DEMpath);
	N     = Tinfo.Width;
	M     = Tinfo.Height;
	dx    = Tinfo.ModelPixelScaleTag(1);
	dy    = Tinfo.ModelPixelScaleTag(2);
	minx  = Tinfo.ModelTiepointTag(4);
	maxy  = Tinfo.ModelTiepointTag(5);

	%Generate vectors
	xdata = minx + dx/2 + ((0:N-1).*dx);
	ydata = maxy - dy/2 - ((M  -1:-1:0).*dy);

	%Read image
	data=double(flipud(imread(DEMpath)));
end

sout = InterpFromGrid(xdata,ydata,data,X,Y);
