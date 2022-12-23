function sout = interpArcticDEM(glacier,X,Y),

dem=['./DEMs/' glacier '_ArcticDEM.tif'];

usemap = 1;
if license('test','map_toolbox')==0,
	disp('WARNING: map toolbox not installed, trying house code');
	usemap = 0;
elseif license('checkout','map_toolbox')==0
	disp('WARNING: map toolbox not available (checkout failed), trying house code');
	usemap = 0;
end

if usemap,
	[data,R] = geotiffread(dem);
	data=double(flipud(data));
	xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
	xdata =(xdata(1:end-1)+xdata(2:end))/2;
	ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
	ydata =(ydata(1:end-1)+ydata(2:end))/2;
else

	%Get image info
	Tinfo = imfinfo(dem);
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

surf = InterpFromGrid(xdata,ydata,data,X,Y);
surf(surf == -9999) = nan;

% Must reference ArcticDEM to geoid
geoid = interpGeoid(X,Y);
geoid(geoid == -9999) = nan;
sout = surf - geoid;

