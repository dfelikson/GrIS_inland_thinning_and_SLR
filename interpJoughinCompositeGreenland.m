function [vxout vyout] = interpLandsatVelocityGreenland(X,Y),

vxout = [];
vyout = [];

vx_raster = '/Users/dfelikso/Research/Data/Velocity/MEaSUREs/greenland_vel_mosaic250_vx_v1.tif';
vy_raster = '/Users/dfelikso/Research/Data/Velocity/MEaSUREs/greenland_vel_mosaic250_vy_v1.tif';

[vx,R] = geotiffread(vx_raster);
[vy,~] = geotiffread(vy_raster);

vx(vx==-2000000000) = nan;
vy(vy==-2000000000) = nan;

vx=double(flipud(vx)); vy=double(flipud(vy));

xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
xdata =(xdata(1:end-1)+xdata(2:end))/2;
ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
ydata =(ydata(1:end-1)+ydata(2:end))/2;

vxout = InterpFromGrid(xdata,ydata,vx,X,Y);
vyout = InterpFromGrid(xdata,ydata,vy,X,Y);

