function [vxout vyout] = interpLandsatVelocityGreenland(X,Y),

vxout = [];
vyout = [];

glacier = evalin('base','glacier');

directory = '~/Research/Data/Velocity/Landsat';
velocity_raster = '';
if exist(['Vel/' evalin('base','glacier') '_vel.tif'],'file')
   velocity_raster = ['Vel/' evalin('base','glacier') '_vel.tif'];
else
   velocity_raster = [directory '/GRL_merged_month_Aug_1985-1999_filtered_vel.tif'];
   fprintf(['\n\033[' '103;30' 'm   WARNING: ' 'Vel/' evalin('base','glacier') '.tif not found; defaulting to: ' velocity_raster '\033[0m \n\n']);
end

velocity_raster = strrep(velocity_raster, '_vel.tif', '');

if ~exist(strcat(velocity_raster,'_vx.tif'),'file') | ~exist(strcat(velocity_raster,'_vy.tif'),'file')
   fprintf(['\n\033[' '103;30' 'm   ERROR: Could not find vx and vy tifs for ' velocity_raster '\033[0m \n\n']);
   return
end

[vx,R] = geotiffread(strcat(velocity_raster,'_vx.tif'));
[vy,~] = geotiffread(strcat(velocity_raster,'_vy.tif'));

vx(vx==-9999) = nan;
vy(vy==-9999) = nan;

vx=double(flipud(vx)); vy=double(flipud(vy));

xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
xdata =(xdata(1:end-1)+xdata(2:end))/2;
ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
ydata =(ydata(1:end-1)+ydata(2:end))/2;

vxout = InterpFromGrid(xdata,ydata,vx,X,Y) * 365.25;
vyout = InterpFromGrid(xdata,ydata,vy,X,Y) * 365.25;

vxout(vxout==-3.652134750000000e+06) = nan;
vyout(vyout==-3.652134750000000e+06) = nan;
