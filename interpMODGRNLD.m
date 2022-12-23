function output = interpMODGRNLD(X,Y),

disp('   -- MODGRNLD: interpolating surface temperature');
ncfile = '/Users/dfelikso/Research/Data/MODGRNLD/MODGRNLD.2001.yearly.v01.1.nc';
disp(['   -- MODGRNLD: ' ncfile]);
xdata = double(ncread(ncfile,'x'));
ydata = double(ncread(ncfile,'y'));
data  = double(ncread(ncfile,'Ice_Surface_Temperature_Mean')');

output = InterpFromGrid(xdata,ydata,data,double(X),double(Y));
%output = InterpFromGrid(xdata,ydata,data,double(X),double(Y),'nearest');

