function veli = interpRosenauLandsatVel(X,Y,opt,optvalue1,optvalue2,optvalue3),

veli = [];

velnc='/Users/denisfelikson/GoogleDrive/Research/Data/Velocity/Landsat/Rosenau/GRL_002_all.EPSG3413.vel_md.nc';
%switch oshostname(),
%	case {'murdo','thwaites','astrid'}
%		morlighem2013nc=['/u/astrid-r1b/morlighe/issmjpl/proj-morlighem/DatasetGreenland/Output/MCdataset-' ncdate '.nc']';
%	case {'ronne'}
%		morlighem2013nc=['/home/ModelData/Greenland/BedMachine/MCdataset-' ncdate '.nc'];
%	otherwise
%		error('machine not supported yet');
%end

dn0  = datenum('1970-01-01 00:00:00');
time = ncread(velnc,'time');
dn   = dn0 + time;
switch opt
   case 'bestcoverage'
      % This will select the valocity field that populates the most grid nodes between two dates
      startdn = datenum(optvalue1);
      enddn   = datenum(optvalue2);
      id1t = find(dn >= startdn, 1, 'first');
      id2t = find(dn <= enddn,   1, 'last');

   case 'bestcoveragemonth'
      startdn = datenum(optvalue1);
      enddn   = datenum(optvalue2);
      month = optvalue3;
      id1t = find(dn >= startdn, 1, 'first');
      id2t = find(dn <= enddn,   1, 'last');

   case 'date'
      selectdn = datenum(optvalue1);
      [~, id1t] = min(abs(dn-selectdn));
      id2t = id1t;
end

disp(['   -- Rosenau Landsat velocity: loading ' datestr(dn(id1t)) ' to ' datestr(dn(id2t))]);
xdata = double(ncread(velnc,'x'));
ydata = double(ncread(velnc,'y'));

offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(xdata<=xmax);
if isempty(posx)
   disp('ERROR: Rosenau velocity doesn''t span glacier domain.')
   return
end
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
posy=find(ydata>=ymin);
if isempty(posy)
   disp('ERROR: Rosenau velocity doesn''t span glacier domain.')
   return
end
id1y=max(1,find(ydata<=ymax,1)-offset);
id2y=min(numel(ydata),posy(end)+offset);

xdata=xdata(id1x:id2x);
ydata=ydata(id1y:id2y);

vx  = double(ncread(velnc,'vx',[id1x id1y id1t],[id2x-id1x+1 id2y-id1y+1 id2t-id1t+1],[1 1 1]));
vx(find(vx==-9999))=NaN;
vy  = double(ncread(velnc,'vy',[id1x id1y id1t],[id2x-id1x+1 id2y-id1y+1 id2t-id1t+1],[1 1 1]));
vy(find(vy==-9999))=NaN;

switch opt
   case 'bestcoverage'
      for i=1:size(vx,3)
         nvx = numel(find(~isnan(vx(:,:,i))));
         nvy = numel(find(~isnan(vy(:,:,i))));
         nvalid(i) = min(nvx,nvy);
      end
      
      [~, loadidx] = max(nvalid);
      disp(['   -- Rosenau Landsat velocity: selecting ' datestr(dn(id1t+loadidx-1))]);
      vx = vx(:,:,loadidx);
      vy = vy(:,:,loadidx);

   case 'bestcoveragemonth'
      [~, months] = datevec(dn);
      checkidx = months == month;
      checkidx = find(checkidx(id1t:id2t));
      for i=checkidx
         nvx = numel(find(~isnan(vx(:,:,i))));
         nvy = numel(find(~isnan(vy(:,:,i))));
         nvalid(i) = min(nvx,nvy);
      end
      
      [~, loadidx] = nanmax(nvalid);
      disp(['   -- Rosenau Landsat velocity: selecting ' datestr(dn(id1t+loadidx-1))]);
      vx = vx(:,:,loadidx);
      vy = vy(:,:,loadidx);


end

vx = vx';
vy = vy';
vel = sqrt( vx.^2 + vy.^2 );

disp(['   -- Rosenau Landsat velocity: interpolating ']);
veli = InterpFromGrid(xdata,ydata,vel,double(X),double(Y));

