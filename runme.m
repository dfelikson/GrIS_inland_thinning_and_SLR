steps = [1];

%% Setup %%
glacier = 'KAK'; glacier_epoch = 1985;
%glacier = 'KLG'; glacier_epoch = 1999; %glacier_epoch = 1981;
prefix = '';

% To launch a job remotely and not wait, set the following:
%  md.cluster.interactive = 0;
%  md.settings.waitonlock = nan; or = 0;
clusterName = ''; % localhost
%clusterName = 'oibserve';
%clusterName = 'discover';
%clusterName = 'ls5';
%clusterName = 'melt';
switch clusterName %%{{{
   case {'','gs15serac'}
      cluster = generic('name', oshostname(), 'np', 2);
      waitonlock = 10; %nan;

   case 'oibserve'
      cluster = generic('name', 'gs615-oibserve.ndc.nasa.gov', 'np', 12, 'interactive', 0, ...
         'login', 'dfelikso', ...
         'codepath', '/home/dfelikso/Software/ISSM/trunk-jpl/bin', ...
         'etcpath', '/home/dfelikso/Software/ISSM/trunk-jpl/etc', ...
         'executionpath', '/home/dfelikso/Projects/GrIS_Calibrated_SLR/ISSM/execution');
      cluster.interactive = 0; %1;
      waitonlock = nan;

   case 'discover'
      cluster=discover;
      cluster.name='discover.nccs.nasa.gov';
      cluster.login='dfelikso';
      cluster.numnodes=1;
      cluster.cpuspernode=16;
      cluster.time=1*60;
      cluster.interactive=0;
      cluster.processor='sand';
      cluster.queue='allnccs';
      cluster.codepath='/discover/nobackup/dfelikso/Software/ISSM/trunk-jpl/bin';
      cluster.executionpath='/discover/nobackup/dfelikso/Software/ISSM/trunk-jpl/execution';
      cluster.email='denis.felikson@nasa.gov';
      waitonlock = nan;

   case 'ls5'
	   % Lonestar5 has 24 cores per node
	   % 'time' is the total run time for the job
      cluster = lonestar('project', 'Central-West-GrIS', ...
                         'login','byaa735', ...
                         'codepath','/home1/00761/byaa735/Software/ISSM/trunk-jpl/bin', ...
                         'executionpath','/scratch/00761/byaa735/ISSM/execution',...
                         'email','denis.felikson@utexas.edu');
      waitonlock = 0;
   
   case 'melt'
      % NOTE: To run interactively, set up, e.g.:
      %     waitonlock = 1;
      %     interactive = 36000;
      % NOTE: On melt, must run long jobs (e.g., 100 yr transient) with:
      %     waitonlock = 0;
      %     interactive = 0;
      % This uploads Greenland.bin and prints the launch command to the Matlab command window. The user
      % can then ssh to melt manually and launch this from within a screen. If I can figure out how to
      % start a screen, launch a command, and detatch from within ssh, this can be automated.
      issm_dir = '/home/student/denis/Software/ISSM/trunk-jpl';
      %cluster=generic('name','melt.ig.utexas.edu','np',4,'shell','csh','login','denis',...
      %   'codepath',[issm_dir '/bin'],...
      %   'executionpath',[issm_dir '/execution'],...
      %   'etcpath',[issm_dir '/etc']);
      cluster=melt('name','melt.ig.utexas.edu','np',6,'login','denis',...
         'codepath',[issm_dir '/bin'],...
         'executionpath','/disk/student/denis/Software/ISSM/trunk-jpl/execution',...
         'etcpath',[issm_dir '/etc']);
      cluster.interactive = 1;
      waitonlock = 36000;
end
%%}}}


%% Model parameters %%
org=organizer('repository',['./Models'],'prefix',[glacier prefix '_'],'steps',steps);

%% Processing %%
set(0,'DefaultFigureWindowStyle' , 'normal')

if ~exist(['Exp/' glacier '.exp'],'file'), % {{{
   disp(['Domain outline needed: Exp/' glacier '.exp.'])
   
   s = input('Is there a shapefile that can be used for the model domain (y/n)?','s');
   if strcmpi(s(1),'y')
      [file,path] = uigetfile('/Users/denisfelikson/Research/Projects/ModeledInlandThinning/Data/Glacier model domains/*.shp');
      shp2exp(fullfile(path,file),['Exp/' glacier '.exp']);
   else
      disp('Exiting ... rerun runme')
      return
   end
   %draw_glacier_domain(glacier, 'GIMP', 'MEaSUREs');
end
%}}}

if perform(org,'Mesh'),% {{{ STEP 1
	
   %% Mesh sizing{{{
   switch glacier
      case 'UMI'
         triangleresolution = 300;
         hmin = 150;
         hmax = 10000;
      case 'HEL'
         triangleresolution = 250; %500;
         hmin = 150; %300;
         hmax = 500; %10000;
      case 'KLG'
         triangleresolution = 1000;
         hmin = 300;
         hmax = 10000;
      case 'KAK'
         triangleresolution = 300;
         hmin = 250;
         hmax = 10000;
      otherwise
         triangleresolution = 300;
         hmin = 150;
         hmax = 10000;
   end
   %%}}}
   fprintf('Creating mesh with the following parameters:\n   %20s = %5d\n   %20s = %5d\n   %20s = %5d\n', 'triangleresolution', triangleresolution, 'hmin', hmin, 'hmax', hmax);
   s = input('-> Continue (y/n)?','s');
   if ~strcmpi(s(1),'y')
      return
   end
	md=triangle(model,['./Exp/' glacier '.exp'],triangleresolution);

   % Read velocity (for mesh refinement)
	%if ~exist('vel','var'),
      disp('Reading Joughin composite velocities');
		[velx, vely] = interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
      %disp(['Reading Joughin year ' num2str(2000) ' velocities']);
		%[velx, vely] = interpJoughin(md.mesh.x,md.mesh.y,2000);
		vel  = sqrt(velx.^2+vely.^2);
	%end

	%Adapt mesh
	disp('Optimizing mesh');
   
   %Refine mesh beyond terminus (because BAMG refinement above uses velocities that -> 0 where there's still ice in 1985)
   filename = ['Exp/' glacier '_refineFront.exp'];
   if ~exist(filename,'file'),
      s = input('Is there a shapefile that can be used to refine the front (y/n)?','s');
      if strcmpi(s(1),'y')
         [file,path] = uigetfile('/Users/denisfelikson/Research/Projects/ModeledInlandThinning/Data/Glacier model domains/*.shp');
         shp2exp(fullfile(path,file),['Exp/' glacier '_refineFront.exp']);
      else
         filename = ['Exp/' glacier '_Front0.5.exp'];
      end
   end
   
   if ~exist(filename,'file'),
      fprintf(['\n\033[' '103;30' 'm   WARNING: mesh has not been refined beyond present-day front! \033[0m \n\n']);
      disp(['  -- Refining mesh using velocities'])
      md=bamg(md,'hmin',hmin,'hmax',hmax,'field',vel,'err',2);
   else
      disp(['  -- Refining mesh using velocities and ' filename])
      hmaxVertices=NaN*ones(md.mesh.numberofvertices,1);
      in=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,filename,'node',1);
      hmaxVertices(find(in))=hmin;
      md=bamg(md,'hmax',hmax,'hmin',hmin,'err',2,'field',vel,'hmaxVertices',hmaxVertices);
      % %Reload velocities (if needed later)
      % [velx, vely] = interpJoughin(md.mesh.x,md.mesh.y,2000);
      % vel  = sqrt(velx.^2+vely.^2);
   end

	%[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,45,70);
	md.mesh.epsg=3413;

	savemodel(org,md);
end %}}}
if perform(org,'Param'),% {{{ STEP 2

	md=loadmodel(org,'Mesh');
	md=parameterize(md,'Par/Greenland.par');
	
   md=setflowequation(md,'SSA','all');
   
   % Weaken shear margins
   filename = ['Exp/' glacier '_shearmargins.exp'];
   if exist(filename, 'file')
      disp(['Weakening shear margins by 60% using: ' filename]);
      pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,1));
      md.materials.rheology_B(pos) = 0.60 .* md.materials.rheology_B(pos);
   end

	savemodel(org,md);
end%}}}
