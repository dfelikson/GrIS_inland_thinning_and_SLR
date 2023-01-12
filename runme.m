steps = [4];

%% Setup %%
glacier = 'KAK'; glacier_epoch = datetime(1985,07,23);
%glacier = 'KLG'; glacier_epoch = 1999; %glacier_epoch = 1981;
prefix = '';

% To launch a job remotely and not wait, set the following:
%  md.cluster.interactive = 0;
%  md.settings.waitonlock = nan; or = 0;
clusterName = ''; % localhost
clusterName = 'oibserve';
%clusterName = 'discover';
%clusterName = 'ls5';
%clusterName = 'melt';
switch clusterName %%{{{
   case {'','gs15serac'}
      cluster = generic('name', oshostname(), 'np', 4);
      waitonlock = Inf; %nan;

   case 'oibserve'
      cluster = generic('name', 'gs615-oibserve.ndc.nasa.gov', 'np', 24, 'interactive', 0, ...
         'login', 'dfelikso', ...
         'codepath', '/home/dfelikso/Software/ISSM/trunk-jpl/bin', ...
         'etcpath', '/home/dfelikso/Software/ISSM/trunk-jpl/etc', ...
         'executionpath', '/home/dfelikso/Projects/GrIS_Calibrated_SLR/ISSM/execution');
      cluster.interactive = 0; %1;
      waitonlock = Inf; %nan;

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
if perform(org,'Lcurve'),% {{{ STEP 3

   md = loadmodel(org,'Param');

   % Control inversion -- general
   md.inversion=m1qn3inversion(md.inversion);
   md.inversion.iscontrol=1;
   md.verbose=verbose('solution',false,'control',true,'convergence',false);

   % Control -- other
   md.transient.issmb = 0;
   md.transient.isthermal = 0;

   % Cost functions
   md.inversion.cost_functions=[101 103 501]; %Abs, Log, reg
   md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,length(md.inversion.cost_functions));
   md.inversion.cost_functions_coefficients(:,1)=2000;
   md.inversion.cost_functions_coefficients(:,2)=40;
   %Computed in a loop below
   %md.inversion.cost_functions_coefficients(:,3)=.2*50^-3;

   % %Remove obs where the front from the velocities are upstream of our current front
   % filename = ['Exp/' glacier '_velfront.exp'];
   % if exist(filename,'file'),
   %  disp(['Correcting cost functions for front inconsistencies']);
   %  pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,2));
   %  md.friction.coefficient(pos)=min(md.friction.coefficient(pos),100);
   % end
   
   % Where vel==0, set coefficients to 0 (i.e., don't try to match this in model)
   disp(['Removing vel==0 obs from inversion']);
   pos = find(md.inversion.vel_obs == 0);
   md.inversion.cost_functions_coefficients(pos,:) = 0;

   % Controls
   md.inversion.control_parameters={'FrictionCoefficient'};
   %md.inversion.maxsteps=50;
   %md.inversion.maxiter =50;
   md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
   md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1);
   md.inversion.control_scaling_factors=1;

   % Set basal friction coefficient initial guess to something low at front %%{{{
   filename = ['Exp/' glacier '_coeffront.exp'];
   %if ~exist(filename,'file'),
   %   plotmodel(md,'data',md.friction.coefficient,'mask',md.mask.ice_levelset<0)
   %   exptool(filename)
   %end
   if exist(filename,'file'),
      disp(['Correcting basal friction coefficient initial guess for front inconsistencies']);
      flags = ContourToNodes(md.mesh.x,md.mesh.y,filename,2);
      %flags = md.inversion.vel_obs == 0;
      pos1 = find(flags); pos2 = find(~flags);
      %md.friction.coefficient(pos1,:) = 50;
      %md.friction.coefficient(pos1,:) = 100;
      md.friction.coefficient(pos1) = 1;
      md.inversion.max_parameters(pos1)= 1;
      %md.friction.coefficient(pos1,:) = 40;
      %md.friction.coefficient(pos1,:) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),md.friction.coefficient(pos2,:),md.mesh.x(pos1),md.mesh.y(pos1));

      % % Special case: KLG
      % % Although observed velocities at the front look reasonable, a "ridge" of high friction comes out of the inversion. This is meant to
      % % neglect velocities along this erroneous ridge in the inversion.
      % if glacier == 'KLG'
      %    md.inversion.cost_functions_coefficients(pos1,1) = 0;
      %    md.inversion.cost_functions_coefficients(pos1,2) = 0;
      % end
   end
   %%}}}

   % %Fix friction coefficient
   % filename = ['Exp/' glacier '_fixFrictionCoefficient.exp'];
   % if exist(filename,'file'),
   %    disp(['Ignoring manually selected velocities in the inversion using ' filename]);
   %    pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,1));
   %    %md.inversion.min_parameters(pos)=md.friction.coefficient(pos);
   %    md.inversion.max_parameters(pos)=10; %md.friction.coefficient(pos);
   %    %md.inversion.cost_functions_coefficients(pos,1) = 0;
   %    %md.inversion.cost_functions_coefficients(pos,2) = 0;
   %    %md.inversion.cost_functions_coefficients(pos,3) = 0;
   % end

   %Additional parameters
   md.stressbalance.restol=0.01;
   md.stressbalance.reltol=0.1;
   md.stressbalance.abstol=NaN;
   %md.stressbalance.requested_outputs={'default','DeviatoricStressxx','DeviatoricStressyy','DeviatoricStressxy'}

   % Go solve
   md.cluster=cluster;
   md.settings.waitonlock = waitonlock;
   % NOTE: For KAK, 2e-06 is the winner!
   % NOTE: For KLG ...
   coeffs = [1e-06, 2e-06:1e-05:8e-05, 8e-05, 4e-03];
   for i = 1:length(coeffs)
      md.inversion.cost_functions_coefficients(:,3)=coeffs(i);
      md=solve(md,'Stressbalance');
      md.results.Lcurve(i) = md.results.StressbalanceSolution;
   end

   savemodel(org,md);
end%}}}
if perform(org,'Inversion'),% {{{ STEP 4

   md = loadmodel(org,'Param');

   % Control inversion -- general
   md.inversion=m1qn3inversion(md.inversion);
   md.inversion.iscontrol=1;
   md.verbose=verbose('solution',false,'control',true,'convergence',false);

   % Control -- other
   md.transient.issmb = 0;
   md.transient.isthermal = 0;

   % Cost functions
   md.inversion.cost_functions=[101 103 501]; %Abs, Log, reg
   md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,length(md.inversion.cost_functions));
   md.inversion.cost_functions_coefficients(:,1)=2000;
   md.inversion.cost_functions_coefficients(:,2)=40;
   % NOTE: For KAK, 2e-06 is the winner
   md.inversion.cost_functions_coefficients(:,3)=2e-06;

   % %Remove obs where the front from the velocities are upstream of our current front
   % filename = ['Exp/' glacier '_velfront.exp'];
   % if exist(filename,'file'),
   %  disp(['Correcting cost functions for front inconsistencies']);
   %  pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,2));
   %  md.friction.coefficient(pos)=min(md.friction.coefficient(pos),100);
   % end
   
   % Where there is ice but the inversion cost function coefficient was zero, use nearest neighbor interpolation for inversion.vel_obs
   pos1 = find(md.mask.ice_levelset<0 & md.inversion.vel_obs==0);
   pos2 = find(md.mask.ice_levelset<0 & md.inversion.vel_obs >0);
   x_valid = md.mesh.x(pos2);
   y_valid = md.mesh.y(pos2);
   vx_obs_valid = md.inversion.vx_obs(pos2);
   vy_obs_valid = md.inversion.vy_obs(pos2);
   vel_obs_valid = md.inversion.vel_obs(pos2);
   for i = 1:length(pos1)
      pos = pos1(i);
      [~, nearest_neighbor_pos] = min( sqrt( (md.mesh.x(pos)-x_valid).^2 + (md.mesh.y(pos)-y_valid).^2 ) );
      md.inversion.vx_obs(pos) = vx_obs_valid(nearest_neighbor_pos);
      md.inversion.vy_obs(pos) = vy_obs_valid(nearest_neighbor_pos);
      md.inversion.vel_obs(pos) = vel_obs_valid(nearest_neighbor_pos);
   end

   % Controls
   md.inversion.control_parameters={'FrictionCoefficient'};
   md.inversion.maxsteps=50;
   md.inversion.maxiter =50;
   md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
   md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1);
   md.inversion.control_scaling_factors=1;

   %Additional parameters
   md.stressbalance.restol=0.01;
   md.stressbalance.reltol=0.1;
   md.stressbalance.abstol=NaN;
   %md.stressbalance.requested_outputs={'default','DeviatoricStressxx','DeviatoricStressyy','DeviatoricStressxy'}

   % Go solve
   md.cluster=cluster;
   md.settings.waitonlock = waitonlock;
   md=solve(md,'Stressbalance');

   savemodel(org,md);
end%}}}

if perform(org,'Transient'),% {{{ STEP 5

   % NOTE: Question to look into: how sensitive is the stress balance to initial velocity?

   md = loadmodel(org,'Inversion');
   md = parameterize(md,'Par/KAK_AERODEM_1985.par');

   % Remove duplicate spclevelset times
   pos = find(diff(md.levelset.spclevelset(end,:))==0);
   md.levelset.spclevelset(:,pos) = [];

   md.timestepping.start_time = year(glacier_epoch) + day(glacier_epoch, 'dayofyear') / day(datetime(year(glacier_epoch), 12, 31), 'dayofyear');
   % To 2015
   md.timestepping.final_time = 2015;
   md.settings.output_frequency = (1/md.timestepping.time_step)/8; % forward run to 2015
   % To 2100
   %md.timestepping.final_time = 2100;
   %md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2100

   % Set the requested outputs
   md.transient.requested_outputs={'default','IceVolume'};
   md.stressbalance.requested_outputs={'default'};

   % Go solve
   %if contains(cluster.name, 'ls5')
   %   cluster.time = 5*3600;
   %   fprintf('Check cluster params:\n')
   %   cluster
   %   s = input('-> Continue (y/n)?','s');
   %   if ~strcmpi(s(1),'y')
   %      return
   %   end
   %end
   md.verbose.solution=1;
   md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md=solve(md,'transient');

   savemodel(org,md);
end%}}}

