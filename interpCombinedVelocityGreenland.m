function [velx_combined, vely_combined, vel_id] = interpCombinedVelocityGreenland(x, y, scalefactor_MEaSUREs, scalefactor_Landsat, order)
   if ~exist('order','var')
      order = 'MEaSUREs+Landsat';
   end

   disp('        -- Loading MEaSUREs velocity composite');
   [velx_MEaSUREs, vely_MEaSUREs] = interpJoughinCompositeGreenland(x,y);

   fprintf('        -- Scaling MEaSUREs velocity by %4.1f\n', scalefactor_MEaSUREs);
   velx_MEaSUREs_scaled = scalefactor_MEaSUREs .* velx_MEaSUREs;
   vely_MEaSUREs_scaled = scalefactor_MEaSUREs .* vely_MEaSUREs;

   disp('        -- Loading Landsat velocity');
   [velx_Landsat, vely_Landsat] = interpLandsatVelocityGreenland(x,y);

   fprintf('        -- Scaling Landsat velocity by %4.1f\n', scalefactor_Landsat);
   velx_Landsat_scaled = scalefactor_Landsat .* velx_Landsat;
   vely_Landsat_scaled = scalefactor_Landsat .* vely_Landsat;

   disp('        -- Combining velocities');
   glacier = evalin('base','glacier');
   if exist(['Exp/' glacier '_useLandsatVel.exp'])
      vel_id = ones(size(x));
      md = evalin('base','md');
      pos = find(ContourToMesh(md.mesh.elements, x, y, ['Exp/' glacier '_useLandsatVel.exp'], 'node', 1));
      velx_combined = velx_MEaSUREs_scaled;
      vely_combined = vely_MEaSUREs_scaled;
      velx_combined(pos) = velx_Landsat_scaled(pos);
      vely_combined(pos) = vely_Landsat_scaled(pos);
      vel_id(pos) = 2;

   else
      switch order
         case 'MEaSUREs+Landsat'
            vel_id = ones(size(x));
            % Where do we not have MEaSUREs velocities (in the ice mask)
            pos = isnan(velx_MEaSUREs);
            velx_combined = velx_MEaSUREs_scaled;
            vely_combined = vely_MEaSUREs_scaled;
            velx_combined(pos) = velx_Landsat_scaled(pos);
            vely_combined(pos) = vely_Landsat_scaled(pos);
            vel_id(pos) = 2;

         case 'Landsat+MEaSUREs'
            vel_id = 2 * ones(size(x));
            % Where do we not have Landsat velocities (in the ice mask)
            pos = isnan(velx_Landsat);
            velx_combined = velx_Landsat_scaled;
            vely_combined = vely_Landsat_scaled;
            velx_combined(pos) = velx_MEaSUREs_scaled(pos);
            vely_combined(pos) = vely_MEaSUREs_scaled(pos);
            vel_id(pos) = 1;

      end
   end

end % main function

