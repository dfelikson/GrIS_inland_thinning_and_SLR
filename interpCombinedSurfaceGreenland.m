function [surf_combined, surf_aero, surf_gimp, surf_id] = interpCombinedSurfaceGreenland(x, y)

   disp('      -- reading Howat surface');
   surf_gimp = interpBedmachineGreenland(x,y,'surface');
   
   disp('      -- reading AERO DEM surface');
   surf_aero = interpAEROdem(x,y);

   filename = ['Exp/' evalin('base','glacier') '_AEROoverlapGIMP.exp'];
   in = ContourToNodes(x,y,filename,2);

   mean_diff = mean(surf_aero(find(in)) - surf_gimp(find(in)),'omitnan');
   surf_gimp_shifted = surf_gimp + mean_diff;
   disp(sprintf('        -- Shifting GIMP surface by %4.1f meters', mean_diff));

   % Where do we not have AERO surface
   surf_id = ones(size(x));
   disp('        -- Combining surfaces');
   pos = isnan(surf_aero);
   surf_combined = surf_aero;
   surf_combined(pos) = surf_gimp_shifted(pos);
   surf_id(pos) = 2;

end % main function

