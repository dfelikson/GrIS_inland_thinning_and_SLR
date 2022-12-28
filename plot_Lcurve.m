function plot_Lcurve(md, regularization_coefficients, fig_handle)

   % if A is md.results.StressbalanceSolution.J
   % if A(:,5) and A(:,6) are the "misfit" cost functions (log and abs for example)
   %    A(:,7) is the Thikhonov regularization term
   %    A(:,1) is the weight applied to A(:,7)


	if ~exist('fig_handle','var')
		fig_handle = figure;
	else
		figure(fig_handle.Number);
		hold on;
	end

   % Pull out the relevant terms from J
   for i = 1:numel(md.results.Lcurve)
      Jres(i) = sum(md.results.Lcurve(i).J(end,1:2));
      Jreg(i) = md.results.Lcurve(i).J(end,3) / regularization_coefficients(i);
   end

   loglog(Jres, Jreg, '-s', 'Color', [.3 .8 .4], 'MarkerSize', 6, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'LineWidth', 2)
   
   voffset=zeros(size(Jres));
   voffset(end)=0;
   hoffset=Jres/15;
   text(Jres+hoffset,Jreg+voffset,[repmat('\gamma_t = ',length(regularization_coefficients),1) num2str(regularization_coefficients(:),'%2.0e')],...
         'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','Middle');

   %xlim([7*10^3 10^5]);
   %ylim([4*10^8 7*10^10]);
   
   xlabel('$\mathcal J_{res} \left({\bf v}\right)$', 'interpreter', 'latex', 'FontSize', 15);
   ylabel('$\mathcal J_{reg} \left(\alpha\right)$',  'interpreter', 'latex', 'FontSize', 15);
   set(gcf,'color','w');
   set(gca,'FontSize',14);

   glacier = evalin('base','glacier');
   export_fig(['Graphics/' glacier '_Lcurve.pdf'], '-nocrop');

end % main function

