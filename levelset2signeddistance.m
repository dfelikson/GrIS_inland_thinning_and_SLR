function levelset_out = levelset2signeddistance(md, levelset_in)

	levelset_out = levelset_in;

	if size(levelset_in,2)>1,
      for i=1:size(levelset_in,2)
         levelset = levelset_in(1:end-1,i);
         pos      = find(levelset<0);

         if exist('TEMP.exp','file'), delete('TEMP.exp'); end
			isoline(md,levelset,'output','TEMP.exp');
         levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'TEMP.exp'));
         delete('TEMP.exp');
         levelset(pos) = - levelset(pos);

         levelset_out(1:end-1,i) = levelset;
      end
	else
		levelset = levelset_in;
		pos      = find(levelset<0);

		if exist('TEMP.exp','file'), delete('TEMP.exp'); end
		isoline(md,levelset,'output','TEMP.exp');
		levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'TEMP.exp'));
		delete('TEMP.exp');
		levelset(pos) = - levelset(pos);

		levelset_out = levelset;
	end

end % main function

