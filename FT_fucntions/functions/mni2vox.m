function VOX=mni2vox(M,MNI)

%function VOX=mni2vox(M,MNI)

%MNI-space to Voxel-space

 T=M(1:3,4);

 M=M(1:3,1:3);

 for i=1:3

MNI(i,:)=MNI(i,:)-T(i);

 end

 VOX=round(inv(M)*MNI);

return