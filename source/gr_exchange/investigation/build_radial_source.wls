base = StringReplace[DirectoryName[$InputFileName], "\\" -> "/"] <> "/";
f=1-2/r;vr[j_]:=f*(j*(j+1)/r^2-6/r^3);
cc[j_]:=(j+2)*(j-1)*((j+2)*(j-1)+2)/6;
hh[j_]:=cc[j]+12*(r-2)/(r^2*((j+2)*(j-1)*r+6));
dr[expr_,j_]:=D[expr,r]+D[expr,y]*yp+D[expr,yp]*((vr[j]-w^2)*y/f^2-D[f,r]*yp/f);
Do[
 If[j==2,fs=Get[base<>"source_220_coeffs.wl"];az=hh[2]/(cc[2]+2*I*w);bz=2*f/(cc[2]+2*I*w);z=az*y+bz*yp;zp=dr[z,2];sz=fs.{z^2,z*zp,zp^2},
 ss=Get[base<>"source_320_specialized.wl"];
 sz=ss/.{\[Psi]1l1o[r]->y,\[Psi]1l2o[r]->y,Derivative[1][\[Psi]1l1o][r]->yp,Derivative[1][\[Psi]1l2o][r]->yp,\[Omega]1->w,\[Omega]2->w}];
 sz=Together[sz];Print["Z source ",j," ",LeafCount[sz]];
 Put[Together/@{Coefficient[sz,y,2],Coefficient[Coefficient[sz,y,1],yp,1],Coefficient[sz,yp,2]},base<>"zsource_rwparents_"<>ToString[j]<>".wl"];
 sr=Together[(hh[4]*sz-2*f*dr[sz,j])/(cc[4]-4*I*w)];
 la=I*w-2*I*w/(r-2)+4*I*w/r;
 ss=Together[(r^3/(r-2))*(sr/.{y->pp,yp->la*pp+(2/r^2)*dp})];
 gg=Together/@({Coefficient[ss,pp,2],Coefficient[Coefficient[ss,pp,1],dp,1],Coefficient[ss,dp,2]}/.r->2/(1-u));
 Put[gg,base<>"leaver_source_"<>ToString[j]<>".wl"];
 Print["Leaver source ",j," ",LeafCount/@gg];
 rec=Together[Get[base<>"reconstruct_"<>If[j==2,"220","320"]<>".wl"]/.{\[Omega]1->w,\[Omega]2->w}];
 Put[rec,base<>"reconstruct_simple_"<>ToString[j]<>".wl"];Print["Reconstruction ",rec];
,{j,{2,3}}];Exit[];
