rwCoefficients[nmax_]:=Module[{w=.37367168441804183579-.08896231568893569828 I,rho,a,al,be,ga},
 rho=-2 I w;a=ConstantArray[0.I,nmax+1];a[[1]]=1.;
 Do[al=n^2+(2 rho+2)n+2 rho+1;be=-(2 n^2+(8 rho+2)n+8 rho^2+4 rho+3);ga=n^2+4 rho n+4 rho^2-4;
 a[[n+2]]=-(be a[[n+1]]+If[n==0,0,ga a[[n]]])/al,{n,0,nmax-1}];a];
polyval[c_,u_]:=Fold[#1 u+#2&,Last[c],Reverse[Most[c]]];
parentZ[u_,a_]:=Module[{w=.37367168441804183579-.08896231568893569828 I,rho,y,yp,ypp,q,qp,la,lap,num,den,t,tp,z,zp},
 rho=-2 I w;y=polyval[a,u];yp=polyval[Range[Length[a]-1]Rest[a],u];ypp=polyval[Range[2,Length[a]-1]Range[Length[a]-2]Drop[a,2],u];
 q=u(1-u)^2/2;qp=(1-u)(1-3 u)/2;la=rho/u-rho/(1-u)-rho/(1-u)^2;lap=-rho/u^2-rho/(1-u)^2-2 rho/(1-u)^3;
 num=6 u(1-u)^2;den=14-6 u;t=4+num/den;tp=(6(1-u)(1-3 u)den+6 num)/den^2;
 z=(t y+2 q(yp+la y))/(4+2 I w);zp=(tp y+t yp+2 qp(yp+la y)+2 q(ypp+lap y+la yp))/(4+2 I w);
 {z,(1-u)^2(zp+la z)/2}];
rwOperator[u_]:=Module[{w=2(.37367168441804183579-.08896231568893569828 I),rho,q,qp,la,lap,invr,vz},
 rho=-2 I w;q=u(1-u)^2/2;qp=(1-u)(1-3 u)/2;la=rho/u-rho/(1-u)-rho/(1-u)^2;lap=-rho/u^2-rho/(1-u)^2-2 rho/(1-u)^3;
 invr=(1-u)/2;vz=u/(18+6 invr)^2(18^2 invr^2(20+6 invr)+36 invr^4(18+2 invr));
 {q^2,2 q^2 la+q qp,q^2(la^2+lap)+q qp la+w^2-vz}];
cauchy[values_,radius_,nmax_]:=Take[Fourier[values,FourierParameters->{-1,-1}],nmax+1]/radius^Range[0,nmax];
rwSeries[p2_,p1_,p0_,source_,nmax_,z0_]:=Module[{z=ConstantArray[0.I,nmax+1],diag,res,j},z[[1]]=z0;
 Do[diag=p2[[3]]k(k-1)+p1[[2]]k+p0[[1]];res=0.I;
 Do[j=k-m;If[j<k,res+=p0[[m+1]]z[[j+1]]];j=k-m+1;If[0<=j<k,res+=p1[[m+1]]j z[[j+1]]];
 j=k-m+2;If[0<=j<k,res+=p2[[m+1]]j(j-1)z[[j+1]]],{m,0,k}];z[[k+1]]=(source[[k+1]]-res)/diag,{k,1,nmax}];z];
QCSRWZ[]:=Module[{a=rwCoefficients[160],u=.95 Exp[2 Pi I Range[0,1023]/1024],zz,dz,source,sc,pc,p,h,z0,z,radial,tt,rows={},parentInf,funcs,rr},
 Block[{$Context="QCS`Private`",$ContextPath={"QCS`Private`","System`"}},Get[path["source/analysis/fctSourceEEE_Z_220x220.wl"]]];
 If[!And@@(NumericQ[N[#/.r->3.]]&/@{fctsource1Z,fctsource2Z,fctsource4Z}),Return[$Failed]];
 funcs=Function[Evaluate[{r}],Evaluate[#]]&/@{fctsource1Z,fctsource2Z,fctsource4Z};rr=2/(1-u);{zz,dz}=parentZ[u,a];
 source=(funcs[[1]]/@rr)zz^2+2(funcs[[2]]/@rr)zz dz+(funcs[[3]]/@rr)dz^2;
 sc=cauchy[source,.95,150];pc=cauchy[#,.95,152]&/@rwOperator[u];parentInf=Total[Take[a,151]];
 Do[p=rwSeries[pc[[1]],pc[[2]],pc[[3]],sc,n,0];h=rwSeries[pc[[1]],pc[[2]],pc[[3]],sc,n,1]-p;
 Do[z0=-(Conjugate[Take[h,-tail]].Take[p,-tail])/(Conjugate[Take[h,-tail]].Take[h,-tail]);z=p+z0 h;radial=Total[z]/parentInf^2;
 tt=Sqrt[10.]/2 radial+.01241195804083564-.1218015501258257 I;
 AppendTo[rows,<|"N"->n,"tail"->tail,"z0"->z0,"radial_ratio"->radial,"normalized_TT"->tt|>],{tail,{8,12,20}}],{n,{60,80,100,120,140}}];
 save["rwz",<|"series_table"->rows,"quoted_N80_tail12"->SelectFirst[rows,#["N"]==80&&#["tail"]==12&]|>];rows];
