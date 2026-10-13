potential[r_,ell_]:=With[{lam=(ell-1)(ell+2)/2},
 (1-2/r)(2 lam^2(lam+1)r^3+6 lam^2 r^2+18 lam r+18)/(r^3(lam r+3)^2)];
tortoise[r_]:=r+2 Log[r/2-1];
radialSegment[w_,ell_,r0_,r1_,initial_,sign_]:=Module[{p,q,u,dr=r1-r0,rr,ff,vv,sol},
 rr=r0+u dr;ff=1-2/rr;vv=potential[rr,ell];
 sol=NDSolveValue[SetPrecision[{p'[u]==dr q[u],q'[u]==dr(vv p[u]-(2 I sign w ff+ff 2/rr^2)q[u])/ff^2,
 p[0]==initial[[1]],q[0]==initial[[2]]},40],{p,q},{u,0,1},WorkingPrecision->40,AccuracyGoal->22,PrecisionGoal->22,MaxSteps->Infinity];
 Through[sol[1]]];
wronskian[w0_,ell_,radius0_]:=Module[{w=SetPrecision[w0,40],radius=SetPrecision[radius0,40],rh=N[2+10^-8,40],rm=N[3,40],ri,h1,h,up,a1,a2,v3,lam,psih,dh,psiu,du},
 lam=(ell-1)(ell+2)/2;h1=(potential[rh,ell]/(1-2/rh)/2)/(1/4-I w);
 h=radialSegment[w,ell,rh,rm,{1+h1(rh-2),h1},-1];
 ri=rm+radius Exp[-I Arg[w]];a1=I ell(ell+1)/(2 w);v3=-4 lam-10-12/lam;
 a2=(a1(2+4 I w-ell(ell+1))-v3)/(4 I w);
 up=radialSegment[w,ell,ri,rm,{1+a1/ri+a2/ri^2,-a1/ri^2-2 a2/ri^3},1];
 psih=Exp[-I w tortoise[rm]]h[[1]];dh=Exp[-I w tortoise[rm]](h[[2]]-I w h[[1]]/(1-2/rm));
 psiu=Exp[I w tortoise[rm]]up[[1]];du=Exp[I w tortoise[rm]](up[[2]]+I w up[[1]]/(1-2/rm));
 (1-2/rm)(psih du-psiu dh)];
quadraticRatio[ell_,ell1_,ell2_,n1_,n2_]:=Module[{row,keys},
 keys=<|"l"->ell,"l1"->ell1,"l2"->ell2,"n1"->n1,"n2"->n2,"p"->0,"p1"->0,"p2"->0,"mirror1"->0,"mirror2"->0|>;
 row=SelectFirst[csvRows["source/external_inputs/QuadraticQNM/QuadraticRatio.csv"],Lookup[#,Keys[keys]]==Values[keys]&];
 row["ReRatioNL"]+I row["ImRatioNL"]];
angular[ell_,ell1_,ell2_,m1_,m2_]:=N[(-1)^(m1+m2) ThreeJSymbol[{ell,-m1-m2},{ell1,m1},{ell2,m2}]];
QCSTwoPole[radii_:{500.,1000.,2000.,5000.}]:=Module[{p,c,zdelta,wp,w4,b4,nu,raw,term0,term1,full,rows={},matrices,common,report},
 Do[Print["Two-pole radial integration: radius ",radius];p=ConstantArray[0.I,2];c=ConstantArray[0.I,{2,2}];zdelta=c;
 Do[wp=First[mode[ell,0,0.]];nu=2 wp;raw=quadraticRatio[4,ell,ell,0,0];
 p[[ell-1]]=raw angular[4,ell,ell,2,2](-nu^2/2)wronskian[nu,4,radius];
 Do[{w4,b4}=mode[4,n,0.];zdelta[[n+1,ell-1]]=(b4/(-w4^2/2))/(nu-w4),{n,0,1}],{ell,{2,3}}];
 Do[wp=First[mode[ell,0,0.]];Do[w4=First[mode[4,n,0.]];nu=w4+wp;
 c[[ell-1,n+1]]=quadraticRatio[ell+4,ell,4,0,n]angular[ell+4,ell,4,ell,4](-nu^2/2)wronskian[nu,ell+4,radius],{n,0,1}],{ell,{2,3}}];
 term0=Outer[Times,c[[All,1]],p zdelta[[1]]];term1=Outer[Times,c[[All,2]],p zdelta[[2]]];full=term0+term1;
 AppendTo[rows,<|"ray_radius"->radius,"isolated"->rankMetrics[term0],"two_pole"->rankMetrics[full],"overtone_fraction"->fro[term1]/fro[term0]|>],{radius,radii}];
 common=term0[[1,1]];matrices=<|"n0_direct_vertices"->term0/common,"n1_direct_vertices"->term1/common,"n01_direct_vertices"->full/common|>;
 report=<|"convergence"->rows,"matrices"->matrices,"reference_relative_matrix_error"->fro[matrices["n01_direct_vertices"]-referenceMatrices[]["n01_direct_vertices"]]/fro[referenceMatrices[]["n01_direct_vertices"]]|>;
 save["two_pole",report];Export[path["results/two_pole_matrices.wxf"],matrices,"WXF"];report];
