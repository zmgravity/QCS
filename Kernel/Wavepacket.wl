relative[a_,b_]:=fro[a-b]/fro[b];
packetSeed[p_,kind_]:=With[{z=p["z"]},If[kind=="original",Join[.25 z^2(1-z)^2,ConstantArray[0.,Length[z]]],
 Join[.25 z^2(1-z)^2 Exp[-((z-.55)/.22)^2+5 I z],.06 z^2(1-z)^2 Cos[3 Pi z]Exp[-2 I z]]]];
trialSpace[p_,order_]:=Module[{q=First[basis[p["z"],order]],zz},q=Transpose[(#/Norm[#])&/@Transpose[q]];
 zz=zero[p["n"],order];ArrayFlatten[{{q,zz},{zz,q}}]];
packetInverse[p_,radiation_,t_,o_,order_:26,window_:60.,channels_:{1,2},noise_:0.]:=Module[{s,ix,a,b,u,sv,v,keep,c,eta,x},
 s=trialSpace[p,order];ix=Flatten[Position[t,z_/;z<=window+10^-10]];a=Flatten[Map[#.s&,o[[ix,channels]]],1];b=Flatten[radiation[[ix,channels]]];
 If[noise>0,SeedRandom[703];eta=normalComplex[Dimensions[b]];b+=noise Norm[b]eta/Norm[eta]];
 {u,sv,v}=SingularValueDecomposition[a,Tolerance->0];sv=Diagonal[sv];keep=Flatten[Position[sv,z_/;z>10^-10 First[sv]]];
 c=v[[All,keep]].((ConjugateTranspose[u[[All,keep]]].b)/sv[[keep]]);x=s.c;
 {x,<|"order"->order,"window"->window,"channels"->channels,"noise"->noise,"retained"->Length[keep],
 "radiation_fit_error"->relative[a.c,b],"radiation_all_error"->relative[Map[#.x&,o],radiation]|>}];
QCSWavepacket[]:=Module[{t=N[Range[0,700]/10],cache=<||>,p,dt,o,f,x,r,refs,ref,pred,initial,coarse,checks={},scans={},saved=<||>,
 rec,diag,n,order,control={},ob,modes,vals,ss,fit,design,werror,report},
 Do[p=QCSProblem[2,nn];dt=If[nn==32,.005,.0025];Print["Wavepacket response operators: n=",nn];
 o=response[p["m4"],boundary[p],t,dt];f=Table[With[{pp=QCSProblem[ell,nn]},response[pp["T"],pp["obs"],t,dt][[All,1,1;;2 nn]]],{ell,{2,3}}];
 AssociateTo[cache,nn-><|"O"->o,"F"->f|>],{nn,{32,48,64}}];
 Do[x=packetSeed[QCSProblem[2,64],kind];r=Map[#.x&,cache[64]["O"]];refs={};
 Do[p=QCSProblem[ell,64];initial=Join[x,ConstantArray[0.I,130]];
 ref=Last[propagate[p["T"],p["obs"],initial,70.,.1,.0025]];AppendTo[refs,ref];
 coarse=Last[propagate[p["T"],p["obs"],initial,70.,.1,.005]];
 AppendTo[checks,<|"seed"->kind,"ell"->ell+4,"row_kernel_vs_direct"->relative[cache[64]["F"][[ell-1]].x,ref],
 "dt_005_vs_0025"->relative[coarse,ref],"radial_32_vs_64"->relative[cache[32]["F"][[ell-1]].packetSeed[QCSProblem[2,32],kind],ref],
 "radial_48_vs_64"->relative[cache[48]["F"][[ell-1]].packetSeed[QCSProblem[2,48],kind],ref]|>],{ell,{2,3}}];
 AssociateTo[saved,kind-><|"radiation"->r,"reference"->Transpose[refs]|>];
 Do[{n,order}=spec;{rec,diag}=packetInverse[QCSProblem[2,n],r,t,cache[n]["O"],order];
 werror=Table[relative[cache[n]["F"][[j]].rec,refs[[j]]],{j,2}];
 AppendTo[scans,Join[diag,<|"seed"->kind,"n_inverse"->n,"waveform_error"->werror|>]];
 If[n==32&&order==26,AssociateTo[saved,kind->Join[saved[kind],<|"prediction"->Transpose[Table[cache[n]["F"][[j]].rec,{j,2}]]|>]]],
 {spec,{{32,18},{32,22},{32,26},{48,26},{48,30},{48,34}}}],{kind,{"original","wavepacket"}}];
 r=saved["wavepacket"]["radiation"];refs=Transpose[saved["wavepacket"]["reference"]];
 Do[{rec,diag}=packetInverse[QCSProblem[2,32],r,t,cache[32]["O"],26,window,channels,noise];
 AppendTo[scans,Join[diag,<|"seed"->"wavepacket","n_inverse"->32,"waveform_error"->Table[relative[cache[32]["F"][[j]].rec,refs[[j]]],{j,2}]|>]],
 {window,{10.,20.,40.,60.,70.}},{channels,{{1,2},{2}}},{noise,{0.,10^-6,10^-4}}];
 report=<|"checks"->checks,"inverse_scans"->scans,"random_generator"->"Wolfram default generator with seed 703"|>;
 save["wavepacket",report];Export[path["results/wavepacket.wxf"],Join[<|"t"->t|>,saved],"WXF"];report];
