QCSPreparePair[]:=Module[{p=QCSProblem[2,64],a,x,c,t,r,initial,tout,y,outputs=<||>,saved},
 a=Last[json["source/reduced_input/pair_results.json"]]["hidden_amplitude"];
 x=Transpose[{.5 baseline[p["z"]]+a hidden[p["z"]],.5 baseline[p["z"]]-a hidden[p["z"]]}];
 c=boundary[p];Print["Propagating boundary radiation..."];
 {t,r}=propagate[p["m4"],Join[c,c.p["m4"]],x,800,.0125,.00125];
 saved=<|"t"->t,"radiation_plus"->r[[All,1;;2,1]],"radiation_minus"->r[[All,1;;2,2]]|>;
 Do[Print["Propagating cubic channel ",ell+4,"..."];p=QCSProblem[ell,64];
 initial=ConstantArray[0.+0.I,{258,2}];initial[[1;;128]]=x;
 {tout,y}=propagate[p["T"],{p["obs"],p["obs"].p["T"]},initial,70,.0125,.00125];
 AssociateTo[outputs,ell->y],{ell,{2,3}}];
 AssociateTo[saved,{"tout"->tout,"strain_plus"->Transpose[Table[outputs[ell][[All,1,1]],{ell,{2,3}}]],
 "strain_minus"->Transpose[Table[outputs[ell][[All,1,2]],{ell,{2,3}}]],
 "news_plus"->Transpose[Table[outputs[ell][[All,2,1]],{ell,{2,3}}]],
 "news_minus"->Transpose[Table[outputs[ell][[All,2,2]],{ell,{2,3}}]]}];
 Export[path["results/paired_preparations.wxf"],saved,"WXF"];saved];

chooseFilter[a_,b_,ar_,br_]:=Module[{u,ss,v,s,bv,va,delta,best=Infinity,chosen,record,f,bias,gain,objective},
 {u,ss,v}=SingularValueDecomposition[pack[a],Min[Dimensions[a]],Tolerance->0];s=Diagonal[ss];
 bv=b.v;va=ConjugateTranspose[u].ar;delta=10^-6 opnorm[ar];
 Do[f=If[kind=="tsvd",Map[If[#>cutoff First[s],1/#,0.]&,s],s/(s^2+(cutoff First[s])^2)];
 bias=opnorm[(bv.DiagonalMatrix[f]).va-br];gain=opnorm[bv.DiagonalMatrix[f]];objective=bias+delta gain;
 If[objective<best,best=objective;chosen=f;record=<|"kind"->kind,"cutoff"->cutoff,"objective"->objective,"bias"->bias,"noise_gain"->gain|>],
 {kind,{"tsvd","tikhonov"}},{cutoff,N[10^Subdivide[-9,-2,28]]}];{record,pack[v.DiagonalMatrix[chosen]],u}];

QCSReconstructPair[pair_:Automatic,regenerate_:False]:=Module[{fine,k=QCSKernels[48,regenerate],kr=QCSKernels[64,regenerate],order=26,end=40,s,sr,a,ar,b,br,best,vf,u,
 count,ti,w,r,estimated,tf,wf,p,f,pred={},truth={},saved,filters},
 fine=If[pair===Automatic,QCSData["source/reduced_input/paired_preparations_fine.mat"],pair];
 s=stateBasis[k["z"],order];sr=stateBasis[kr["z"],order];
 a=obsMatrix[k,s,{1,2},end];ar=obsMatrix[kr,sr,{1,2},end];b=outMatrix[k,s,0,70];br=outMatrix[kr,sr,0,70];
 Print["Selecting the reconstruction filter from operators..."];{best,vf,u}=chooseFilter[a,b,ar,br];
 count=Round[end/.1]+1;ti=N[Range[0,count-1]/10];w=Last[weights[ti]];
 r=Transpose[Table[Flatten[scale[w,fine["radiation_"<>sign][[1;;8(count-1)+1;;8]]]],{sign,{"plus","minus"}}]];
 estimated=vf.(ConjugateTranspose[u].r);tf=N[Range[0,2800]/40];wf=Last[weights[tf]];
 Do[p=QCSProblem[ell,48];f=response[p["T"],p["obs"],tf,.00125][[All,1,1;;96]].s;
 AppendTo[pred,f.estimated];AppendTo[truth,Transpose[Table[fine["strain_"<>sign][[1;;-1;;2,ell-1]],{sign,{"plus","minus"}}]]],{ell,{2,3}}];
 saved=<|"t"->tf,"prediction_plus"->Transpose[pred[[All,All,1]]],"prediction_minus"->Transpose[pred[[All,All,2]]],
 "reference_plus"->Transpose[truth[[All,All,1]]],"reference_minus"->Transpose[truth[[All,All,2]]],"filter"->best|>;
 Export[path["results/recovery_waveforms.wxf"],saved,"WXF"];save["reconstruction_filter",best];saved];
