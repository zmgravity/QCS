BeginPackage["QCS`"];
QCSRoot::usage="QCSRoot is the directory containing this computational release.";
QCSData::usage="QCSData[relative] imports a retained JSON or MAT data file.";
QCSProblem::usage="QCSProblem[ell,n] constructs the Schwarzschild linear and quadratic cascade operators.";
QCSKernels::usage="QCSKernels[n,regenerate] loads or regenerates boundary and cubic response operators.";
QCSPreparePair::usage="QCSPreparePair[] propagates the retained pair natively in Wolfram Language.";
QCSReconstructPair::usage="QCSReconstructPair[pair,regenerate] reconstructs cubic waveforms from boundary radiation. True regenerates its response operators.";
QCSCheckData::usage="QCSCheckData[] recomputes energies and figure diagnostics from the retained arrays.";
QCSFigures::usage="QCSFigures[] makes numerical figures using the retained reference arrays.";
QCSKerr::usage="QCSKerr[] calculates Kerr pole factors and scaling exponents from the public mode tables.";
QCSGR::usage="QCSGR[] independently evaluates the Schwarzschild production and cascade residues.";
QCSTwoPole::usage="QCSTwoPole[] rebuilds the two-pole exchange matrix with native radial integration.";
QCSCausal::usage="QCSCausal[] performs native causal waveform fits and noise forecasts.";
QCSWavepacket::usage="QCSWavepacket[] runs forward propagation and boundary reconstruction for smooth packets.";
QCSRWZ::usage="QCSRWZ[] recomputes the quadratic RWZ benchmark from the retained Wolfram source expressions.";
QCSFlux::usage="QCSFlux[] verifies the classical two-boundary spectral identity using native Jost integration.";
QCSScattering::usage="QCSScattering[] performs the continuum scattering synthesis and refinement cases.";
QCSPoleScattering::usage="QCSPoleScattering[] checks residue sewing and the frequency-dependent radiation contour.";
QCSRun::usage="QCSRun[name] runs a documented Wolfram workflow and writes results under results/.";
Begin["`Private`"];
QCSRoot=DirectoryName[DirectoryName[$InputFileName]];
If[!DirectoryQ[FileNameJoin[{QCSRoot,"results"}]],CreateDirectory[FileNameJoin[{QCSRoot,"results"}]]];
path[s_String]:=FileNameJoin[Prepend[StringSplit[s,"/"],QCSRoot]];
json[s_String]:=Import[path[s],"RawJSON"];
pack[a_]:=Developer`ToPackedArray[N[a]];
cp[a_List]:=a[[1]]+I a[[2]];
ca[a_List]:=Map[cp,a,{ArrayDepth[a]-1}];
fro[a_]:=Norm[Flatten[a]];
opnorm[a_]:=First[SingularValueList[pack[a],1]];
scale[v_,m_]:=MapThread[Times,{v,m}];
zero[n_,m_]:=ConstantArray[0.,{n,m}];
toJSON[a_]:=a/.Complex[x_,y_]:>{x,y};
save[name_String,value_]:=(If[!DirectoryQ[path["results"]],CreateDirectory[path["results"]]];
 Export[path["results/"<>name<>".json"],toJSON[value],"RawJSON"]);
QCSData[relative_String]:=If[StringEndsQ[relative,".mat"],
 Module[{f=path[relative],keys,data,shapes},keys=Import[f,"Labels"];data=Import[f,"Data"];
 shapes=json[StringDrop[relative,-4]<>".shapes.json"];
 AssociationThread[keys,MapThread[If[shapes[#1]==={},First[Flatten[#2]],
    pack[ArrayReshape[Flatten[#2],shapes[#1]]]]&,{keys,data}]]],json[relative]];

cheb[z_]:=Module[{n=Length[z],w,d,at0},
 w=(-1.)^Range[0,n-1] Sin[Pi (Range[0,n-1]+.5)/n];
 d=Table[If[i==j,0.,w[[j]]/w[[i]]/(z[[i]]-z[[j]])],{i,n},{j,n}];
 Do[d[[i,i]]=-Total[d[[i]]],{i,n}];at0=w/(-z);{pack[d],pack[at0/Total[at0]]}];
operators[z_,ell_,d_]:=Module[{lam=(ell-1)(ell+2)/2.,v,l1,l2},
 v=(8 lam^2(lam+1)+12 lam^2 z+18 lam z^2+9 z^3)/(2 lam+3 z)^2;
 l1=scale[z^2(1-z)/(16(1+z)),d.d]+scale[z(2-3 z)/(16(1+z)),d]-DiagonalMatrix[v/(16(1+z))];
 l2=scale[(1-2 z^2)/(4(1+z)),d]-DiagonalMatrix[z/(2(1+z))];{pack[l1],pack[l2]}];
derivative[z_,coef_,d_]:=Module[{res=zero[Length[z],Length[z]],power=IdentityMatrix[Length[z]]},
 Do[res+=scale[coef[[All,k]],power];power=power.d,{k,Length[First[coef]]}];pack[res]];
alpha[ell_]:=Sqrt[(ell+2)(ell+1)ell(ell-1)]/2.;
QCSProblem[ell_Integer,n_Integer]:=QCSProblem[ell,n]=Module[
 {channel,sample,z,d,at0,l1,l2,o1,o2,m4,mo,sd,coupling,total,obs,row,fz,fq},
 channel=SelectFirst[json["source/gr_exchange/theory_direction1/cascade_evolution_inputs.json"]["rows"],#["spectator_l"]==ell&];
 sample=SelectFirst[channel["samples"],#["n"]==n&];z=pack[sample["sigma"]];{d,at0}=cheb[z];
 {l1,l2}=operators[z,4,d];{o1,o2}=operators[z,ell+4,d];
 m4=pack[ArrayFlatten[{{zero[n,n],IdentityMatrix[n]},{l1,l2}}]];
 mo=pack[ArrayFlatten[{{zero[n,n],IdentityMatrix[n]},{o1,o2}}]];
 sd=-I cp[channel["omega"]];
 coupling=Join[derivative[z,ca[sample["qcoeff"]],d],derivative[z,ca[sample["vcoeff"]],d],2];
 total=ConstantArray[0.+0.I,{4 n+2,4 n+2}];
 total[[1;;2 n,1;;2 n]]=m4+sd IdentityMatrix[2 n];
 total[[2 n+1;;4 n,2 n+1;;4 n]]=mo;total[[3 n+1;;4 n,1;;2 n]]=coupling;
 total[[4 n+1,4 n+1]]=sd;total[[4 n+1,1;;n]]=at0;total[[4 n+2,1;;n]]=at0;
 fz=alpha[ell+4]/(alpha[2]^2 alpha[ell]);fq=alpha[4]/alpha[2]^2;
 obs=ConstantArray[0.+0.I,4 n+2];obs[[2 n+1;;3 n]]=fz at0;row=Join[at0,ConstantArray[0.,n]];
 Do[obs[[1;;2 n]]+=fq coefficient row;row=row.m4,{coefficient,ca[channel["output_polynomial"]]}];
 obs[[4 n+1]]=fq cp[channel["output_primitive"]];obs[[4 n+2]]=fq cp[channel["output_product_primitive"]];
 <|"n"->n,"z"->z,"d"->d,"at0"->at0,"m4"->m4,"mo"->mo,"T"->pack[total],"obs"->pack[obs],
 "coupling"->pack[coupling],"l1"->l1,"l2"->l2,"sd"->sd|>];
boundary[p_]:=Module[{z=p["z"],n=p["n"],w,at1},
 w=(-1.)^Range[0,n-1] Sin[Pi(Range[0,n-1]+.5)/n];at1=w/(1-z);at1/=Total[at1];
 pack[{Join[at1,ConstantArray[0.,n]],Join[p["at0"],ConstantArray[0.,n]]}]];
rk[a_,dt_]:=Module[{e=N[IdentityMatrix[Length[a]]],z=dt a},pack[e+z.(e+z.(e/2+z.(e/6+z/24)))]];
response[a_,observation_,times_,dt_]:=Module[{step,row},
 step=MatrixPower[rk[a,dt],Round[(times[[2]]-times[[1]])/dt]];
 row=If[VectorQ[observation],{observation},observation];
 pack[NestList[#.step&,pack[row],Length[times]-1]]];
propagate[a_,observation_,initial_,end_,sample_,dt_]:=Module[{step,state,times,values},
 step=MatrixPower[rk[a,dt],Round[sample/dt]];state=pack[initial];times=N[Range[0,Round[end/sample]] sample];
 values=Table[With[{v=observation.state},state=step.state;v],{Length[times]}];{times,pack[values]}];

gauss[n_Integer]:=gauss[n]=Module[{b,vals,vec,order},
 b=N[Table[k/Sqrt[4 k^2-1],{k,n-1}]];
 {vals,vec}=Eigensystem[Normal[SparseArray[{Band[{1,2}]->b,Band[{2,1}]->b},{n,n}]]];
 order=Ordering[vals];{(vals[[order]]+1)/2,vec[[order,1]]^2}];
basis[z_,order_]:=Module[{x=2 z-1,e=z^2(1-z)^2,de=2 z(1-z)(1-2 z),t,dt},
 t=Transpose[Table[ChebyshevT[k,x],{k,0,order-1}]];
 dt=Transpose[Prepend[Table[2 k ChebyshevU[k-1,x],{k,1,order-1}],ConstantArray[0.,Length[z]]]];
 {pack[scale[e,t]],pack[scale[de,t]+scale[e,dt]]}];
gram[order_]:=gram[order]=Module[{z,w,q,d,lam=9.,u,gq,gv,g,wh},
 {z,w}=gauss[Max[160,4 order]];{q,d}=basis[z,order];
 u=(8 lam^2(lam+1)+12 lam^2 z+18 lam z^2+9 z^3)/(2 lam+3 z)^2;
 gq=Transpose[d].scale[w z^2(1-z)/4,d]+Transpose[q].scale[w u/4,q];
 gv=Transpose[q].scale[w 4(1+z),q];g=ArrayFlatten[{{gq,zero[order,order]},{zero[order,order],gv}}];
 g=(g+Transpose[g])/2;
 wh=LinearSolve[CholeskyDecomposition[g],N[IdentityMatrix[2 order]]];{pack[wh],pack[g]}];
stateBasis[z_,order_,wh_:Automatic]:=Module[{q=First[basis[z,order]],zz},zz=zero[Length[z],order];
 pack[ArrayFlatten[{{q,zz},{zz,q}}].If[wh===Automatic,First[gram[order]],wh]]];
energy[c_,ng_:320]:=Module[{z,w,q,d,order=Length[c]/2,lam=9.,u,qq,dq,v},
 {z,w}=gauss[ng];{q,d}=basis[z,order];qq=q.Take[c,order];dq=d.Take[c,order];v=q.Drop[c,order];
 u=(8 lam^2(lam+1)+12 lam^2 z+18 lam z^2+9 z^3)/(2 lam+3 z)^2;
 Total[w(z^2(1-z)Abs[dq]^2/4+u Abs[qq]^2/4+4(1+z)Abs[v]^2)]];
baseline[z_]:=Module[{c=ConstantArray[0.,36]},c[[1]]=.25;c/=Sqrt[energy[c,256]];stateBasis[z,18,IdentityMatrix[36]].c];
hiddenCoefficients[]:=hiddenCoefficients[]=QCSData["source/reduced_input/ambiguity_coefficients.mat"]["p18_infinity_400_strain_0_70"];
hidden[z_]:=stateBasis[z,18,IdentityMatrix[36]].hiddenCoefficients[];
weights[t_,lo_:0,hi_:Infinity]:=Module[{ix,w},ix=Flatten[Position[t,x_/;lo-10^-8<=x<=hi+10^-8]];
 w=ConstantArray[t[[2]]-t[[1]],Length[ix]];w[[1]]/=2;w[[-1]]/=2;{ix,Sqrt[w]}];
l2[y_,t_]:=Sqrt[Total[MapThread[Times,{Differences[t]/2,Most[Abs[y]^2]+Rest[Abs[y]^2]}]]];
obsMatrix[k_,s_,channels_,end_]:=Module[{ix,w,vals},{ix,w}=weights[k["t"],0,end];
 vals=Map[#.s&,k["O"][[ix,channels,All]]];pack[Flatten[scale[w,vals],1]]];
outMatrix[k_,s_,lo_,hi_,news_:False]:=Module[{ix,w},{ix,w}=weights[k["tout"],lo,hi];
 pack[Join@@Table[scale[w,k["F"<>If[news,"d",""]<>ToString[j]][[ix]].s],{j,{6,7}}]]];
QCSKernels[n_Integer,False]:=QCSKernels[n,False]=QCSData["source/reduced_input/kernels_n"<>ToString[n]<>"_end400_sample0.1.mat"];
QCSKernels[n_Integer]:=QCSKernels[n,False];
QCSKernels[n_Integer,True]:=Module[{dt=If[n==32,.005,.0025],t=N[Range[0,4000]/10],tout=N[Range[0,1000]/10],p,c,result,pp},
 p=QCSProblem[2,n];c=boundary[p];result=<|"t"->t,"z"->p["z"],"tout"->tout,
 "O"->response[p["m4"],c,t,dt],"Od"->response[p["m4"],c.p["m4"],t,dt]|>;
 Do[pp=QCSProblem[ell,n];AssociateTo[result,{
 "F"<>ToString[ell+4]->response[pp["T"],pp["obs"],tout,dt][[All,1,1;;2 n]],
 "Fd"<>ToString[ell+4]->response[pp["T"],pp["obs"].pp["T"],tout,dt][[All,1,1;;2 n]]}],{ell,{2,3}}];result];
Get[FileNameJoin[{DirectoryName[$InputFileName],"Pair.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"Diagnostics.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"GR.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"Kerr.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"TwoPole.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"Causal.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"Wavepacket.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"RWZ.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"Scattering.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"PoleScattering.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"Figures.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"Run.wl"}]];
End[];
EndPackage[];
