integrateAt[freq_,start_,targets0_,initial_,sign_,step_]:=Module[{targets=Sort[DeleteDuplicates[targets0]],x=start,p=initial,result=<||>,n,h,vv},
 If[Last[targets0]<start,targets=Reverse[targets]];
 Do[n=Max[1,Ceiling[Abs[target-x]/step]];h=(target-x)/n;vv=pack[potentialX[x+Range[0,2 n]h/2]];
 p=profileRK[{N[freq+0.I]},h,vv,{N[p[[1]]+0.I]},{N[p[[2]]+0.I]},sign,{n}][[1,All,1]];
 AssociateTo[result,target->p];x=target,{target,targets}];result];
seriesValues[freq_,x_,order_]:=Module[{r,f,w,a},{r,f,w}=geometry[x];a=outgoingSeries[{freq},order][[All,1]];
 {polyval[a,1/r],-f/r^2 polyval[Range[order]Rest[a],1/r]}];
advancedFactors[freq_,sigma_,step_:.00125,xright_:60.,order_:20]:=Module[{r=2/sigma,x,height,core,wl,h1,ur,vr,pu,qu,pv,qv,wa,
 po,qo,pi,qi,a,b,hu,hv,dhu,dhv,up,uq,vp,vq,fac,eu,ev,dx},
 x=r+2 Log[r/2-1];height=r+4 Log[r/2]-2 Log[(r-2)/2];core=Sort[DeleteDuplicates[Join[{0.,xright},Select[x,#<=xright&]]]];
 wl=Last[geometry[-50.]];h1=((16 9.^3+40 9.^2+36 9.+18)/(8(2 9.+3)^2))/(.25+I freq);
 ur=integrateAt[freq,-50.,core,{1+h1 wl,h1 wl/(2(1+wl))},1,step];
 vr=integrateAt[freq,xright,Reverse[core],seriesValues[-freq,xright,order],-1,step];
 {pu,qu}=ur[0.];{pv,qv}=vr[0.];wa=pu qv-qu pv-2 I freq pu pv;
 {pu,qu}=ur[xright];{po,qo}=seriesValues[freq,xright,order];{pi,qi}=seriesValues[-freq,xright,order];
 {a,b}=LinearSolve[{{po,pi},{qo+I freq po,qi-I freq pi}},{pu,qu+I freq pu}];
 hu=hv=dhu=dhv=ConstantArray[0.I,Length[x]];
 Do[dx=-2/(sigma[[j]]^2(1-sigma[[j]]));
 If[x[[j]]<=xright,{up,uq}=ur[x[[j]]];{vp,vq}=vr[x[[j]]],
 {po,qo}=seriesValues[freq,x[[j]],order];{vp,vq}=seriesValues[-freq,x[[j]],order];fac=b Exp[-2 I freq(x[[j]]-xright)];up=a po+fac vp;uq=a qo+fac(vq-2 I freq vp)];
 eu=Exp[I freq(x[[j]]-height[[j]])];ev=Exp[-I freq(x[[j]]+height[[j]])];hu[[j]]=eu up;hv[[j]]=ev vp;
 dhu[[j]]=dx eu(uq+2 I freq sigma[[j]]^2 up);dhv[[j]]=dx ev(vq-2 I freq(1-sigma[[j]]^2)vp),{j,Length[x]}];
 <|"k"->-{hv,hu}/wa,"dk"->-{dhv,dhu}/wa,"at_infinity"->{0.,-a/wa}|>];
productionResidues[n_,radius_:.004,count_:48,fieldChange_:False,wrongData_:False]:=Module[{prod,row,p,z,d,at0,at1,l1,l2,eye,seeds,b,g,omega,ss,sp,delta,s,rhs,profile,field,amplitudes},
 prod=json["source/gr_exchange/theory_direction1/schwarzschild_samples.json"];row=SelectFirst[prod["rows"],#["n"]==n&];p=QCSProblem[2,n];z=p["z"];d=p["d"];at0=p["at0"];at1=boundary[p][[1,1;;n]];
 l1=p["l1"];l2=p["l2"];eye=IdentityMatrix[n];b=ca[row["source"]];g=.2 ca[row["z_parent"]]^2;
 omega=wlcp[json["source/gr_exchange/investigation/overlap_inputs.json"]["modes"]["440"]["omega"]];ss=-I omega;sp=-2 I cp[prod["omega220"]];
 seeds=Transpose[{ConstantArray[0.,n],.25 z^2(1-z)^2}];field=ConstantArray[0.I,{n,2}];amplitudes=ConstantArray[0.I,{2,2}];
 Do[delta=radius Exp[2 Pi I(j+.5)/count];s=ss+delta;rhs=Transpose[{b,b}]/(s-sp)+(s eye-l2).seeds;
 If[fieldChange,rhs-=Transpose[ConstantArray[(sp^2 eye-sp l2-l1).g/(s-sp),2]];
 If[!wrongData,rhs-=Transpose[ConstantArray[(s eye-l2).g+sp g,2]]]];
 profile=LinearSolve[s^2 eye-s l2-l1,rhs];If[fieldChange,profile+=Transpose[{g,g}]/(s-sp)];
 field+=delta profile/count;amplitudes+=delta {Exp[-4 s](at1.profile),at0.profile}/count,{j,0,count-1}];
 <|"z"->z,"field"->field,"derivative"->d.field,"amplitudes"->amplitudes,"omega"->omega|>];
convertProfile[channel_,sample_,profile_,derivative_,atInfinity_,ss_]:=Module[{z=sample["sigma"],n,d,at0,o1,o2,sd,star,aq,bq,powers,source,output,fz,fq,rp,rho},
 n=Length[z];{d,at0}=cheb[z];{o1,o2}=operators[z,channel["spectator_l"]+4,d];sd=-I cp[channel["omega"]];star=ss+sd;
 aq=ca[sample["time_qcoeff"]];bq=ca[sample["time_dqcoeff"]];powers=ss^Range[0,Length[First[aq]]-1];source=scale[aq.powers,profile]+scale[bq.powers,derivative];
 output=LinearSolve[star^2 IdentityMatrix[n]-star o2-o1,source];fz=alpha[channel["spectator_l"]+4]/(alpha[2]^2 alpha[channel["spectator_l"]]);fq=alpha[4]/alpha[2]^2;
 rp=ca[channel["output_polynomial"]];rho=rp.(ss^Range[0,Length[rp]-1])+cp[channel["output_primitive"]]/ss+cp[channel["output_product_primitive"]]/star;
 fz(at0.output)+fq rho atInfinity];
poleSewCase[n_,step_,xright_:60.,order_:20,radius_:.004,count_:48,fieldChange_:False,frequencyDependent_:False]:=Module[
 {prod,omega,ss,adv,sewn,dsewn,pinf,channels={},ch,sample,value,p,row,source,sp,particular,initial,reference,
 z,eye,seeds,b,at0,at1,delta,s,freq,profile,amp,retarded},
 prod=productionResidues[n,radius,count,fieldChange];omega=prod["omega"];ss=-I omega;z=prod["z"];p=QCSProblem[2,n];
 If[!frequencyDependent,adv=advancedFactors[omega,z,step,xright,order];
 sewn=2 I omega Transpose[adv["k"]].prod["amplitudes"];dsewn=2 I omega Transpose[adv["dk"]].prod["amplitudes"];pinf=2 I omega adv["at_infinity"].prod["amplitudes"],
 sewn=dsewn=retarded=ConstantArray[0.I,{n,2}];pinf=ConstantArray[0.I,2];eye=IdentityMatrix[n];at0=p["at0"];at1=boundary[p][[1,1;;n]];
 row=SelectFirst[json["source/gr_exchange/theory_direction1/schwarzschild_samples.json"]["rows"],#["n"]==n&];b=ca[row["source"]];
 sp=-2 I cp[json["source/gr_exchange/theory_direction1/schwarzschild_samples.json"]["omega220"]];seeds=Transpose[{ConstantArray[0.,n],.25 z^2(1-z)^2}];
 Do[delta=radius Exp[2 Pi I(j+.5)/count];s=ss+delta;freq=I s;profile=LinearSolve[s^2 eye-s p["l2"]-p["l1"],Transpose[{b,b}]/(s-sp)+(s eye-p["l2"]).seeds];
 amp={Exp[-4 s](at1.profile),at0.profile};adv=advancedFactors[freq,z,step,xright,order];
 sewn+=delta/count 2 I freq Transpose[adv["k"]].amp;dsewn+=delta/count 2 I freq Transpose[adv["dk"]].amp;pinf+=delta/count 2 I freq adv["at_infinity"].amp,{j,0,count-1}]];
 Do[sample=SelectFirst[ch["samples"],#["n"]==n&];value=convertProfile[ch,sample,sewn,dsewn,pinf,ss];p=QCSProblem[ch["spectator_l"],n];
 row=SelectFirst[json["source/gr_exchange/theory_direction1/schwarzschild_samples.json"]["rows"],#["n"]==n&];source=ca[row["source"]];sp=-2 I cp[json["source/gr_exchange/theory_direction1/schwarzschild_samples.json"]["omega220"]];
 particular=LinearSolve[sp^2 IdentityMatrix[n]-sp p["l2"]-p["l1"],source];initial=ConstantArray[0.I,{4 n+2,2}];
 initial[[1;;n,1]]=-particular;initial[[1;;n,2]]=-particular+.25 z^2(1-z)^2;initial[[n+1;;2 n,1]]=-sp particular;initial[[n+1;;2 n,2]]=-sp particular;
 reference=Quiet[contour[p["T"],p["obs"],initial,ss+p["sd"],radius,count],LinearSolve::luc];
 AppendTo[channels,<|"output"->ch["spectator_l"]+4,"sewn"->value,"reference"->reference,"relative_error"->Abs[value/reference-1]|>],
 {ch,json["source/gr_exchange/theory_direction1/cascade_evolution_inputs.json"]["rows"]}];
 <|"n"->n,"step"->step,"xright"->xright,"order"->order,"radius"->radius,"count"->count,"field_change"->fieldChange,
 "frequency_dependent"->frequencyDependent,"profile_relative_error"->relative[sewn,prod["field"]],"channels"->channels|>];
QCSPoleScattering[full_:True]:=Module[{cases,rows,contours={},maximum,report,bad,base},
 cases=If[full,{{32,.00125},{48,.00125},{64,.00125},{48,.000625},{64,.000625},{48,.00125,50.},{48,.00125,80.},
 {48,.00125,60.,14},{48,.00125,60.,26},{48,.00125,60.,20,.002,32},{48,.00125,60.,20,.008,64},{48,.00125,60.,20,.004,48,True}},{{48,.00125}}];
 rows=Table[Print["Pole scattering: ",config];poleSewCase@@config,{config,cases}];
 If[full,contours=Table[Print["Frequency-dependent pole contour: ",config];poleSewCase@@config,
 {config,{{32,.005,60.,20,.004,32,False,True},{32,.0025,60.,20,.004,32,False,True},{32,.00125,60.,20,.008,32,False,True},{48,.00125,60.,20,.004,32,False,True}}}]];
 maximum=Max[Flatten[Lookup[Flatten[Lookup[rows,"channels"]],"relative_error"]]];
 base=productionResidues[48];bad=productionResidues[48,.004,48,True,True];
 report=<|"refinements"->rows,"contour_sewing"->contours,"maximum_sewn_relative_error"->maximum,"original_uniform_target"->10^-7,
 "original_uniform_target_passed"->(maximum<10^-7),"incorrect_reset_data_amplitude_errors"->Abs[bad["amplitudes"]/base["amplitudes"]-1]|>;
 save["pole_scattering",report];report];
