geometry[x_]:=Module[{rhs=x/2-1,y,w},y=Map[If[#>1,Log[Max[#,1]],#]&,Flatten[{rhs}]];
 Do[w=Exp[y];y-=(y+w-Flatten[{rhs}])/(1+w),{15}];w=Exp[y];
 If[ListQ[x],{2(1+w),w/(1+w),w},{2(1+w[[1]]),w[[1]]/(1+w[[1]]),w[[1]]}]];
potentialX[x_]:=With[{g=geometry[x],lam=9.},g[[2]](2 lam^2(lam+1)g[[1]]^3+6 lam^2 g[[1]]^2+18 lam g[[1]]+18)/(g[[1]]^3(lam g[[1]]+3)^2)];
outgoingSeries[freq_,order_]:=Module[{lam=9.,u,vc,aa,n,q,total,j},
 vc=N[CoefficientList[Normal[Series[u^2(1-2 u)(2 lam^2(lam+1)+6 lam^2 u+18 lam u^2+18 u^3)/(lam+3 u)^2,{u,0,order+1}]],u]];
 aa=ConstantArray[0.I,{order+1,Length[freq]}];aa[[1]]=ConstantArray[1.I^0,Length[freq]];
 Do[q=n+1;total=ConstantArray[0.I,Length[freq]];
 Do[j=q-term[[1]]+2;If[0<=j<n,total+=term[[2]]j(j-1)aa[[j+1]]],{term,{{4,1},{5,-4},{6,4}}}];
 Do[j=q-term[[1]]+1;If[0<=j<n,total+=term[[2]]j aa[[j+1]]],{term,{{3,2+4 I freq},{4,-10},{5,12}}}];
 Do[total-=vc[[power+1]]aa[[q-power+1]],{power,2,q}];aa[[n+1]]=total/(2 I freq n),{n,1,order}];aa];
profileRK=Compile[{{freq,_Complex,1},{h,_Real},{vv,_Real,1},{initialP,_Complex,1},{initialQ,_Complex,1},{sign,_Integer},{targets,_Integer,1}},
 Module[{p=initialP,q=initialQ,out=Table[0.I,{Length[targets]},{2},{Length[freq]}],k,j,next=1,count=(Length[vv]-1)/2,
 p2=initialP,q2=initialQ,p3=initialP,q3=initialQ,p4=initialP,q4=initialQ,k1=initialQ,k2=initialQ,k3=initialQ,k4=initialQ},
 For[k=0,k<=count,k++,If[next<=Length[targets],If[k==targets[[next]],out[[next,1]]=p;out[[next,2]]=q;next++]];
 If[k<count,k1=vv[[2 k+1]]p-2 I sign freq q;p2=p+h q/2;q2=q+h k1/2;
 k2=vv[[2 k+2]]p2-2 I sign freq q2;p3=p+h q2/2;q3=q+h k2/2;
 k3=vv[[2 k+2]]p3-2 I sign freq q3;p4=p+h q3;q4=q+h k3;
 k4=vv[[2 k+3]]p4-2 I sign freq q4;p+=h(q+2 q2+2 q3+q4)/6;q+=h(k1+2 k2+2 k3+k4)/6]];out],RuntimeOptions->"Speed"];
jost[freq0_List,step_:.02,xright_:60.,order_:20]:=Module[{freq=pack[freq0+0.I],points={-6.,-2.,2.,6.,10.},xleft=-50.,wl,vh,h1,ini,aa,rr,ur,n,h,vv,
 targets,raw,up,dUp,vp,dVp,w,ww,green,kernel,variation},
 wl=Last[geometry[xleft]];vh=(16 9.^3+40 9.^2+36 9.+18)/(8(2 9.+3)^2);h1=vh/(.25-I freq);ini={1+h1 wl,h1 wl/(2(1+wl))};
 n=Round[(Last[points]-xleft)/step];h=(Last[points]-xleft)/n;vv=pack[potentialX[xleft+Range[0,2 n]h/2]];targets=Round[(points-xleft)/h];
 raw=profileRK[freq,h,vv,pack[ini[[1]]],pack[ini[[2]]],-1,targets];
 up=Transpose[MapThread[Exp[-I freq #1]#2[[1]]&,{points,raw}]];
 dUp=Transpose[MapThread[Exp[-I freq #1](#2[[2]]-I freq #2[[1]])&,{points,raw}]];
 rr=First[geometry[xright]];ur=1/rr;aa=outgoingSeries[freq,order];ini={polyval[aa,ur],(-ur^2+2 ur^3)polyval[scale[Range[order],Rest[aa]],ur]};
 n=Round[(xright-First[points])/step];h=(First[points]-xright)/n;vv=pack[potentialX[xright+Range[0,2 n]h/2]];targets=Round[(Reverse[points]-xright)/h];
 raw=Reverse[profileRK[freq,h,vv,pack[ini[[1]]],pack[ini[[2]]],1,targets]];
 vp=Transpose[MapThread[Exp[I freq #1]#2[[1]]&,{points,raw}]];
 dVp=Transpose[MapThread[Exp[I freq #1](#2[[2]]+I freq #2[[1]])&,{points,raw}]];
 ww=up dVp-dUp vp;w=ww[[All,3]];variation=MapThread[Max[Abs[#1-#2]]/Abs[#2]&,{ww,w}];
 green=Table[-up[[k,Min[i,j]]]vp[[k,Max[i,j]]]/w[[k]],{k,Length[freq]},{i,5},{j,5}];
 kernel=Table[{-vp[[k]]/w[[k]],-up[[k]]/w[[k]]},{k,Length[freq]}];
 <|"G"->green,"K"->kernel,"W"->w,"u"->up,"v"->vp,"du"->dUp,"dv"->dVp,"wronskian_variation"->variation|>];
fluxReal[step_:.01,xright_:60.,order_:20]:=Module[{freq={.4,.6,.8,1.},data,g,k,spectral,horizon,infinity,flux,ev,sv,rows={}},
 data=jost[freq,step,xright,order];Do[g=data["G"][[j]];k=data["K"][[j]];spectral=(g-ConjugateTranspose[g])/(2 I);
 horizon=freq[[j]]Outer[Times,Conjugate[k[[1]]],k[[1]]];infinity=freq[[j]]Outer[Times,Conjugate[k[[2]]],k[[2]]];flux=horizon+infinity;
 sv=SingularValueList[spectral,Tolerance->0];ev=Eigenvalues[spectral];
 AppendTo[rows,<|"omega"->freq[[j]],"step"->step,"relative_error"->relative[flux,spectral],"horizon_omission_error"->relative[infinity,spectral],
 "infinity_omission_error"->relative[horizon,spectral],"rank_two_singular_ratio"->sv[[3]]/sv[[1]],"minimum_eigenvalue"->Min[Re[ev]],
 "wronskian_relative_variation"->data["wronskian_variation"][[j]]|>],{j,4}];{rows,data}];
QCSFlux[]:=Module[{rows={},data,ret,adv,omega,delta,z,difference,cut,literal,res,cutres,sv,poles={},report},
 omega=wlcp[json["source/gr_exchange/investigation/overlap_inputs.json"]["modes"]["440"]["omega"]];
 Do[Print["Flux factorization: radial step ",step];{data,ret}=fluxReal[step];rows=Join[rows,data];
 delta=.006 Exp[2 Pi I(Range[0,31]+.5)/32];z=omega+delta;ret=jost[z,step];adv=jost[-z,step];difference=ret["G"]-adv["G"];
 cut=Table[2 I z[[j]]Transpose[adv["K"][[j]]].ret["K"][[j]],{j,32}];
 literal=Table[2 I z[[j]]ConjugateTranspose[ret["K"][[j]]].ret["K"][[j]],{j,32}];
 res=Total[scale[delta,ret["G"]]]/32;cutres=Total[scale[delta,cut]]/32;sv=SingularValueList[res,Tolerance->0];
 AppendTo[poles,<|"step"->step,"analytic_cut_relative_error"->relative[cut,difference],"residue_cut_relative_error"->relative[cutres,res],
 "naive_conjugation_relative_error"->relative[literal,difference],"rank_one_singular_ratio"->sv[[2]]/sv[[1]]|>],{step,{.04,.02,.01}}];
 report=<|"real_frequency_checks"->rows,"pole_checks"->poles|>;save["flux",report];report];
interpolationRows[z_,sigma_]:=Module[{n=Length[z],w,row},w=(-1.)^Range[0,n-1]Sin[Pi(Range[0,n-1]+.5)/n];
 Table[row=w/(x-z);Join[row/Total[row],ConstantArray[0.,n]],{x,sigma}]];
scatteringCase[length_Integer,dt_,step_,xright_,cutoff_]:=Module[{points={-6.,-2.,2.,6.,10.},t,p,x,c,rad,truth,sigma,r,height,freq,ix,positive,ret,kern,factors,rf,spec,recovered,mask,errs,obs},
 t=N[Range[0,length-1]dt];p=QCSProblem[2,64];x=packetSeed[p,"wavepacket"];c=boundary[p];
 rad=Map[#.x&,response[p["m4"],c,t,.0025]];r=First[geometry[points]];sigma=2/r;height=r+4 Log[r/2]-2 Log[(r-2)/2];
 obs=interpolationRows[p["z"],sigma];truth=Map[#.x&,response[p["m4"],obs,t,.0025]];
 freq=N[2 Pi Join[Range[0,length/2-1],Range[-length/2,-1]]/(length dt)];ix=Flatten[Position[freq,f_/;0<f<=cutoff+10^-10]];positive=freq[[ix]];
 ret=jost[positive,step,xright,20];kern=Table[2 I positive[[j]]scale[Exp[-I positive[[j]]height],Transpose[Conjugate[ret["K"][[j]]]]],{j,Length[ix]}];
 factors=ConstantArray[0.I,{length,5,2}];Do[factors[[ix[[j]]]]=kern[[j]].DiagonalMatrix[{Exp[4 I positive[[j]]],1.}];factors[[length-ix[[j]]+2]]=Conjugate[factors[[ix[[j]]]]],{j,Length[ix]}];
 rf=Transpose[InverseFourier[#,FourierParameters->{1,-1}]&/@Transpose[rad]];
 spec=MapThread[Dot,{factors,rf}];recovered=Transpose[Fourier[#,FourierParameters->{1,-1}]&/@Transpose[spec]];
 mask=Flatten[Position[t,v_/;.5<=v<=70.]];errs=Table[relative[recovered[[mask,j]],truth[[mask,j]]],{j,5}];
 <|"length"->length,"sample_dt"->dt,"step"->step,"xright"->xright,"cutoff"->cutoff,"window_errors"->errs,"maximum_error"->Max[errs]|>];
QCSScattering[cases_:Automatic]:=Module[{configs,rows},configs=If[cases===Automatic,
 {{1024,.2,.04,1200.,4.},{1024,.2,.02,1200.,4.},{2048,.2,.02,1200.,4.},{2048,.1,.02,1200.,8.},
 {4096,.05,.01,1200.,16.},{4096,.05,.01,1600.,16.},{8192,.025,.005,1200.,32.},{8192,.05,.01,1600.,16.}},cases];
 rows=Table[Print["Continuum scattering: ",config];scatteringCase@@config,{config,configs}];save["scattering",<|"cases"->rows|>];rows];
