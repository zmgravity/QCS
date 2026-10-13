fresponse[w_,p_,t_]:=If[Abs[w-p]<10^-10,-I t Exp[-I p t],(Exp[-I w t]-Exp[-I p t])/(w-p)];
causalDesign[t_,drive_,internal_,spectator_,output_,overtone_,free_:False]:=Module[{total=drive+spectator,cols},
 cols=If[free,Join[{Exp[-I total t]},Exp[-I(# +spectator)t]&/@Take[internal,If[overtone,2,1]]],
 (fresponse[total,output,t]-fresponse[#+spectator,output,t])&/@Take[internal,If[overtone,2,1]]];
 Transpose[Join[cols,{Exp[-I output t],Exp[-I(Re[total]+.23-I/1.7)t],Exp[-.12 I t]/(1+t/4)^4}]]];
boundedBatch[m_]:=Module[{x=m[[All,1,1]]m[[All,2,2]],y=m[[All,1,2]]m[[All,2,1]]},2 Abs[x-y]/(Abs[x]+Abs[y])];
normalComplex[dims_]:=(RandomVariate[NormalDistribution[],dims]+I RandomVariate[NormalDistribution[],dims])/Sqrt[2.];
writeCSV[name_,rows_]:=Export[path["results/"<>name<>".csv"],Prepend[Values/@rows,Keys[First[rows]]],"CSV"];
QCSCausal[samples_Integer:100000,matrices_:Automatic]:=Module[
 {m=If[matrices===Automatic,referenceMatrices[],matrices],gains,a0,a1,asum,wp,wi,wo,t=N[Range[0,1400]/20],
 starts={0.,2.,5.,8.,12.},specs,exchange,waves,exNorm,base,nu,rec,res,cov,ix,design,coef,selector,
 recovery={},power={},mediator={},contact={},covariances=<||>,fits=<||>,full,null,basis441,orth,nc,threshold,
 std,nullDraw,altDraw,nullStat,altStat,prob,dirs,stat,bias,baseline,report,fit,y,free,over,tag,ratio,rhos},
 gains=Outer[Times,{1.,.73 Exp[.41 I]},{1.,1.24 Exp[-.27 I]}];a0=gains m["n0_direct_vertices"];a1=gains m["n1_direct_vertices"];asum=a0+a1;
 wp=First[mode[#,0,0.]]&/@{2,3};wi=First[mode[4,#,0.]]&/@{0,1};wo={First[mode[6,0,0.]],1.40973512-.09550963 I};
 exchange=Table[base=causalDesign[t,2 wp[[b]],wi,wp[[a]],wo[[a]],True];base[[All,1;;2]].{a0[[a,b]],a1[[a,b]]},{a,2},{b,2}];
 waves=Table[base=causalDesign[t,2 wp[[b]],wi,wp[[a]],wo[[a]],True];
 exchange[[a,b]]+Max[Abs[exchange[[a,b]]]]Exp[I(.37+.51(2(a-1)+b-1))](base[[All,3;;5]].{.55,.35,.15}),{a,2},{b,2}];
 exNorm=Map[Norm,exchange,{2}];specs={{"causal_full",False,True},{"causal_omit_441",False,False},{"free_full",True,True},{"free_omit_441",True,False}};
 Do[ix=Flatten[Position[t,x_/;x>=start]];
 Do[{tag,free,over}=spec;rec=ConstantArray[0.I,{2,2}];res=zero[2,2];cov=zero[2,2];
 Do[design=causalDesign[t[[ix]],2 wp[[b]],wi,wp[[a]],wo[[a]],over,free];y=waves[[a,b,ix]];
 coef=LeastSquares[design,y];fit=design.coef;res[[a,b]]=Norm[y-fit]/Norm[y];
 selector=ConstantArray[0.,Length[coef]];If[free,selector[[1]]=2 wp[[b]]+wp[[a]]-wo[[a]],selector[[1;;If[over,2,1]]]=1.];
 rec[[a,b]]=selector.coef;cov[[a,b]]=Re[selector.LinearSolve[ConjugateTranspose[design].design,Conjugate[selector]]];
 If[start==0.&&a==1&&b==1,AssociateTo[fits,tag->fit]],{a,2},{b,2}];
 AssociateTo[covariances,ToString[{tag,start},InputForm]->cov];
 AppendTo[recovery,Join[<|"start_time_M"->start,"temporal_model"->tag,"relative_matrix_error"->fro[rec-asum]/fro[asum],"max_waveform_residual"->Max[res]|>,rankMetrics[rec]]],{spec,specs}],{start,starts}];
 writeCSV["causal_recovery",recovery];Export[path["results/causal_waveforms.wxf"],<|"t"->t,"waveforms"->waves,"fits"->fits|>,"WXF"];
 SeedRandom[20260814,Method->"MersenneTwister"];rhos={10,15,20,30,40,50,75,100,150,200,300,500,750,1000};
 Do[Print["Causal noise forecast: ",tag,", start ",start];
 Do[std=exNorm/rho Sqrt[covariances[ToString[{tag,start},InputForm]]];
 nullDraw=(a0+std #)&/@normalComplex[{samples,2,2}];altDraw=(asum+std #)&/@normalComplex[{samples,2,2}];
 nullStat=boundedBatch[nullDraw];altStat=boundedBatch[altDraw];threshold=Quantile[nullStat,.99,{{0,0},{1,0}}];
 prob=N[Count[altStat,x_/;x>threshold]/samples];AppendTo[power,<|"temporal_model"->tag,"start_time_M"->start,"rho"->rho,"threshold"->threshold,"power"->prob,"power_mc_se"->Sqrt[prob(1-prob)/samples]|>],{rho,rhos}],{tag,{"causal_full","free_full"}},{start,starts}];
 threshold=Quantile[ChiSquareDistribution[8],.99];
 Do[ix=Flatten[Position[t,x_/;x>=start]];nc=0.;
 Do[full=causalDesign[t[[ix]],2 wp[[b]],wi,wp[[a]],wo[[a]],True];null=causalDesign[t[[ix]],2 wp[[b]],wi,wp[[a]],wo[[a]],False];
 basis441=full[[All,2]];orth=basis441-null.LeastSquares[null,basis441];nc+=2(Norm[a1[[a,b]]orth]/exNorm[[a,b]])^2,{a,2},{b,2}];
 Do[AppendTo[mediator,<|"start_time_M"->start,"rho"->rho,"chi2_threshold"->threshold,"noncentrality"->nc rho^2,
 "power"->SurvivalFunction[NoncentralChiSquareDistribution[8,nc rho^2],threshold]|>],{rho,rhos}],{start,starts}];
 SeedRandom[20260816,Method->"MersenneTwister"];dirs=normalComplex[{samples,2,2}];dirs=(#/fro[#])&/@dirs;baseline=rankMetrics[asum]["bounded_closure"];
 Do[stat=boundedBatch[(asum+fraction fro[asum]#)&/@dirs];bias=Abs[stat-baseline]/baseline;
 AppendTo[contact,<|"contact_frobenius_fraction"->fraction,"samples"->samples,"D_median"->Median[stat],"D_q05"->Quantile[stat,.05],"D_q95"->Quantile[stat,.95],
 "relative_bias_median"->Median[bias],"prob_relative_bias_gt_0p25"->N[Count[bias,x_/;x>.25]/samples],"prob_relative_bias_gt_0p50"->N[Count[bias,x_/;x>.5]/samples]|>],{fraction,N[Range[0,20]/200]}];
 writeCSV["causal_detection_power",power];writeCSV["mediator_likelihood_power",mediator];writeCSV["contact_sensitivity",contact];
 report=<|"samples"->samples,"random_generator"->"Wolfram MersenneTwister", "input_closure"->baseline,
 "full_fit_max_matrix_error"->Max[Lookup[Select[recovery,#["temporal_model"]=="causal_full"&],"relative_matrix_error"]],"recovery"->recovery|>;
 save["causal",report];report];
