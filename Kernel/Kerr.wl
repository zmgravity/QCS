csvRows[relative_]:=Module[{a=Import[path[relative],"CSV"]},AssociationThread[First[a],#]&/@Rest[a]];
modeTable[ell_,n_]:=modeTable[ell,n]=SortBy[Import[path["source/external_inputs/KerrQNMEFs-2/l"<>ToString[ell]<>"/s-2l"<>ToString[ell]<>"m"<>ToString[ell]<>"n"<>ToString[n]<>".dat"],"CSV"],First];
mode[ell_,n_,spin_]:=Module[{d=modeTable[ell,n],w,b},
 w=Interpolation[Transpose[{d[[All,1]],(d[[All,2]]+I d[[All,3]])/2}],InterpolationOrder->1];
 b=Interpolation[Transpose[{d[[All,1]],(d[[All,8]]+I d[[All,9]])/16}],InterpolationOrder->1];
 If[ListQ[spin],{w/@spin,b/@spin},{w[spin],b[spin]}]];
rankMetrics[m_]:=Module[{singular=SingularValueList[m,Tolerance->0],x,y,z},
 x=m[[1,1]]m[[2,2]];y=m[[1,2]]m[[2,1]];z=x/y;
 <|"s2_over_s1"->singular[[2]]/singular[[1]],"lambda"->Log[Abs[z]],"phi"->Arg[z],
 "closure_radius"->Sqrt[Log[Abs[z]]^2+Arg[z]^2],"bounded_closure"->2 Abs[x-y]/(Abs[x]+Abs[y])|>];
referenceMatrices[]:=Module[{rows=csvRows["source/output/schwarzschild_direct_exchange_rank_audit/normalized_matrices.csv"],result=<||>,m},
 Do[m=zero[2,2];Do[m[[r["alpha"],r["beta"]]]=r["normalized_real"]+I r["normalized_imag"],{r,Select[rows,#["matrix"]==name&]}];
 AssociateTo[result,name->m],{name,DeleteDuplicates[Lookup[rows,"matrix"]]}];result];
QCSKerr[]:=Module[{maxSpin,spins,kappa,w2,w4,w6,b4,b6,dp,dc,da,values,fits,at,ww41,ww61,report,rows},
 maxSpin=Min[Table[Last[modeTable[ell,n]][[1]],{ell,{2,4,6}},{n,0,7}]];
 spins=Select[modeTable[2,0][[All,1]],.9999<=#<=maxSpin&];
 kappa=Sqrt[1-spins^2]/(2(1+Sqrt[1-spins^2]));
 w2=First[mode[2,0,spins]];{w4,b4}=mode[4,0,spins];{w6,b6}=mode[6,0,spins];
 dp=2 w2-w4;dc=w4+w2-w6;da=3 w2-w6;
 values=<|"Delta_p"->dp,"Delta_c"->dc,"Delta_A"->da,"B_440"->b4,"B_660"->b6,"pole_factor"->b4 b6/(dp da)|>;
 fits=AssociationMap[Last[LeastSquares[Transpose[{ConstantArray[1.,Length[kappa]],Log[kappa]}],Log[Abs[values[#]]]]]&,Keys[values]];
 w2=First[mode[2,0,.999999]];{w4,b4}=mode[4,0,.999999];{w6,b6}=mode[6,0,.999999];
 ww41=First[mode[4,1,.999999]];ww61=First[mode[6,1,.999999]];dp=2 w2-w4;dc=w4+w2-w6;da=3 w2-w6;
 report=<|"fit_points"->Length[spins],"exponents_chi_ge_0_9999"->fits,"at_chi_0_999999"-><|"abs_K0"->Abs[b4 b6/(dp da)],
 "production_isolation_ratio"->Abs[dp]/Abs[ww41-w4],"conversion_isolation_ratio"->Abs[dc]/Abs[ww61-w6],
 "total_isolation_ratio"->Abs[da]/Abs[ww61-w6]|>|>;
 save["kerr",report];Export[path["results/kerr_scalings.csv"],Prepend[Transpose[Join[{spins,kappa},Abs/@Values[values]]],Join[{"chi","kappa"},Keys[values]]],"CSV"];
 writeCSV["kerr_target_spin_diagnostics",kerrTargets[]];report];
kerrTargets[]:=Module[{spins={0.,.2,.4,.6,.68,.8,.9,.95,.97,.98,.99,.995,.999,.9995,.9999,.99999,.999999},w2,w4,b4,w6,b6,t4,t6,dp,dc,da,d4,d6},
 Table[w2=First[mode[2,0,a]];{w4,b4}=Transpose[Table[mode[4,n,a],{n,0,7}]];{w6,b6}=Transpose[Table[mode[6,n,a],{n,0,7}]];
 t4=b4/(2 w2-w4);t6=b6/(3 w2-w6);dp=2 w2-w4[[1]];dc=w4[[1]]+w2-w6[[1]];da=3 w2-w6[[1]];
 d4=Abs[2 w2-w4];d6=Abs[3 w2-w6];
 <|"a"->a,"M_kappa"->Sqrt[1-a^2]/(2(1+Sqrt[1-a^2])),"Delta_p_abs"->Abs[dp],"Delta_c_abs"->Abs[dc],"Delta_A_abs"->Abs[da],
 "B440_abs"->Abs[b4[[1]]],"B660_abs"->Abs[b6[[1]]],"pole_factor_abs"->Abs[t4[[1]]t6[[1]]],"pole_factor_phase"->Arg[t4[[1]]t6[[1]]],
 "spacing_440_441_abs"->Abs[w4[[2]]-w4[[1]]],"spacing_660_661_abs"->Abs[w6[[2]]-w6[[1]]],
 "Delta_p_over_440_spacing"->Abs[dp/(w4[[2]]-w4[[1]])],"Delta_c_over_660_spacing"->Abs[dc/(w6[[2]]-w6[[1]])],"Delta_A_over_660_spacing"->Abs[da/(w6[[2]]-w6[[1]])],
 "Delta_p_over_440_linewidth"->Abs[dp/Im[w4[[1]]]],"Delta_A_over_660_linewidth"->Abs[da/Im[w6[[1]]]],
 "nearest_44_overtone_over_fund_distance"->Min[Rest[d4]]/First[d4],"nearest_66_overtone_over_fund_distance"->Min[Rest[d6]]/First[d6],
 "n44_within_twice_fund_distance"->Count[d4,x_/;x<=2 First[d4]],"n66_within_twice_fund_distance"->Count[d6,x_/;x<=2 First[d6]],
 "term441_over_440_abs"->Abs[t4[[2]]/t4[[1]]],"term661_over_660_abs"->Abs[t6[[2]]/t6[[1]]]|>,{a,spins}]];
