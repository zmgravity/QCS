wlcp[p_]:=cp[N[ToExpression/@p]];
contour[t_,obs_,initial_,s_,radius_:.004,count_:48]:=Module[{eye=N[IdentityMatrix[Length[t]]]},
 Total[Table[With[{delta=radius Exp[2 Pi I(k+.5)/count]},
 (obs.LinearSolve[(s+delta)eye-t,initial])delta/count],{k,0,count-1}]]];
production[row_,omega_,omegaS_,conversion_]:=Module[
 {z=row["sigma"],n,d,at0,l1,l2,eye,m,ss,u,sv,v,right,left,b,parent,norm,sc=-2 I omega,
 residue,forced,g,bp,shifted,bad,ctr=0.I,ctrp=0.I,field=0.,omit=0.,delta,s,x,xp,reconstruction,factor=Sqrt[90.]/6},
 n=Length[z];{d,at0}=cheb[z];{l1,l2}=operators[z,4,d];eye=N[IdentityMatrix[n]];
 m=ArrayFlatten[{{zero[n,n],eye},{l1,l2}}];ss=First[MinimalBy[Eigenvalues[m],Abs[#+I omegaS]&]];
 {u,sv,v}=SingularValueDecomposition[pack[ss^2 eye-ss l2-l1],Tolerance->0];
 right=v[[All,-1]];right/=at0.right;left=Conjugate[u[[All,-1]]];norm=left.((2 ss eye-l2).right);
 b=ca[row["source"]];parent=ca[row["z_parent"]];residue=left.b/((ss-sc)norm);
 forced=LinearSolve[sc^2 eye-sc l2-l1,b];g=.2 parent^2;bp=b-(sc^2 eye-sc l2-l1).g;
 shifted=left.(bp/(ss-sc)-(ss eye-l2).g-sc g)/norm;bad=left.bp/((ss-sc)norm);
 Do[delta=.008 Exp[2 Pi I(k+.5)/48];s=ss+delta;
 x=LinearSolve[s^2 eye-s l2-l1,b/(s-sc)];xp=LinearSolve[s^2 eye-s l2-l1,bp/(s-sc)-(s eye-l2).g-sc g];
 ctr+=(at0.x)delta/48;ctrp+=(at0.(xp+g/(s-sc)))delta/48;
 field=Max[field,Norm[x-xp-g/(s-sc)]/Norm[x]];omit=Max[omit,Abs[at0.(x-xp)]/Abs[at0.x]],{k,0,47}];
 reconstruction=Sqrt[7.](-279 I omega+285 omega^2-4 I omega^3)/(9072 Sqrt[Pi]);
 <|"n"->n,"omega440"->I ss,"omega_error"->Abs[I ss-omegaS],"production_strain"->factor residue,
 "cubic_shifted_strain"->factor residue conversion,"forced_quadratic_strain"->factor(at0.forced)+reconstruction,
 "residue_contour_error"->Abs[ctr-residue]/Abs[residue],"field_change_residue_error"->Abs[shifted-residue]/Abs[residue],
 "field_change_contour_error"->Abs[ctrp-ctr]/Abs[ctr],"full_second_order_field_match"->field,
 "incorrect_zero_initial_data_error"->Abs[bad-residue]/Abs[residue],"omitted_output_full_response_error"->omit|>];
QCSGR[]:=Module[{prod,old,conv,omega,omegaS,rows,cascade={},p,n,particular,initial,ss,cc,direct,u,sv,v,
 right,left,norm,prediction,source,channel,sample,report},
 prod=json["source/gr_exchange/theory_direction1/schwarzschild_samples.json"];
 old=json["source/gr_exchange/investigation/overlap_inputs.json"];conv=json["source/gr_exchange/investigation/conversion_results.json"];
 omega=cp[prod["omega220"]];omegaS=wlcp[old["modes"]["440"]["omega"]];ss=-I omegaS;
 cc=cp[SelectFirst[conv,#["parent_l"]==2&&#["intermediate"]=="440"&&#["resolution"]==320&]["strain"]];
 rows=production[#,omega,omegaS,cc]&/@prod["rows"];
 Do[Print["GR cascade: ell=",ell,", n=",nn];p=QCSProblem[ell,nn];n=nn;
 source=ca[SelectFirst[prod["rows"],#["n"]==n&]["source"]];
 particular=LinearSolve[(-2 I omega)^2 IdentityMatrix[n]+2 I omega p["l2"]-p["l1"],source];
 initial=ConstantArray[0.+0.I,{4 n+2,2}];
 initial[[1;;n,1]]=-particular;initial[[1;;n,2]]=-particular+.25 p["z"]^2(1-p["z"])^2;
 initial[[n+1;;2 n,1]]=2 I omega particular;initial[[n+1;;2 n,2]]=2 I omega particular;
 direct=contour[p["T"],p["obs"],initial,ss+p["sd"]];
 {u,sv,v}=SingularValueDecomposition[pack[ss^2 IdentityMatrix[n]-ss p["l2"]-p["l1"]],Tolerance->0];
 right=v[[All,-1]];right/=p["at0"].right;left=Conjugate[u[[All,-1]]];norm=left.((2 ss IdentityMatrix[n]-p["l2"]).right);
 cc=cp[SelectFirst[conv,#["parent_l"]==ell&&#["intermediate"]=="440"&&#["resolution"]==320&]["strain"]];
 prediction=(left.((ss IdentityMatrix[n]-p["l2"]).initial[[1;;n]]+initial[[n+1;;2 n]]))/norm Sqrt[90.]/6 cc;
 AppendTo[cascade,<|"spectator_l"->ell,"n"->n,"direct"->direct,"prediction"->prediction,
 "relative_error"->Abs[direct-prediction]/Abs[prediction]|>],{ell,{2,3}},{nn,{32,48,64}}];
 report=<|"production"->rows,"cascade"->cascade,"maximum_factorization_error"->Max[Flatten[Lookup[cascade,"relative_error"]]]|>;
 save["gr",report];report];
