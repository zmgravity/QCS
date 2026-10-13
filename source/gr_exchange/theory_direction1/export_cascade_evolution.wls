base = StringReplace[DirectoryName[DirectoryName[$InputFileName]], "\\" -> "/"] <> "/";
$MaxExtraPrecision=500;
inp=Import[base<>"investigation/overlap_inputs.json","RawJSON"];
cx[z_]:=SetPrecision[ToExpression[z[[1]]]+I ToExpression[z[[2]]],70];
enc[z_]:={Re[N[z,17]],Im[N[z,17]]};
poly[a_,x_]:=Fold[#1*x+#2&,Last[a],Reverse[Most[a]]];
JS[a_,b_,c_]:=ThreeJSymbol[a,b,c];
cc[j_]:=(j+2)*(j-1)*((j+2)*(j-1)+2)/6;
hh[j_]:=cc[j]+3*x^2*(1-x)/((j+2)*(j-1)+3*x);
az[j_,ww_]:=(hh[j]+2 I ww*(1-2*x^2))/(cc[j]+2 I ww);
bz[j_,ww_]:=-x^2*(1-x)/(cc[j]+2 I ww);
profile[key_,j_,grid_]:=Module[{ww,aa,da,dda,p,dp,ddp,z,dz},
 ww=cx[inp["modes"][key]["omega"]];aa=cx/@inp["modes"][key]["a"];
 da=Range[Length[aa]-1] Rest[aa];dda=Range[Length[da]-1]Rest[da];
 p=(poly[aa,1-#]&/@grid)/Total[aa];
 dp=-(poly[da,1-#]&/@grid)/Total[aa];
 ddp=(poly[dda,1-#]&/@grid)/Total[aa];
 z=(az[j,ww]/.x->grid)*p+(bz[j,ww]/.x->grid)*dp;
 dz=(D[az[j,ww],x]/.x->grid)*p+
     ((az[j,ww]+D[bz[j,ww],x])/.x->grid)*dp+
     (bz[j,ww]/.x->grid)*ddp;
 {z,dz}
];
lam=9;
pot=(8 lam^2(lam+1)+12 lam^2*x+18 lam*x^2+9*x^3)/(2 lam+3*x)^2;
Lone[e_]:=x^2*(1-x)/(16*(1+x))*D[e,{x,2}]+
 x*(2-3*x)/(16*(1+x))*D[e,x]-pot/(16*(1+x))*e;
Ltwo[e_]:=(1-2*x^2)/(4*(1+x))*D[e,x]-x/(2*(1+x))*e;
td={q[x],v[x],Lone[q[x]]+Ltwo[v[x]],
 Lone[v[x]]+Ltwo[Lone[q[x]]]+Ltwo[Ltwo[v[x]]]};
rows={};
Do[
 M=1;l1=j;m1=j;l2=4;m2=4;l=j+4;m=l;
 Clear[w1,w2,s];
 Get[base<>"investigation/public_inputs/SourceEEE.wl"];
 ss=SourceEEE/.{\[Psi]1l1e[r]->z1,\[Psi]1l2e[r]->z2,
 Derivative[1][\[Psi]1l1e][r]->dz1,Derivative[1][\[Psi]1l2e][r]->dz2,
 \[Omega]1->w1,\[Omega]2->I*s};
 hp=1+2*x-x/(1-x);
 ex=Together[-(ss/.{r->2/x,z1->p1,dz1->I*w1*hp*p1-x^2*dp1/2,
                          z2->q0,dz2->-s*hp*q0-x^2*q1/2})/(4*x^2*(1-x^2))];
 aq=Table[Together[Coefficient[Coefficient[ex,q0],s,k]],{k,0,3}];
 bq=Table[Together[Coefficient[Coefficient[ex,q1],s,k]],{k,0,3}];
 Print["Source temporal degree ",j,": ",Exponent[ex,s]];
 ee=Sum[aq[[k+1]] td[[k+1]]+bq[[k+1]]D[td[[k+1]],x],{k,0,3}];
 qvars=Table[Derivative[k][q][x],{k,0,4}];
 vvars=Table[Derivative[k][v][x],{k,0,3}];
 qc=Together/@(Coefficient[ee,#]&/@qvars);
 vc=Together/@(Coefficient[ee,#]&/@vvars);
 Put[{aq,bq,qc,vc},base<>"theory_direction1/cascade_operator_"<>ToString[j]<>".wl"];
 Get[base<>"investigation/public_inputs/AsymptoticRelations.wl"];
 rec=Together[(REEE/.{NUMERIC->0,\[Omega]1->w1,\[Omega]2->I*s})*
               (-1)^m*ThreeJSymbol[{j,j},{4,4},{j+4,-j-4}]];
 pk=If[j==2,"220","320"];ww=cx[inp["modes"][pk]["omega"]];sd=-I*ww;
 Print["Output poles ",j,": ",Factor[Denominator[rec]]];
 ra0=Together[Cancel[s*rec]/.s->0];
 rb0=Together[Cancel[(s-I*w1)*rec]/.s->I*w1];
 rp=PolynomialQuotient[Numerator[rec],Denominator[rec],s];
 Print["Output decomposition exact: ",FullSimplify[rec-rp-ra0/s-rb0/(s-I*w1)]==0];
 ra=ra0/.w1->ww;rb=rb0/.w1->ww;
 rpol=CoefficientList[rp,s]/.w1->ww;
 Print["Output polynomial degree ",j,": ",Length[rpol]-1];
 rpol=PadRight[rpol,Max[3,Length[rpol]]];
 samples={};
 Do[
  grid=N[(1-Cos[Pi*(Range[n]-1/2)/n])/2,70];
  {p,dp}=profile[pk,j,grid];
  qcc=MapThread[(qc/.{x->#1,p1->#2,dp1->#3,w1->ww})&,{grid,p,dp}];
  vcc=MapThread[(vc/.{x->#1,p1->#2,dp1->#3,w1->ww})&,{grid,p,dp}];
  ac=MapThread[(aq/.{x->#1,p1->#2,dp1->#3,w1->ww})&,{grid,p,dp}];
  bc=MapThread[(bq/.{x->#1,p1->#2,dp1->#3,w1->ww})&,{grid,p,dp}];
  AppendTo[samples,<|"n"->n,"sigma"->N[grid,17],
    "qcoeff"->Map[enc,qcc,{2}],"vcoeff"->Map[enc,vcc,{2}],
    "time_qcoeff"->Map[enc,ac,{2}],"time_dqcoeff"->Map[enc,bc,{2}]|>];
  Print["Grid ",j," ",n," coefficient max ",N[Max[Abs[Flatten[{qcc,vcc}]]],6]];
 ,{n,{24,32,40,48,64}}];
 AppendTo[rows,<|"spectator_l"->j,"omega"->enc[ww],
    "output_polynomial"->(enc/@rpol),"output_primitive"->enc[ra],
    "output_product_primitive"->enc[rb],"samples"->samples|>];
,{j,{2,3}}];
Export[base<>"theory_direction1/cascade_evolution_inputs.json",<|"rows"->rows|>,"RawJSON"];
Exit[];
