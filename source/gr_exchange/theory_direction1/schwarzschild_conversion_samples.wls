base = StringReplace[DirectoryName[DirectoryName[$InputFileName]], "\\" -> "/"] <> "/";
$MaxExtraPrecision=500;
inp=Import[base<>"investigation/overlap_inputs.json","RawJSON"];
cx[z_]:=SetPrecision[ToExpression[z[[1]]]+I ToExpression[z[[2]]],60];
enc[z_]:={Re[N[z,17]],Im[N[z,17]]};
poly[a_,x_]:=Fold[#1*x+#2&,Last[a],Reverse[Most[a]]];
JS[a_,b_,c_]:=ThreeJSymbol[a,b,c];
M=1;l1=2;m1=2;l2=4;m2=4;l=6;m=6;
Get[base<>"investigation/public_inputs/SourceEEE.wl"];
ss=SourceEEE/.{\[Psi]1l1e[r]->z1,\[Psi]1l2e[r]->z2,
 Derivative[1][\[Psi]1l1e][r]->dz1,Derivative[1][\[Psi]1l2e][r]->dz2};
freq1=cx[inp["modes"]["220"]["omega"]];
freq2=cx[inp["modes"]["440"]["omega"]];
ss=ss/.{\[Omega]1->freq1,\[Omega]2->freq2};
Get[base<>"investigation/public_inputs/AsymptoticRelations.wl"];
output=REEE/.{NUMERIC->0,\[Omega]1->freq1,\[Omega]2->freq2};
output=output*ThreeJSymbol[{2,2},{4,4},{6,-6}];
cc[j_]:=(j+2)*(j-1)*((j+2)*(j-1)+2)/6;
hh[j_]:=cc[j]+3*sigma^2*(1-sigma)/((j+2)*(j-1)+3*sigma);
az[j_,ww_]:=(hh[j]+2 I ww*(1-2*sigma^2))/(cc[j]+2 I ww);
bz[j_,ww_]:=-sigma^2*(1-sigma)/(cc[j]+2 I ww);
prof[key_,grid_]:=Module[{j,ww,aa,da,dda,p,dp,ddp,zz,dzz},
 j=inp["modes"][key]["l"];ww=cx[inp["modes"][key]["omega"]];
 aa=cx/@inp["modes"][key]["a"];
 da=Range[Length[aa]-1] Rest[aa];dda=Range[Length[da]-1]Rest[da];
 p=(poly[aa,1-#]&/@grid)/Total[aa];
 dp=-(poly[da,1-#]&/@grid)/Total[aa];
 ddp=(poly[dda,1-#]&/@grid)/Total[aa];
 zz=(az[j,ww]/.sigma->grid)*p+(bz[j,ww]/.sigma->grid)*dp;
 dzz=(D[az[j,ww],sigma]/.sigma->grid)*p+
       ((az[j,ww]+D[bz[j,ww],sigma])/.sigma->grid)*dp+
       (bz[j,ww]/.sigma->grid)*ddp;
 {zz,I ww*(1+2*grid-grid/(1-grid))*zz-grid^2*dzz/2}
];
rows={};
Do[
 grid=N[(1-Cos[Pi*(Range[n]-1/2)/n])/2,60];
 {p1,d1}=prof["220",grid];{p2,d2}=prof["440",grid];
 b=MapThread[-(ss/.{r->2/#1,z1->#2,z2->#3,dz1->#4,dz2->#5})/
                  (4 #1^2*(1-#1^2))&,{grid,p1,p2,d1,d2}];
 AppendTo[rows,<|"n"->n,"sigma"->N[grid,17],
   "source"->(enc/@b),"parent_product"->(enc/@(p1*p2))|>];
 Print["Conversion N=",n," endpoints ",N[b[[{1,-1}]],10]];
,{n,{32,48,64,80}}];
Export[base<>"theory_direction1/schwarzschild_conversion_samples.json",
 <|"output_reconstruction"->enc[output],"rows"->rows|>,"RawJSON"];
Exit[];
