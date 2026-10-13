Get[FileNameJoin[{DirectoryName[DirectoryName[$TestFileName]],"Kernel","QCS.wl"}]];
SetSystemOptions["ParallelOptions"->"ParallelThreadNumber"->1];
checks=QCSCheckData[];
VerificationTest[Abs[checks["infinity_mismatch_to_midpoint_0_400M"]-9.893781599611428*^-7]<10^-16,True,TestID->"Published infinity mismatch"]
VerificationTest[Max[Abs[checks["preparation_energies"]-{.3284468105600639,.3515542044860984}]]<10^-11,True,TestID->"Canonical preparation energies"]
VerificationTest[Max[Abs[checks["pointwise_difference_peak_percent"]-{8.10947781837,11.1231719292}]]<10^-8,True,TestID->"Published cubic separation"]
VerificationTest[Max[Abs[checks["relative_complex_L2_difference_errors"]-{1.0210320029301056*^-5,7.301918567353557*^-6}]]<10^-13,True,TestID->"Retained reconstruction diagnostics"]
VerificationTest[Abs[checks["hidden_direction_energy"]-1]<10^-11,True,TestID->"Unit energy of hidden direction"]
VerificationTest[Length[FileNames["*.mat",FileNameJoin[{QCSRoot,"source"}],Infinity]],13,TestID->"All thirteen reference array archives retained"]
kerr=QCSKerr[];
VerificationTest[Abs[kerr["at_chi_0_999999"]["abs_K0"]-1.5072533194642]<10^-11,True,TestID->"Published near-extremal pole factor"]
VerificationTest[kerr["fit_points"],70,TestID->"Kerr fit window"]
VerificationTest[Max[Abs[Values[kerr["exponents_chi_ge_0_9999"]]-{.997427542990024,.996554526226830,.996990565104261,1.002283931687785,.999866120814652,.007731944408151}]]<10^-10,True,TestID->"Six Kerr scaling exponents"]
nativeTargets=QCS`Private`kerrTargets[];paperTargets=QCS`Private`csvRows["figure_data/kerr_target_spin_diagnostics.csv"];
VerificationTest[Max[Abs[Flatten[Values/@nativeTargets]-Flatten[Lookup[#,Keys[First[nativeTargets]]]&/@paperTargets]]]<10^-9,True,TestID->"All target-spin Kerr and overtone diagnostics"]
pair=QCSPreparePair[];retained=QCSData["source/reduced_input/paired_preparations_fine.mat"];
VerificationTest[Max[Table[QCS`Private`relative[pair[key],retained[key]],{key,{"radiation_plus","radiation_minus","strain_plus","strain_minus","news_plus","news_minus"}}]]<10^-8,True,TestID->"Native pair propagation matches retained data"]
recovery=QCSReconstructPair[pair];retainedRecovery=QCSData["source/reduced_input/completion/recovery_waveforms.mat"];
VerificationTest[Max[Table[QCS`Private`relative[recovery[key],retainedRecovery[key]],{key,{"prediction_plus","prediction_minus"}}]]<10^-8,True,TestID->"Native reconstruction matches retained data"]
generated=QCSKernels[32,True];retained=QCSKernels[32];
VerificationTest[Max[Table[QCS`Private`relative[generated[key],retained[key]],{key,{"O","Od","F6","F7","Fd6","Fd7"}}]]<10^-8,True,TestID->"Regenerated kernels match retained kernels"]
gr=Quiet[QCSGR[],LinearSolve::luc];
VerificationTest[gr["maximum_factorization_error"]<10^-8,True,TestID->"GR isolated-pole factorization"]
VerificationTest[Max[Lookup[gr["production"],"field_change_residue_error"]]<10^-7,True,TestID->"Field redefinition with transported initial data"]
VerificationTest[Min[Lookup[gr["production"],"incorrect_zero_initial_data_error"]]>.1,True,TestID->"Incorrect initial-data reset is detected"]
