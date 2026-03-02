Journal of Exposure Science & Environmental Epidemiology
www.nature.com/jes
ARTICLE
Compare the marginal effects for environmental exposure and
biomonitoring data with repeated measurements and values
below the limit of detection
✉
I-ChenChen 1 , Stephen J.Bertke1 and CherylFairfield Estill1
ThisisaU.S.GovernmentworkandnotundercopyrightprotectionintheUS;foreigncopyrightprotectionmayapply2024
BACKGROUND: Environmental exposure andbiomonitoring datawithrepeated measurements fromenvironmentaland
occupationalstudiesarecommonlyright-skewedandinthepresenceoflimitsofdetection(LOD).However,existingmodelhasnot
been discussed forsmall-sample properties and highlyskewed datawith non-detects and repeatedmeasurements.
OBJECTIVE: Marginal modelingprovides analternative toanalyzinglongitudinal and cluster data,in whichtheparameter
interpretations arewith respect to marginalorpopulation-averaged means.
METHODS: Weoutlined thetheories of threemarginalmodels,i.e., generalizedestimating equations (GEE),quadratic inference
functions (QIF), andgeneralized methodof moments(GMM). With these approaches, weproposedto incorporate thefill-in
methods,includingsingleandmultiplevalueimputationtechniques,suchthatanymeasurementslessthanthelimitofdetection
are assigned values.
RESULTS: We demonstrated thatthe GEEmethod workswell intermsof estimating theregression parametersin smallsample
sizes, whiletheQIF and GMMoutperform inlarge-samplesettings, asparameter estimates areconsistent and haverelatively
smaller mean squarederror. Nospecific fill-in methodcan bedeemed superior aseach hasits ownmerits.
IMPACT:
● Marginalmodelingisfirstlyemployedtoanalyzerepeatedmeasuresdatawithnon-detects,inwhichonlythemeanstructure
needs tobe correctly provided toobtain consistent parameter estimates. Afterreplacing non-detects through substitution
methodsandutilizingsmall-samplebiascorrections,inasimulationstudywefoundthattheestimatingapproachesusedinthe
marginalmodels have correspondingadvantages under a widerangeofsample sizes. Wealso applied themodels to
longitudinal andcluster working examples.
Keywords: Marginalanalysis; Left censoring; Right skewness; Limitofdetection; Repeated measures; Environmentalexposure
Journal of Exposure Science& Environmental Epidemiology (2024)34:1018–1027; https://doi.org/10.1038/s41370-024-00640-7
INTRODUCTION the concentration of an analyte in a biological urine or serum
In environmental and occupational studies, repeated concentra- sample, or an environmental hand wipe or personal breathing
tion measurements not detected or falling below limits of zoneairsample,areoftensubjecttonon-detectsorleftcensoring
detection (LOD) of laboratory instruments are called left- and the data are skewed to the right. These measurements are
censored repeated measures data. The unquantified non-detects usuallycollectedfromthesamesubjectorthesamestudysite.In
are generally low-level concentrations with values between zero such cases, statistical modeling of exposure and biomonitoring
and LOD. Statistical models continue to arise for industrial data can be complicated when repeated measurements are
hygienists to analyze left-censored environmental exposure and collectedina cluster or longitudinalstudy.
biomonitoring data with repeated measures in cluster and Statistical methods have been continuously proposed to
longitudinal studies because the estimation of the effect of analyze left-censored data sets. The substitution or single value
exposure on risk of disease and the importance of within- and imputation method, e.g., assigning a value (LOD/2 or LOD/ 2)
between-worker variability in occupational exposure have been [1,2]formeasurementslessthantheLOD,iscommonlyadopted
increasingly acknowledged. Analytical results from laboratories, byindustrialhygienists.Unfortunately,thereisnouniquereplaced
environmental contaminants, and occupational exposures, e.g., valueforthissubstitutionorsinglevalueimputationmethod,and
1DivisionofFieldStudiesandEngineering,NationalInstituteforOccupationalSafetyandHealth,CentersforDiseaseControlandPrevention,Cincinnati,OH,USA.
✉
email:okv0@cdc.gov
Received:14March2023Revised:3January2024Accepted:4January2024
Publishedonline:22January2024

regression parameter estimation of the substitution method can such that any measurements less than the LOD are assigned
be biased for high censoring proportion [1]. This substitution values.Wealsowillconsidersmall-samplebiascorrectionsforthe
approachisalsonotadvisableunlesslessthan10%ofvaluesare empiricalcovarianceestimatorsofregressionparameterestimates
[25–31].
below the LOD [3]. Comparatively, the multiple random value Therefore, the resulting approaches will have the
imputation technique, e.g., creating imputed values based on potential to perform better than the existing methods in small-
throughoutscatterofthedataset,hasbeenadvocated[3–5].The
use of a maximum likelihood (ML) estimation approach has also compare the proposed methods and evaluate how well they
been shown to outperform other methods when the working estimate regression parameters for exposure and biomonitoring
modeliswellspecifiedandthesamplesizeislarge[6–8,20],but
datawithrepeatedmeasurementsunderarangeofsamplesizes
bias and imprecision are expected with the ML approach when andLODproportions.Finally,weillustratetheproposedmethods
sample size is small, i.e., fewer than 50 detectable values, and/or using a longitudinal chlorpyrifos exposure dataset and a cluster
flameretardant
censoring proportion is high, even though the distribution is biomonitoring dataset.
correctlyspecified[7,9,10].Inaddition,theMLmethodbasedon
whenthetruedatadistributionismis-specified[11].Recently,the
METHODS
β-substitutionmethodderivingthecalculationofaβfactorbased In this section, we described the marginal models popularly used for
comparable to the ML method, even when sample sizes were approaches proposed for filling in left-censored data without repeated
The paper evaluates these methods using simulation and two real-world examples (one longitudinal and one clustered occupational dataset).

### Simulation design (best-effort reconstruction)

- Sample sizes (`N`): 30, 100, and 500.
- Repeated measures per subject (`M`): 3.
- Censoring proportions: 10%, 20%, and 30%.
- Estimation approaches compared primarily: `GEE` and `QIF` (with `GMM` discussed separately).
- Substitution/imputation methods: `LOD/2`, `β`-substitution, multiple random imputation (with/without covariates), and `QQ`-plot ordered imputation.

### Reported simulation findings (best-effort reconstruction)

- `GEE` generally performs better in small samples.
- `QIF` shows advantages in larger samples.
- No single fill-in method is uniformly superior across all settings.
- Coverage probability and relative-efficiency patterns vary by censoring level and sample size.

Figure 3. Relative efficiency (`RE`) comparisons across estimation and fill-in methods.

Figure 4. Coverage probability comparisons across estimation and fill-in methods.

### Real-world applications (best-effort reconstruction)

The paper presents:

1. A longitudinal chlorpyrifos exposure example with repeated air measurements.
2. A clustered nail-salon biomonitoring example (`DPhP`), with repeated measurements within salons.

For these applications, the article reports generally similar effect directions across methods, with method-dependent differences in standard errors and efficiency.
Even when distributional assumptions are met, maximum-likelihood approaches can be biased or imprecise in small samples [20]. In contrast, marginal models can be used as an alternative for repeated measures and clustered data because they target population-averaged effects and can remain consistent when the mean structure is correctly specified, even if the working correlation is misspecified [21].

## Marginal models

This paper discusses three estimating approaches for correlated data:

- `GEE` (generalized estimating equations)
- `QIF` (quadratic inference functions)
- `GMM` (generalized method of moments)

The authors combine these approaches with practical fill-in methods for left-censored observations (values below the `LOD`), including:

- single-value substitution (`LOD/2`)
- `β`-substitution
- multiple random-value imputation (with and without covariates)
- multiple ordered-value imputation (`QQ-plot` based)

They also discuss small-sample covariance/bias corrections for robust standard error estimation [25–31].

Figure 1. Empirical bias comparison across estimation approaches (`GEE`, `QIF`) and substitution methods.

Figure 2. Empirical mean squared error (`MSE`) comparison across estimation approaches and substitution methods.


In the nail-salon example, workers who had not worked the previous day had lower pre-shift `DPhP` concentrations than those who worked the prior day. The `GEE` approach is emphasized for small cluster counts, and method selection is discussed using `CIC` (correlation information criterion), where smaller values indicate lower variability of regression estimates.

## Discussion

For repeated-measures exposure and biomonitoring data with non-detects, incorporating fill-in methods into marginal models (`GEE`, `QIF`, `GMM`) is feasible. The simulation results indicate:

- `GEE` tends to perform better in smaller samples.
- `QIF` tends to perform better in moderate-to-large samples.
- `GMM` can be unstable for inference in smaller samples when many moment conditions are used.

## Table 3

Parameter estimates, bias-corrected standard errors (in parentheses), and correlation-related summary values from the flame-retardant analysis.

| Model | Covariate contrast | LOD/2 | Beta-substitution | Imputation (no covariates) | Imputation (with covariates) | QQ-plot |
|---|---|---:|---:|---:|---:|---:|
| GEE | Two or more days ago vs previous day | −0.76 (0.17) | −0.72 (0.11) | −0.80 (0.09) | −0.59 (0.12) | −0.36 (0.43) |
| QIF | Two or more days ago vs previous day | −0.73 (0.08) | −0.69 (0.08) | −0.59 (0.10) | −0.91 (0.17) | −0.71 (0.21) |

Estimated correlation parameter in `Rᵢ` (GEE; as extracted):

| LOD/2 | Beta-substitution | Imputation (no covariates) | Imputation (with covariates) | QQ-plot |
|---:|---:|---:|---:|---:|
| −0.03 | −0.21 | −0.03 | −0.25 | −0.29 |

`CIC` (GEE; as extracted):

| LOD/2 | Beta-substitution | Imputation (no covariates) | Imputation (with covariates) | QQ-plot |
|---:|---:|---:|---:|---:|
| 1.67 | 0.95 | 1.35 | 0.89 | 1.18 |

Abbreviations: `GEE` = generalized estimating equations; `QIF` = quadratic inference functions; `CIC` = correlation information criterion.

## Conclusions

Only the mean structure and a working covariance/correlation structure are needed to fit these marginal approaches. After handling values below `LOD` and applying small-sample bias corrections, the extracted article indicates that:

- `GEE` is generally preferable for small samples.
- `QIF` can be preferable for larger samples.
- `LOD/2`, `β`-substitution, and multiple-imputation variants can all be reasonable, with performance depending on sample size, censoring level, and distributional features.
## Data availability

Detailed information for the two working examples is reported in the cited articles [32, 33]. The original paper indicates that simulation/application code and helper functions are provided in supplementary material.

## References
substitutionaredifficulttoperformstandardnormalitytests,they
areeasiertoimplementandcalculate.Anotherlimitationisthatall 1. HornungRW,ReedLD.Estimationofaverageconcentrationinthepresenceof
nondetectablevalues.ApplOccupEnvironHyg.1990;5:46–51.
2. BurstynI,TeschkeK.Studyingthedeterminantsofexposure:areviewofmeth-
assumedtobeindependent.Thereisnodifferenceintheestimate
ods.AmIndHygAssocJ.1999;60:57–72.
3. LubinJH,ColtJS,CamannD,DavisS,CerhanJR,SeversonRK,etal.Epidemiologic
disregard of the correlation will result in positively biased SE evaluation of measurement data in the presence of detection limits. Environ
estimates, i.e., incorrect estimates of the sampling variability. In HealthPerspect.2004;112:1691–6.

4. HuybrechtsT,ThasO,DewulfJ,VanLangenhovH.Howtoestimatemoments 33. EstillCF,SloneJ,MayerAC,ChenIC,ZhouM,LaGuardiaMJ,etal.Assessmentof
andquantilesofenvironmentaldatasetswithnondetectedobservations?Acase TriphenylPhosphate(TPhP)exposuretonailsalonworkersbyair,handwipe,and
studyonvolatileorganiccompoundsinmarinewatersamples.JChromatogrA. urineanalysis.IntJHygEnvironHealth.2021;231:113630.
2002;975:123–33. 34. WindmeijerF.Afinitesamplecorrectionforthevarianceoflinearefficienttwo-
5. BaccarelliA,PfeifferR,ConsonniD,PesatoriAC,BonziniM,PattersonDGJr,etal. stepGMMestimators.JEcon.2005;126:25–51.
Handlingofdioxinmeasurementdatainthepresenceofnondetectablevalues: 35. KauermannG,CarrollRJ.Anoteontheefficiencyofsandwichcovariancematrix
overview of available methods and their application in the Seveso chloracne estimation.JAmStatAssoc.2001;96:1387–96.
study.Chemosphere.2005;60:898–906. 36. LittleRJA,RubinDB.StatisticalAnalysiswithMissingData.2nded.Hoboken,New
6. AmemiyaT.Regressionanalysiswhenthedependentvariableistruncatednor- Jersey:JohnWileyandSons;2002.
mal.Econometrica.1973;41:997–1016. 37. Hargarten PM, Wheeler DC. miWQS: Multiple Imputation Using Weighted
7. Helsel DR. Fabricating data: how substituting values for nondetects can ruin Quantile Sum Regression. R package version 0.4.4; 2021. https://CRAN.R-
results,andwhatcanbedoneaboutit.Chemosphere.2006;65:2434–9. roject.org/package=miWQS
8. HewettP,GanserGH.Acomparisonofseveralmethodsforanalyzingcensored 38. RCoreTeam.R:alanguageandenvironmentforstatisticalcomputing.Vienna,
data.AnnOccupHyg.2007;51:611–32. Austria:RFoundationforStatisticalComputing;2022.https://www.R-project.org
9. GilliomRJ,HelselDR.Estimationofdistributionalparametersforcensoredtracelevel 39. Jacqmin-GaddaH,Thi´ebautR.Analysisofleft-censoredlongitudinaldatawith
waterqualitydata1.estimationtechniques.WaterResourRes.1986;22:135–46. applicationtoviralloadinHIVinfection.Biostatistics.2000;1:355–68.
10. Helsel DR, Cohn TA. Estimation of descriptive statistics for multiply censored 40. Hin LY, Wang YG. Working-correlation-structure identification in generalized
waterqualitydata.WaterResourRes.1988;24:1997–2004. estimatingequations.StatMed.2009;28:642–58.
11. ShoariN,DubéJS,ChenouriS.Estimatingthemeanandstandarddeviationof 41. WestgatePM.Criterionforthesimultaneousselectionofaworkingcorrelation
environmentaldatawithbelowdetectionlimitobservations:Consideringhighly structureandeithergeneralizedestimatingequationsorthequadraticinference
skeweddataandmodelmisspecification.Chemosphere.2015;138:599–608. functionapproach.BiometricalJ.2014;56:461–76.
12. GanserGH,HewettP.Anaccuratesubstitutionmethodforanalyzingcensored 42. DigglePJ,HeagertyPJ,LiangKY,ZegerSL.TheAnalysisofLongitudinalData,2nd
data.JOccupEnvironHyg.2010;7:233–44. ed.NewYork:OxfordUniversityPress;2002.
13. Pleil JD. QQ-plots for assessing distributions of biomarker measurements and 43. NeweyWK,SmithRJ.HigherorderpropertiesofGMMandgeneralizedempirical
generatingdefensiblesummarystatistics.JBreathRes.2016;10:035001. likelihoodestimators.Econometrica.2004;72:219–55.
14. PleilJD.Imputingdefensiblevaluesforleft-censored‘belowlevelofquantitation’ 44. Westgate PM. A bias-corrected covariance estimator for improved inference
(LoQ)biomarkermeasurements.JBreathRes.2016;10:045001. whenusinganunstructuredcorrelationwithquadraticinferencefunctions.Stat
15. Thi´ebaut R, Jacqmin-Gadda H. Mixed models for longitudinal left-censored ProbabLett.2013;83:1553–8.
repeatedmeasures.ComputMethodsProgBiomed.2004;74:255–60. 45. SASInstituteInc.SAS/STAT9.3UsersGuide.SASInstituteInc.,Cary,NC;2011.
16. Thi´ebaut R, Guedj J, Jacqmin-Gadda H, Chenê G, Trimoulet P, Neau D, et al. 46. ChenIC,WestgatePM.Anovelapproachtoselectingclassificationtypesfortime-
Estimation of dynamical model parameters taking into account undetectable dependentcovariatesforthemarginalanalysisoflongitudinaldata.StatMethods
markervalues.BMCMedResMethodol.2006;6:38. MedRes.2018;28:3176–86.
17. VaidaF,LiuL.Fastimplementationfornormalmixedeffectsmodelswithcen- 47. ChenLS,PrenticeRL,WangP.ApenalizedEMalgorithmincorporatingmissing
soredresponse.JComputGraphStat.2009;18:797–817. datamechanismforgaussianparameterestimation.Biometrics.2014;70:312–22.
18. JinY,HeinMJ,DeddensJA,HinesCJ.Analysisoflognormallydistributedexpo-
suredatawithrepeatedmeasuresandvaluesbelowthelimitofdetectionusing
SAS.AnnOccupHyg.2011;55:97–112. ACKNOWLEDGEMENTS
19. LeidelNA,BuschKA,LynchJR.Occupationalexposuresamplingstrategymanual WewouldliketothankthepeoplefromtheDivisionofFieldStudiesandEngineering
(DHEW [NIOSH] publication no. 77-173). Cincinnati, OH: National Institute for atCDC’sNationalInstituteforOccupationalSafetyandHealthwhoassistedinthe
OccupationalSafetyandHealth;1977. study.Thefindingsandconclusionsinthismanuscriptarethoseoftheauthorsand
20. HelselDR.Lessthanobvious:statisticaltreatmentofdatabelowthedetection do not necessarily represent the official position of the National Institute for
limit.EnvironSciTechnol.1990;24:1766–74.
OccupationalSafetyandHealth,CentersforDiseaseControlandPrevention.
21. LiangKY,ZegerSL.Longitudinaldataanalysisusinggeneralizedlinearmodels.
Biometrika.1986;73:13–22.
22. WangYG,CareyV.Workingcorrelationstructuremisspecification,estimationand
AUTHORCONTRIBUTIONS
covariatedesign:implicationsforgeneralisedestimatingequationsperformance.
Biometrika.2003;90:29–41. ICCwasresponsiblefordesigningstatisticalmethods,conductingasimulationstudy,
23. Hansen LP. Large sample properties of generalized method of moments esti- analyzing two real-world datasets, interpreting simulation and application results,
mators.Econometrica.1982;50:1029–54. producingtablesandfigures,draftingtheinitialmanuscript,revisingthemanuscript,
24. Qu A, Lindsay BG, Li B. Improving generalised estimating equations using
andapprovingthefinalversionofmanuscript.SJBcontributedtointerpretationsof
quadraticinferencefunctions.Biometrika.2000;87:823–36. simulation and application results, revised manuscript, provided feedback, and
25. Mancl LA, DeRouen TA. A covariance estimator for GEE with improvedsmall-
approvedthefinalversion.CFEcontributedtodatacurationandextraction,revised
sampleproperties.Biometrics.2001;57:126–34.
manuscript,providedfeedback,andapprovedthefinalversion.Additionally,Whitney
26. WestgatePM.Abias-correctedcovarianceestimateforimprovedinferencewith F.TannerandYu-ChengChenreviewedthepaperandprovidedhelpfulfeedback.
quadraticinferencefunctions.StatMed.2012;31:4003–22.
27. WestgatePM.Abiascorrectionforcovarianceestimatorstoimproveinference
with generalized estimating equations that use an unstructured correlation COMPETINGINTERESTS
matrix.StatMed.2013;32:2850–8. Theauthorsdeclarenocompetinginterests.
28. WestgatePM.Acovariancecorrectionthataccountsforcorrelationestimationto
improvefinite-sampleinferencewithgeneralizedestimatingequations:Astudy
on its applicability with structured correlation matrices. J Stat Comput Simul. ADDITIONAL INFORMATION
2016;86:1891–1900.
Supplementary information The online version contains supplementary material
29. ChenIC,WestgatePM.Improvedmethodsforthemarginalanalysisoflongitudinal availableathttps://doi.org/10.1038/s41370-024-00640-7.
datainthepresenceoftime-dependentcovariates.StatMed.2017;36:2533–46.
30. FordWP,WestgatePM.Improvedstandarderrorestimatorformaintainingthe
CorrespondenceandrequestsformaterialsshouldbeaddressedtoI-ChenChen.
validityofinferenceinclusterrandomizedtrialswithasmallnumberofclusters.
BiometricalJ.2017;59:478–95.
Reprints and permission information is available at http://www.nature.com/
31. Ford WP, Westgate PM. A comparison of bias-corrected empirical covariance
reprints
estimators withgeneralizedestimatingequations insmall-samplelongitudinal
studysettings.StatMed.2018;37:4318–29. Publisher’snoteSpringerNatureremainsneutralwithregardtojurisdictionalclaims
32. HinesCJ,DeddensJA.Determinantsofchlorpyrifosexposuresandurinary3,5,6-trichloro- inpublishedmapsandinstitutionalaffiliations.
2-pyridinollevelsamongtermiticideapplicators.AnnOccupHyg.2001;45:309–21.
