
* Survival analysis resource : https://communities.sas.com/t5/SAS-Communities-Library/Kaplan-Meier-Survival-Plotting-Macro-NEWSURV/ta-p/479747;
* The macro available from the tutorial above has to be run before running this program;
FILENAME FILE '/home/u35553484/sasuser.v94/TMB - survival analysis/NCI_NeverSmoker_n131_20210812_TMB_SAS.xlsx';

PROC IMPORT DATAFILE=FILE
	DBMS=XLSX replace
	OUT=NSLC;
	GETNAMES=YES;
	LABEL Ftime="Days";
	RANGE="Sheet 1$A1:AO93";
RUN;

PROC MEANS data=NSLC median MaxDec=3;
var complete_WGS_TMB;
run;

data NSLC_form;
set NSLC;
length sex_recode $6 Path_stage $8 TMB_class_med $4 TMB_class_you $4 Histo $30;
* Sex variable conversion;
if      sex = 2 then sex_recode = "Female";
else if sex = 1 then sex_recode = "Male";
* Pathological stage variable conversion;
if pathological_stage = "1A1" then Path_stage = "I";
else if pathological_stage = "1A2" then Path_stage = "I";
else if pathological_stage = "1A3" then Path_stage = "I";
else if pathological_stage = "1B" then Path_stage = "I";
else if pathological_stage = "2B" then Path_stage = "II";
else if pathological_stage = "2A" then Path_stage = "II";
else if pathological_stage = "3A" then Path_stage = "III & IV";
else if pathological_stage = "3B" then Path_stage = "III & IV";
else if pathological_stage = "4A" then Path_stage = "III & IV";

* Histology conversion;
if histology = 2 then Histo = "squamous cell carcinoma";
else if histology = 3 then Histo = "adenocarcinoma";
else if histology = 4 then Histo = "carcinoid tumor";
else if histology = 6 then Histo = "adenosquamous carcinoma";
else if histology = 7 then Histo = "sarcomatoid carcinoma";

* WGS TMB median cutoff is set at 1.57 mutation/Mb;
if complete_WGS_TMB >= 1.57 then TMB_class_med = "high";
else if complete_WGS_TMB < 1.57 then TMB_class_med = "low";
* WGS TMB Youden's cutoff is set at 1.70 mutation/Mb;
if complete_WGS_TMB >= 1.70 then TMB_class_you = "high";
else if complete_WGS_TMB < 1.70 then TMB_class_you = "low";
* WGS TMB Youden's cutoff is set at 1.70 mutation/Mb;
if complete_WES_TMB >= 1.20 then TMB_class_you_WES = "high";
else if complete_WES_TMB < 1.20 then TMB_class_you_WES = "low";

run;

data NSLC;
   set NSLC_form (drop=sex rename=(sex_recode=Sex));
run;

ods listing gpath="/home/u35553484/sasuser.v94/TMB - survival analysis";
ods graphics / reset=all border=off;

/* Skipping code block. Remove the "/" (in this line and in the line at the end of the block) to run code.


* Kaplan-Meier (KM) survival plots 
* 1) sex based KM;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=Sex, CLASSREF=Female,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=60, TIMEDX=days, DISPLAY=legend EV_N MEDIAN timelist PVAL, CLASSORDER = 2 1,
RISKLIST=0 to 250 by 20,RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_sex, DPI=300, PTABSIZE=11,
HEIGHT = 6in, WIDTH = 6in, STATCOLOR = 1, PLOTPVAL = wald, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, AUTOALIGN=bottomleft);

* 2) pathological stage based KM;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=Path_stage, CLASSREF=I,
COLOR=(CX2A9D8F CXF50025 CXA63A9B CXE9C46A), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, 
XLABEL=Time after surgery (days), YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=60, TIMEDX=days,
DISPLAY=legend EV_N MEDIAN timelist PVAL, RISKLIST=0 to 250 by 20, RISKCOLOR=1, REFLINES=medians, PTABSIZE=11,
REFLINEAXIS=both, PLOTNAME=Surv_Analysis_stage, DPI=300, HEIGHT = 6in, WIDTH = 6in, PLOTPVAL = wald, STATCOLOR = 1,
UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, AUTOALIGN=bottomleft);

* Versions with Patients-at-risk table inside the plot
* 1) sex based KM
* 2) pathological stage based KM;
* 3) histology based KM;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=Sex, CLASSREF=Female,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=60, TIMEDX=days, DISPLAY=legend EV_N MEDIAN timelist PVAL, CLASSORDER = 2 1,
RISKLIST=0 to 250 by 20,RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_sex_v2, DPI=300, PTABSIZE=11, 
HEIGHT = 6in, WIDTH = 6in, STATCOLOR = 1, PLOTPVAL = wald, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, AUTOALIGN=bottomleft,
LOCATION=outside, RISKLOCATION=inside, BORDER = 0);
*/

/*
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=Path_stage, CLASSREF=I,
COLOR=(CX2A9D8F CXF50025 CXA63A9B CXE9C46A), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, 
XLABEL=Time after surgery (days), YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=1826, TIMEDX=days, CLASSDESC = Pathological stage,
DISPLAY=legend EV_N MEDIAN timelist PVAL, RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, AUTOALIGN=bottomleft,
REFLINEAXIS=both, PLOTNAME=Surv_Analysis_stage_OS, DPI=300, HEIGHT = 6in, WIDTH = 6in, PLOTPVAL = wald, STATCOLOR = 1, 
UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);
*/

/*
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=Histo, CLASSREF=adenocarcinoma,
COLOR=(CX2A9D8F CXF50025 CXA63A9B CXE9C46A), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, 
XLABEL=Time after surgery (days), YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=60, TIMEDX=days, CLASSDESC = Histology Type,
DISPLAY=legend EV_N MEDIAN timelist PVAL, RISKLIST=0 to 250 by 20, RISKCOLOR=1, REFLINES=medians, PTABSIZE=11, AUTOALIGN=bottomleft,
REFLINEAXIS=both, PLOTNAME=Surv_Analysis_histology_v2, DPI=300, HEIGHT = 8in, WIDTH = 8in, PLOTPVAL = wald, STATCOLOR = 1, 
UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

* Another type of KM, but evaluating low vs high TMB according to WGS TMB median cutoff;
* 1) Patients-at-risk table outside the plot;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_med, CLASSREF=low,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=60, TIMEDX=days, DISPLAY=legend EV_N MEDIAN timelist PVAL, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 250 by 20,RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_med, DPI=300, PTABSIZE=11,
HEIGHT = 6in, WIDTH = 6in, STATCOLOR = 1, PLOTPVAL = wald, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, AUTOALIGN=bottomleft, BORDER = 0);

* 2) Patients-at-risk table inside the plot;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_med, CLASSREF=low,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, 
XLABEL=Time after surgery (days), YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=60, TIMEDX=days, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
DISPLAY=legend EV_N MEDIAN timelist PVAL, RISKLIST=0 to 250 by 20, RISKCOLOR=1, REFLINES=medians, PTABSIZE=11, AUTOALIGN=bottomleft,
REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_med_v2, DPI=300, HEIGHT = 6in, WIDTH = 6in, PLOTPVAL = wald, STATCOLOR = 1,
UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

* KM, but evaluating low vs high TMB according to WGS TMB Youden index cutoff;
* 1) Patients-at-risk table outside the plot;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=60, TIMEDX=days, DISPLAY=legend EV_N MEDIAN timelist PVAL, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 250 by 20,RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_you, DPI=300, PTABSIZE=11,
HEIGHT = 6in, WIDTH = 6in, STATCOLOR = 1, PLOTPVAL = wald, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, AUTOALIGN=bottomleft, BORDER = 0); 
*/

/*
* 2) Patients-at-risk table inside the plot;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, 
XLABEL=Time after surgery (days), YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=1826, TIMEDX=days, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
DISPLAY=legend EV_N MEDIAN timelist PVAL, RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, PTABSIZE=11, AUTOALIGN=bottomleft,
REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_you_KM, DPI=300, HEIGHT = 7in, WIDTH = 7in, PLOTPVAL = wald, STATCOLOR = 1,
UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);
*/

/* This section is not used anymore, because there is no generalized cutoff anymore.
* KM, but evaluating low vs high TMB according to WGS TMB generalized cutoff (combination of median and Youden's); ** NOT USED ANYMORE **
* 1) Patients-at-risk table outside the plot;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_gen, CLASSREF=low,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=60, TIMEDX=days, DISPLAY=legend EV_N MEDIAN timelist PVAL, CLASSDESC = TMB Classification,
CLASSORDER = 2 1, RISKLIST=0 to 250 by 20,RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_gen, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 6in, WIDTH = 6in, STATCOLOR = 1, PLOTPVAL = wald, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, AUTOALIGN=bottomleft, BORDER = 0);

* 2) Patients-at-risk table inside the plot;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_gen, CLASSREF=low,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, 
XLABEL=Time after surgery (days), YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=60, TIMEDX=days, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
DISPLAY=legend EV_N MEDIAN timelist PVAL, RISKLIST=0 to 250 by 20, RISKCOLOR=1, REFLINES=medians, PTABSIZE=11, AUTOALIGN=bottomleft,
REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_gen_v2, DPI=300, HEIGHT = 6in, WIDTH = 6in, PLOTPVAL = wald, STATCOLOR = 1,
UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);
(end of skipping block) */
/*
* Cox Proportional-Hazards Model plots using OS data. Covariates are TMB threshold (high or low TMB) and pathological stage (1, 2 or 3);
* 1) Patients-at-risk table outside the plot;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low, PLOTPVALMV = wald,
HRTIES=EFRON, CLASSCOV=Path_stage,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 250 by 20, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox, DPI=300, PTABSIZE=11,
HEIGHT = 6in, WIDTH = 6in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, AUTOALIGN=bottomleft, BORDER = 0);
*/

/*
* 2) Patients-at-risk table inside the plot;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
PLOTPVALMV = wald, PLOTPVAL = wald, HRTIES=EFRON, CLASSCOV=Path_stage, LSIZE = 12, PARSIZE = 11, XTICKVALSIZE = 10, YTICKVALSIZE = 12,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVAL PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_OS, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 7in, WIDTH = 7in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

* Cox Proportional-Hazards Model plots using OS data with median as the cutoff. 
* Covariates are TMB threshold (high or low TMB) and pathological stage (1, 2 or 3);
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_med, CLASSREF=low,
PLOTPVALMV = wald, PLOTPVAL = wald, HRTIES=EFRON, CLASSCOV=Path_stage, LSIZE = 12, PARSIZE = 11, XTICKVALSIZE = 10, YTICKVALSIZE = 12,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVAL PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_OS_med, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 7in, WIDTH = 7in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

* Cox Proportional-Hazards Model plots using RpFS data. Covariates are TMB threshold (high or low TMB) and pathological stage (1, 2 or 3);
* 	Patients-at-risk table inside the plot;
%newsurv(DATA=NSLC, TIME=time_RpFS, CENS=RpFS_indicator, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low, 
PLOTPVALMV = wald, PLOTPVAL = wald, HRTIES=EFRON, CLASSCOV=Path_stage, LSIZE = 12, PARSIZE = 11, XTICKVALSIZE = 10, YTICKVALSIZE = 12,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Relapse-free probability, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVAL PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 6000 by 500, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_RpFS, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 7in, WIDTH = 7in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

* Kaplan-Meier plot using RpFS data. Covariates are TMB threshold (high or low TMB) and pathological stage (1, 2 or 3);
* 	Patients-at-risk table inside the plot;
%newsurv(DATA=NSLC, TIME=time_RpFS, CENS=RpFS_indicator, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, 
XLABEL=Time after surgery (days), YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=1826, TIMEDX=days, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
DISPLAY=legend EV_N MEDIAN timelist PVAL, RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, PTABSIZE=10.8, AUTOALIGN=bottomleft,
REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_KM_RpFS, DPI=300, HEIGHT = 7in, WIDTH = 7in, PLOTPVAL = wald, STATCOLOR = 1,
UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);
*/

/* (Skipping block)
* Cox Proportional-Hazards Model plots. Histology is added as a covariate compared to last Cox plots;
* 1) Patients-at-risk table outside the plot;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
HRTIES=EFRON, CLASSCOV=Path_stage,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 250 by 20, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_histology, DPI=300, PTABSIZE=11,
HEIGHT = 6in, WIDTH = 6in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, AUTOALIGN=bottomleft, BORDER = 0);
* 2) Patients-at-risk table inside the plot;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
HRTIES=EFRON, CLASSCOV=Path_stage Sex histology, CONTCOV=age, 
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 250 by 20, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_histology_v2, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 6in, WIDTH = 6in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

* Cox Proportional-Hazards Model plots based on histology instead of TMB classification;
* Patients-at-risk table inside the plot;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=Histo, CLASSREF=adenocarcinoma, PLOTPVALMV = wald,
HRTIES=EFRON, CLASSCOV=Path_stage,
COLOR=(CX2A9D8F CXF50025 CXA63A9B CXE9C46A CX3D87B8), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVALMV, CLASSDESC = Histology Type,
RISKLIST=0 to 250 by 20, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_histology_Cox, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 7in, WIDTH = 7in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

(end of skipping block) */

* Plots for WES classification.
* 1) KM;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you_WES, CLASSREF=low,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, 
XLABEL=Time after surgery (days), YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=1826, TIMEDX=days, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
DISPLAY=legend EV_N MEDIAN timelist PVAL, RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, PTABSIZE=10, AUTOALIGN=bottomleft,
REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_KM_OS_WES, DPI=300, HEIGHT = 7in, WIDTH = 7in, PLOTPVAL = wald, STATCOLOR = 1,
UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);
* 2) Cox;
%newsurv(DATA=NSLC, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you_WES, CLASSREF=low,
PLOTPVALMV = wald, PLOTPVAL = wald, HRTIES=EFRON, CLASSCOV=Path_stage, LSIZE = 12, PARSIZE = 11, XTICKVALSIZE = 10, YTICKVALSIZE = 12,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVAL PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_OS_WES, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 7in, WIDTH = 7in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);




* Survival analysis on adenocarcinoma and carcinoid hystology patients. Cox regression models;
FILENAME A_C_file '/home/u35553484/sasuser.v94/TMB - survival analysis/NCI_NeverSmoker_n131_20210812_TMB_adeno.xlsx';

PROC IMPORT DATAFILE=A_C_file
	DBMS=XLSX replace
	OUT=A_C_TMB;
	GETNAMES=YES;
	LABEL Ftime="Days";
	RANGE="Sheet 1$A1:AO83";
RUN;

* Adenocarcinoma dataset;
data A_TMB_form;
set A_C_TMB;
length sex_recode $6 Path_stage $8 TMB_class_you $4;
* Only keep adeno histology (3);
if histology = 4 then DELETE;
* Sex variable conversion;
if      sex = 2 then sex_recode = "Female";
else if sex = 1 then sex_recode = "Male";
* Pathological stage variable conversion;
if pathological_stage = "1A1" then Path_stage = "I";
else if pathological_stage = "1A2" then Path_stage = "I";
else if pathological_stage = "1A3" then Path_stage = "I";
else if pathological_stage = "1B" then Path_stage = "I";
else if pathological_stage = "2B" then Path_stage = "II";
else if pathological_stage = "2A" then Path_stage = "II";
else if pathological_stage = "3A" then Path_stage = "III & IV";
else if pathological_stage = "3B" then Path_stage = "III & IV";
else if pathological_stage = "4A" then Path_stage = "III & IV";

* WGS TMB Youden's cutoff for Adeno; 
* This cutoff is used instead of a generalized cutoff because of the lack of an approximated number (e.g. 1 mut/Mb for all histological types);
if complete_WGS_TMB >= 1.7 then TMB_class_you = "high";
else if complete_WGS_TMB < 1.7 then TMB_class_you = "low";

run;

data A_TMB;
   set A_TMB_form (drop=sex rename=(sex_recode=Sex));
run;

/* (Skipping block)
* Carcinoid dataset;
data C_TMB_form;
set A_C_TMB;
length sex_recode $6 Path_stage $8 TMB_class_you $4;
* Only keep adeno histology (3);
if histology = 3 then DELETE;
* Sex variable conversion;
if      sex = 2 then sex_recode = "Female";
else if sex = 1 then sex_recode = "Male";
* Pathological stage variable conversion;
if pathological_stage = "1A1" then Path_stage = "I";
else if pathological_stage = "1A2" then Path_stage = "I";
else if pathological_stage = "1A3" then Path_stage = "I";
else if pathological_stage = "1B" then Path_stage = "I";
else if pathological_stage = "2B" then Path_stage = "II";
else if pathological_stage = "2A" then Path_stage = "II";
else if pathological_stage = "3A" then Path_stage = "III & IV";
else if pathological_stage = "3B" then Path_stage = "III & IV";
else if pathological_stage = "4A" then Path_stage = "III & IV";

* WGS TMB Youden's cutoff for Carcino; 
* This cutoff is used instead of a generalized cutoff because of the lack of an approximated number (e.g. 1 mut/Mb for all histological types);
if complete_WGS_TMB >= 0.45 then TMB_class_you = "high";
else if complete_WGS_TMB < 0.45 then TMB_class_you = "low";

run;

data C_TMB;
   set C_TMB_form (drop=sex rename=(sex_recode=Sex));
run;

ods listing gpath="/home/u35553484/sasuser.v94/TMB - survival analysis";
ods graphics / reset=all border=off;


* KM for adenocarcinomas;
%newsurv(DATA=A_TMB, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, 
XLABEL=Time after surgery (days), YTYPE=PPT, YLABEL=Proportion surviving, TIMELIST=1826, TIMEDX=days, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
DISPLAY=legend EV_N MEDIAN timelist PVAL, RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, PTABSIZE=11, AUTOALIGN=bottomleft,
REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_KM_adeno_OS, DPI=300, HEIGHT = 7in, WIDTH = 7in, PLOTPVAL = wald, STATCOLOR = 1,
UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);
(end of skipping block)*/

/*
* Cox regression models
* 1) Adenocarcinoma patients-at-risk table inside the plot;
%newsurv(DATA=A_TMB, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
HRTIES=EFRON, CLASSCOV=Path_stage, PLOTPVAL = wald, PLOTPVALMV = wald, LSIZE = 12, PARSIZE = 11, XTICKVALSIZE = 10, YTICKVALSIZE = 12,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVAL PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_adeno_OS, DPI=300, PTABSIZE=11pt, AUTOALIGN=bottomleft,
HEIGHT = 7in, WIDTH = 7in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

* TITLE = Cox regression model for adenocarcinoma `with Youden index (1.50 mutations/MB));

/* (Skipping block)
* 2) Carcinoid patients-at-risk table inside the plot;
%newsurv(DATA=C_TMB, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
HRTIES=EFRON, CLASSCOV=Path_stage Sex histology, CONTCOV=age, 
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=10, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 250 by 20, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_carcino_v2, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 6in, WIDTH = 6in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

* TITLE = Cox regression model for carcinoid `with Youden index (0.45 mutations/MB));
(end of skipping block)*/

/*
* Cox Proportional-Hazards Model plots using RpFS data. Covariates are TMB threshold (high or low TMB) and pathological stage (1, 2 or 3);
* 	Patients-at-risk table inside the plot;
%newsurv(DATA=A_TMB, TIME=time_RpFS, CENS=RpFS_indicator, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low, 
PLOTPVALMV = wald, PLOTPVAL = wald, HRTIES=EFRON, CLASSCOV=Path_stage, LSIZE = 12, PARSIZE = 11, XTICKVALSIZE = 10, YTICKVALSIZE = 12,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=500, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Relapse-free probability, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVAL PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 6000 by 500, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_adeno_RpFS, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 7in, WIDTH = 7in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);
*/

* Cox analyses for EGFR positive and driver positive adjustements.;
FILENAME drivers '/home/u35553484/sasuser.v94/TMB - survival analysis/NCI_NeverSmoker_n131_20210812_TMB_92_drivers.xlsx';

PROC IMPORT DATAFILE=drivers
	DBMS=XLSX replace
	OUT=NSLC_drivers;
	GETNAMES=YES;
	LABEL Ftime="Days";
	RANGE="Sheet1$A1:AX93";
RUN;

data NSLC_form_drivers;
set NSLC_drivers;
length sex_recode $6 Path_stage $8 TMB_class_med $4 TMB_class_you $4 Histo $30;
* Sex variable conversion;
if      sex = 2 then sex_recode = "Female";
else if sex = 1 then sex_recode = "Male";
* Pathological stage variable conversion;
if pathological_stage = "1A1" then Path_stage = "I";
else if pathological_stage = "1A2" then Path_stage = "I";
else if pathological_stage = "1A3" then Path_stage = "I";
else if pathological_stage = "1B" then Path_stage = "I";
else if pathological_stage = "2B" then Path_stage = "II";
else if pathological_stage = "2A" then Path_stage = "II";
else if pathological_stage = "3A" then Path_stage = "III & IV";
else if pathological_stage = "3B" then Path_stage = "III & IV";
else if pathological_stage = "4A" then Path_stage = "III & IV";

* Histology conversion;
if histology = 2 then Histo = "squamous cell carcinoma";
else if histology = 3 then Histo = "adenocarcinoma";
else if histology = 4 then Histo = "carcinoid tumor";
else if histology = 6 then Histo = "adenosquamous carcinoma";
else if histology = 7 then Histo = "sarcomatoid carcinoma";

* WGS TMB Youden's cutoff is set at 1.70 mutation/Mb;
if complete_WGS_TMB >= 1.70 then TMB_class_you = "high";
else if complete_WGS_TMB < 1.70 then TMB_class_you = "low";
run;

data NSLC_drivers;
   set NSLC_form_drivers (drop=sex rename=(sex_recode=Sex));
run;

ods listing gpath="/home/u35553484/sasuser.v94/TMB - survival analysis";
ods graphics / reset=all border=off;

/*
* Model 1 : OS = TMB + EGFR;
%newsurv(DATA=NSLC_drivers, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
PLOTPVALMV = wald, PLOTPVAL = wald, HRTIES=EFRON, CLASSCOV=EGFR, LSIZE = 12, PARSIZE = 11, XTICKVALSIZE = 10, YTICKVALSIZE = 12,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVAL PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_OS_EGFR, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 7in, WIDTH = 7in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

* Model 2 : OS = TMB + stade + EGFR;
%newsurv(DATA=NSLC_drivers, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
PLOTPVALMV = wald, PLOTPVAL = wald, HRTIES=EFRON, CLASSCOV=Path_stage EGFR, LSIZE = 12, PARSIZE = 11, XTICKVALSIZE = 10, YTICKVALSIZE = 12,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVAL PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_OS_stage_EGFR, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 7in, WIDTH = 7in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

* Model 3 : OS = TMB + drivers;
%newsurv(DATA=NSLC_drivers, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
PLOTPVALMV = wald, PLOTPVAL = wald, HRTIES=EFRON, CLASSCOV=driver, LSIZE = 12, PARSIZE = 11, XTICKVALSIZE = 10, YTICKVALSIZE = 12,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVAL PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_OS_driver, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 7in, WIDTH = 7in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);

* Model 4 : OS = TMB + stade + drivers;
%newsurv(DATA=NSLC_drivers, TIME=time_os, CENS=VitalStatus, CEN_VL=0, SUMMARY=0, CLASS=TMB_class_you, CLASSREF=low,
PLOTPVALMV = wald, PLOTPVAL = wald, HRTIES=EFRON, CLASSCOV=Path_stage driver, LSIZE = 12, PARSIZE = 11, XTICKVALSIZE = 10, YTICKVALSIZE = 12,
COLOR=(CX2A9D8F CXF50025), PATTERN=solid, SYMBOLSIZE=6pt, symbolweight=bold, XINCREMENT=250, XLABEL=Time after surgery (days),
YTYPE=PPT, YLABEL=Proportion surviving, TIMEDX=days, DISPLAY=legend EV_NMV MEDIAN HRMV PVAL PVALMV, CLASSDESC = TMB Classification, CLASSORDER = 2 1,
RISKLIST=0 to 8000 by 500, RISKCOLOR=1, REFLINES=medians, REFLINEAXIS=both, PLOTNAME=Surv_Analysis_TMB_Cox_OS_stage_drivers, DPI=300, PTABSIZE=11, AUTOALIGN=bottomleft,
HEIGHT = 7in, WIDTH = 7in, STATCOLOR = 1, UNIFORMHEIGHT = 1, RISKROWWEIGHTS = 0.04, LOCATION=outside, RISKLOCATION=inside, BORDER = 0);
*/