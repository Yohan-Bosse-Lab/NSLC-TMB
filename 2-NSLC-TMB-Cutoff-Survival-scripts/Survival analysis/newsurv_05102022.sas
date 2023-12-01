/*------------------------------------------------------------------*
| MACRO NAME  : newsurv
| SHORT DESC  : Creates Kaplan-Meier survival plots with summary
|               built into plot.  Can also output summary table.
*------------------------------------------------------------------*
| CREATED BY  : Meyers, Jeffrey                 (01/20/2013 3:00)
*------------------------------------------------------------------*
| VERSION UPDATES:
| 10.1 - 04/28/2022
|  Added DRAW_UNDERLINES option to turn the underlines on or off
|  Added UL_SIZE, UL_PATTERN and UL_COLOR to modify underline attributes
|  Program corrections
| 10.0 - 08/27/2021
|  Added BY parameter to overlay subset KM curves in graph window and
    subset statistics within plot/table output.  This is only available
    when CLASS is not missing.
    Also added BYORDER to control order of BY values and BYLABEL to 
    modify label of BY variable
|  Added PVAL_INTER and PVAL_INTERMV to print the interaction p-value
    between BY variable and CLASS variable.  Options to print have been
    added to DISPLAY and TABLEDISPLAY
|  Added MODEL_STATS_DISPLAY parameter to determine whether model level
    statistics are printed in columns or rows.
| 9.6 - 02/11/2021
|  Made corrections to adjusted survival calculations
| 9.5 - 01/25/2021
|  Made updates for RMST_PVAL calculations.  Can only be calculated 
|    when CLASS levels don't differ across strata levels.
| 9.4 - 12/22/2020
|  Corrected issue with reference text for hazard ratios
| 9.3 - 12/17/2020
|  Corrected issue with colors when METHOD=INVWTS.
|  Updated coding for efficiency
|  Added option LOGRANK_ONESIDED for PLOTPVAL
| 9.2 - 12/08/2020
|  Corrected issue with adjusted medians for METHOD=DIRECT.
| 9.1 - 11/19/2020
|  Corrected issue with RISKLIST when METHOD=CIF and multiple numbers
|    are entered
| 9.0 - 10/27/2020
|  Added Restricted Means Survival Time options
|     RMST, RMTL (restricted means time lost), and RMST_PVAL are available in
|       DISPLAY and TABLEDISPLAY parameters
|     TAU option is added to specify time-point
|     Options for headers/widths available
|  Added numeric values for all statistics to dataset created by OUT
|    Gives raw values for statistics (estimates and upper/lower CI)
|  Restructured how statistics are merged together within the macro
| 8.8 - 09/28/2020
|  Fixed issue: TABLEFOOTNOTE now accepts line breaks
| 8.7 - 06/29/2020 
|  Updated logic for discrete attribute maps when there is no class variable
| 8.6 - 06/08/2020 
|  Updated logic for discrete attribute maps
| 8.5 - 05/22/2020
|  Corrected an issue when class values end in a % symbol
| 8.4 - 04/13/2020
|  Corrected an error of grabbing the wrong timelist values to be shown
     in the plot when covariates in CLASSCOV or CONTCOV
|  Corrected some display issues in footnotes when using INVWTS or DIRECT methods
|  Corrected display issues within LISTING output.
|  Added DEBUG option to keep temporary datasets and see notes.
| 8.3 - 01/15/2020
|  Corrected TIMELIST duplicating estimates across class values
| 8.2 - 01/13/2020
|  Corrected TIMELIST showing the correct time-points in table
|     when no CLASS is specified.
| 8.1 - 12/10/2019
|  Updated LINEOBSMAX to be 1000000 within the macro
|  Made changes to allow macro to work better in Windows SAS
| 8.0 - 10/21/2019
|  Added two additional methods to METHOD parameter:
      INVWTS: Adjusted survival using inverse weights methods
      DIRECT: Adjusted survival using direct methods
|  Added PLOT_UNADJUST option to hide/display unadjusted curves
      when using METHOD=INVWTS/DIRECT
|  Added PATTERN_ADJUST to set patterns for adjusted curves
|  Changed appearance of TIMELIST in table
|  LISTTIMEPOINTS now can also remove time points from TIMELIST in table
|  When changing ALPHA the correct percentage CI is listed in the headers
|  Updated graphing methods.  Output table (OUTP) has been modified
| 7.4 - 08/06/2019
|  Added ability to add gridlines to the graph.
|  Added CENSORMARKERS=2 to draw the censor markers but not include the
      symbol in the legend.
| 7.3 - 04/12/2019
|  Updated documentation for readability
| 7.2 - 08/28/2018
|  Changed dataset for stratified Gray test in CIF methods
| 7.1 - 05/23/2018
|  Fixed logic that cause the where clause to not be applied properly.
| 7.0 - 02/11/2018
|  Macro is no longer compatible with SAS 9.4M2 or earlier versions
|  Added error check for when WHERE parameter fails
|  Modified how temporary dataset is put together to avoid errors when 
     two submitted variables have the same name
|  Updated CIF section to use the LIFETEST calculation instead of 
     data step calculations
|  Multivariate counts (events/totals) now pulled from a LIFETEST procedure
| 6.62 - 2/05/2018
|  Fixed error within new LISTING report section.
|  Added the ALPHA parameter to change confidence limit types.
| 6.61 - 12/05/2017
|  Fixed error within new LISTING report section.
| 6.6 - 12/02/2017
|  Updated summary table section.  Now has the same appearance across 
     ODS destinations and works in the LISTING output window.
|  Added TSUBTITLEWIDTH option to control the width of the subtitle 
     column.
| 6.5 - 09/26/2017
|  Fixed logic issue when using CIF methods
| 6.4 - 08/09/2017
|  Fixed issue in 9.3 when using class variable.
| 6.3 - 05/29/2017
|  Logic corrections
| 6.3 - 05/25/2017
|  Updated documentation notes
| 6.2 - 05/23/2017
|  New parameter added: SYMBOLWEIGHT. Allows the censor symbols to be
      set to bold making them thicker.
| 6.1 - 05/21/2017
|  New parameters added to change graph styles:
      TRANSPARENT (9.4M3+): Changes the background of the plot to be transparent
      AXISCOLOR: Changes the color of the plot axes and border lines
      BACKGROUND: Changes the background color of the plot
      FONTCOLOR: Changes the color of the text within the plot
|  RISKLABELDLM default changed to missing to match new patients-at-risk format
|  SAS 9.4 now allows the destination to be EXCEL and POWERPOINT
|  SAS 9.4: Reverts to new EMF format when any transparency is specified
|  Added CINDEX and CINDEXMV as options to the DISPLAY and TABLEDISPLAY 
      parameters.  These are only allowed when METHOD=KM.
|  ODS LISTING is disabled for the Report Table.  The format of the table is not designed
      to be displayed by the output window, but rather the results window.
| 6.0 - 05/14/2017
|  Multiple changes to patients-at-risk:
     Values will no longer be cut-off by sizing/offsets when SAS 9.3+ 
       and RISKLOCATION=BOTTOM
     PARDISPLAY added to allow the display of patients at risk (PAR), 
       number of cumulative censors (NCENS), and number of cumulative events (NEVENTS).
       Combinations of PAR_NCENS and PAR_NEVENTS are allowed
     Headers for number of cumulative censors and number of cumulative events added
     RISKLOCATION default changed to BOTTOM from nothing.
|  Warning when using cumulative incidence removed. Will monitor for errors.
*------------------------------------------------------------------*
| PURPOSE
|
| This macro runs survival analysis on a time variable with or 
| without a class variable.  The analyses that are run are:
| number of patients, number of events, median time to event,
| hazard ratios, and p-values (logrank, score, and likelihood-
| ratio).  These analyses are stored in a summary dataset that
| is then be used to output a Kaplan-Meier survival plot
| with the summary information listed on the plot (which output
| appears in the plot is customizable), in a summary table, or
| both.  The plot is extremely customizable and the time can
| be transformed by a factor.  The summary table can be stored
| and added to with further calls of the macro (when newtable=0),
| allowing for the output of multiple models.
| 
| Multiple models can be ran in one one macro call, either to add
| multiple models statistics to the summary table or to plot 
| multiple Kaplan-Meier curves into a lattice plot.  Each of these
| plots and lattices are customizable individually by separating
| each model's attributes with a | delimiter.  The number of models
| run is determined by the number of time variables given.  For Example
| to run two models through the macro enter TIME=var1|var2.  See example
| 6 for more details.  The same time variable can be listed multiple 
| times.  For other parameters, if the macro parameter will not change
| across models, the parameter only needs to be listed once and will
| be assumed to be the same across all models.  For example, if 
| TIME=var1|var2 and CLASS=class1, the macro will assume that class1 
| will be the class variable for both models.  If instead the class 
| parameter is: CLASS=class1|class2, then class1 will be used for 
| model 1 and class2 will be used for model 2.  Blank values are also
| accepted.
|
| 1.0: REQUIRED PARAMETERS
|     DATA = dataset that contains the time variable, censor variable,
|           and optional class, adjusting variables, and stratification variables.
|     TIME = Variable containing time to event information
|     CENS = Numeric variable containing event (coded as a binary variable)
|     CEN_VL = Numeric constant representing value of a non-event in CENS.
|              Default = 0.
|     METHOD = Determines the method for calculating the survival estimates.  Options are
|              KM for Kaplan-Meier, CIF for cumulative incidence function (Competing Risks),
|              INVWTS for inverse-weights adjusted survival, and DIRECT for direct adjusted survival.
|              Default is KM.
|     EV_VL = Numeric constant representing the event of interest in a competing risks
|             analysis.  There is no default for this variable.  It is only required when
|             method = CIF.
|
| 2.0: Optional PARAMETERS
|   2.1: Output controlling options
|       2.1.1: Output datasets
|           OUT = A name for the output dataset for the statistical report table
|           OUTP = A name for the output dataset for the plot dataset.
|       2.1.2: Plot and Table switches
|           PLOT = A flag variable to turn printing the plot off (0) or on (1)
|                  Default = 1 (on), options = 1 (on) or 0 (off).
|           PLOT_UNADJUST = A flag variable to turn printing the unadjusted curves off (0) or on (1)
|                           when METHOD=INVWTS or DIRECT.   Default = 1 (on), options = 1 (on) or 0 (off).
|                           NOTE: 0 is only available when METHOD=INVWTS or DIRECT
|           SUMMARY = A flag variable to display a report table at the end of
|                     the macro off (0) or on (1)
|                     Default = 1 (on), options = 1 (on) or 0 (off).
|           NEWTABLE = A flag variable to determine if a new summary output table is created
|                      or if the results will be saved into a previously created dataset.
|                      1 = New table made, 0= previous dataset.  Default=1.
|       2.1.3: Outputting to a document options
|           DESTINATION = Type of ODS output when creating a document. 
|                         Default is RTF, options are RTF, PDF, HTML.
|           OUTDOC = Filepath with name at end to send the output.
|                    Example: ~/ibm/example.doc
|       2.1.4: Image file options
|         GPATH = Determines the path that the image is saved to.  Defaults to the path 
|                 the SAS session is opened to. 
|                 NOTE: ODS LISTING has to be enabled to create an image file.
|         PLOTNAME = Names the image file.  Plot will be saved per the GPATH parameter.  
|                    Default is _surv.
|         PLOTTYPE = Determines the image type of the plot.  Default is png, options
|                    are png, tiff, jpeg, emf, gif.  
|                    NOTE: There is special code added for TIFF plots in order to make 
|                          appropriately sized image files.  If DPI settings are too high
|                          depending on operating system and SAS version, this may not work.
|                    NOTE2: Special code is made for SAS 9.4 for EMF files.  This is due to SAS
|                           changing the default drivers for EMF files in 9.4, but this change
|                           causes the EMF files to not build properly when trying to convert to
|                           Windows objects.  Code given by Martin Mincey to add registry keys
|                           is used to temporarily add registry keys, then these are removed at
|                           the end of the macro.  This only occurs when making EMF files in SAS 9.4
|         DPI = Determines the dots per inch of the image file.  Default=200.
|         HEIGHT = Sets the height of plot window.  Default is 5in.  Set by a
|                  numeric constant followed by px or in.  Must be in for TIFF.
|         WIDTH = Sets the width of plot window.  Default is 7in.  Set by a
|                 numeric constant followed by px or in. Must be in for TIFF.
|         SVG = Turns scalable vector graphics on or off.  
|               Possible Scalable Vector Graphics formats are EMF within or not within RTF, 
|               PDF, and HTML.  In order to activate the scalable vector graphics, the 
|               DESTINATION must be used in conjunction with the SVG parameter.  To create
|               SVG EMF files use DESTINATION=RTF and PLOTTYPE=EMF.  To create SVG PDF files
|               use DESTINATION=PDF.  To create SVG HTML files use DESTINATION=HTML.
|               Default is 0 (off).  options are 1 or 0
|         TIFFDEVICE = Determines the GDEVICE to use when making TIFF plots.  Default is TIFFP.
|                      Options can be found with PROC GDEVICE catalogue=sashelp.devices;run;
|         ANTIALIASMAX = Maximum threshold to keep line smoothing activated.  Set to
|                        arbitrarily large number if large file.
|   2.2: Input Dataset Modifiers
|       XDIVISOR = Numeric constant to transform the time variable into other units. 
|                  Will divide the time variable before Kaplan-Meier survival estimates are
|                  computed.
|       WHERE = Gives a where clause to subset the DATA dataset down.  Type
|               exactly as would be in a procedure step. Does not change original dataset
|               Example: where=age>70
|   2.3.1: Class Variable Options
|       CLASS = Variable used to subset patients for comparison (Character or
|               Numeric).  Can be formatted; numeric or character.
|       CLASSORDER = List of numbers corresponding to the preferred order of the
|                    alphabetically sorted class variable formatted values.
|                    Example: Values = A, B, C.  CLASSORDER = 2 1 3 would cause
|                    order to be B, A, C.
|                    NOTE: When BY variable is enabled the ` symbol can be used to deliminate
|                          between BY levels
|       DESC = Reverses the order that the class variable is displayed.
|              Default is 0. Can be 0 or 1.  Enter desc to reverse order. 
|   2.3.2: By Variables (Subset CLASS Variables)
|       BY = Variable used to subset the CLASS variable further.  Will calculate statistics and 
|            draw curves within each level of the BY variable.  CLASS levels do not have to match
|            between BY levels.  Opens ability to calculate interaction p-values between CLASS and BY
|            variables.  CLASSORDER and CLASSREF can be specified separately between BY levels.
|       BYORDER = List of numbers corresponding to the preferred order of the
|                 alphabetically sorted class variable formatted values.
|                 Example: Values = A, B, C.  BYORDER = 2 1 3 would cause
|                 order to be B, A, C.
|       PVAL_INTER = Designates the p-value to use for interaction p-value tests.  Options are WALD,
|                    SCORE and LR (likelihood ratio).  Default is LR.
|       PVAL_INTERMV = Designates the adjusted p-value to use for interaction p-value tests.  Options are WALD,
|                      SCORE and LR (likelihood ratio).  Default is LR.
|   2.4: Statistical Modeling Options
|       2.4.1: General Modeling Options
|           ALPHA = Determines the alpha for confidence intervals.  Default is 0.05.
|                   Note that this does not change any headers in the output.
|           STRATA = Variable(s) separated by spaces to use as a stratification    
|                    in Cox models.  Included in the STRATA statement in PHREG and
|                    creates stratified LOGRANK and WILCOXON p-values when
|                    PLOTPVAL=LOGRANK or WILCOXON
|           NMODELS = Sets the number of models computed within the macro. Default=1.
|       2.4.2: Rounding Options
|           HRDIGITS = Number of significant digits to show for the hazard ratios and confidence bounds.  Default is 2.
|           KMESTDIGITS = Number of significant digits to show for the hazard ratios and confidence bounds.  AUTO
|                         will determine the number of significant digits based on the YTYPE parameter.  Default is AUTO.
|           MEDIANDIGITS = Number of significant digits to show for the median time-to-event and confidence bounds. 
|                          Default is 1.
|           PVALDIGITS = Number of significant digits to show for all p-values.  Default is 4.
|           CINDEXDIGITS = Number of significant digits to show for all c-indexes.  Default is 2.
|           RMSTDIGITS = Number of significant digits to show for all restricted means survival time (RMST) and 
|                        restricted means time lost (RMTL).  Default is 1.
|       2.4.1: Kaplan-Meier Modeling
|           SREVERSE = Flag variable to model 1-Survival instead of Survival
|                      Options are 1 (yes) and 0 (no).  Default=0.
|                      NOTE: This is not the same as cumulative incidence.  Set METHOD=CIF
|                            for cumulative incidence/competing risks.
|           TIMELIST = Numeric time-points to collect Kaplan-Meier survival
|                      estimates w/confidence intervals.  Can be entered as 
|                      numeric values separated by spaces, or in a list format
|                      (example: 0 to 60 by 10).  Time-points should match
|                      the TRANSFORMED time-scale of time variable if XDIVISOR is used.
|           TAU = Determines the tau (time-point) for computing Restricted means survival time and
|                 restricted means time lost.  Defaults to largest observed time. Time units should match
|                 the TRANSFORMED time-scale of time variable if XDIVISOR is used.
|           RMST_CI = Determines if CI is included with RMST/RMTL values.  Default is 1 (yes).  Options are 1 or 0.
|           CONFTYPE = Method of computing confidence intervals for Kaplan-Meier
|                      survival estimates. Default is LOG, options are: LOG,
|                      ASINSQRT, LOGLOG, LINEAR, LOGIT.
|           LANDMARK = Gives either a variable or a number to landmark the TIME
|                      variable by for the analysis.  Number must be greater than 0.
|       2.4.1: Cumulative Incidence Modeling
|           CIFVAR = Sets the methods for CIF calculation of the variance.  Options are COUNT (Counting method)
|                    and DELTA (Marubini's delta method)  Default is COUNT.
|       2.4.1: P-values (Requires CLASS Variable)
|           PLOTPVAL = Type of P-value to display in plot and summary document. 
|                      Default is logrank, options are: score, logrank, logrank_onesided, lr, wilcoxon, and wald.
|                      LR stands for Likelihood-ratio. logrank_onesided is logrank/2
|           PLOTPVALMV = Type of adjusted P-value to display in plot and summary document. 
|                        Default is Score, options are: score, wald, and lr.
|                        LR stands for Likelihood-ratio. Requires CLASSCOV or CONTCOV
|       2.4.1: Cox Modeling (Requires CLASS variable)
|           CLASSREF = Value to use as a reference for hazard ratios.  Must match
|                      exact value of the class variable after formatting.
|                      NOTE: When BY variable is being used the ` symbol is used to deliminate between
|                            BY Levels.
|           HRTIES = Determines the method for dealing with ties in the PROC PHREG model 
|                    statement.  Default is BRESLOW.  Options are: BRESLOW, DISCRETE, EFRON, 
|                    and EXACT.
|           CLASSCOV = List of discrete variables to be used as adjusting covariates in
|                      multivariate hazard ratio models.  These are included in the PROC
|                      PHREG CLASS statement. Must be variable names in list separated by 
|                      spaces.
|           CONTCOV = List of continuous variables to be used as adjusting covariates in
|                     multivariate hazard ratio models.  These are not included in the PROC
|                     PHREG CLASS statement. Must be variable names in list separated by 
|                     spaces.
|           REFHRTEXT = Text to be shown for reference level within the hazard ratio column.  
|                       Default is REF.
|           REFPTEXT = Text to be shown for reference level within the covariate p-value column.  
|                      Default is --.
|   2.5: Plot Modification Options
|       2.5.1: Axis Options
|           2.5.1.1: Label Options
|               XLABEL = Sets a label for the x-axis.
|               YLABEL = Sets a label for the y-axis.
|               LFAMILY = Sets the font family for the x/y labels. Default is Albany AMT.   
|               LSIZE = Sets the font size for the text in the x/y labels. Default=10pt.
|                       Must be followed by pt.
|               LWEIGHT = Sets the weight of the text in the x/y labels. Default=bold.
|                         Options = medium or bold.
|           2.5.1.2: Tick-value Options
|               YTYPE = Determines whether Kaplan-Meier survival estimates are in
|                       percentages or proportions.  Default is pct, options are: pct and ppt.
|               XMIN/YMIN = Designates the minimum for the axis. Default=0.
|               XMAX = Designates the maximum for the x-axis. If left missing this will 
|                      calculated as the maximum time value rounded up to the next 
|                      number divisible by 5.
|               YMAX = Designates the maximum for the y-axis. If left missing will be set
|                      to either 1 or 100 depending on YTYPE.
|               XINCREMENT = Designates the distance between tick marks on the X-axis. If 
|                            left missing this will be calculated as (XMAX-XMIN)/5.
|               YINCREMENT = Designates the distance between tick marks on the Y-axis.  If
|                            left missing will be set to 0.1 or 10 depending on YTYPE.
|               XMINOFFSET/YMINOFFSET = Designates the amount of space the plot cannot take up at
|                                       the minimum side of the axis. Range between [0-1). 
|                                       Default=blank.  Blank will be automatically calculated.
|               XMAXOFFSET/YMAXOFFSET = Designates the amount of space the plot cannot take up at
|                                       the maximum side of the axis. Range between [0-1). 
|                                       Default=blank.  Blank will be automatically calculated.
|               XTICKVALFAMILY/YTICKVALFAMILY = Sets the font for the axis tick values. Default=Albany AMT 
|               XTICKVALSIZE/YTICKVALSIZE = Sets the font size for axis tick values. Default=8pt
|               XTICKVALWEIGHT/YTICKVALWEIGHT = Sets the font weight for the axis tick values. Default=normal, bold=bold.
|       2.5.2: Lines and Symbol Options
|           2.5.2.1: Survival/CIF Curves
|               COLOR = A list of colors separated by spaces to color lines in the plot.
|                       Default is black.  If only one color is listed, then the lines will
|                       change in pattern.  If multiple colors are listed, then all the
|                       lines will be solid unless PATTERN is specified.
|               PATTERN = A list of line patterns separted by spaces to set the line types
|                         in the plots.  Default is AUTO (picks numbers if only one color,
|                         does solid if multiple colors).  Options are to do numbers between
|                         1 and 46, or: SOLID, SHORTDASH, MEDIUMDASH, LONGDASH,
|                         MEDIUMDASHSHORTDASH, DASHDASHDOT, DASHDOTDOT, DASH, LONGDASHSHORTDASH,
|                         DOT, THINDOT, SHORTDASHDOT, and MEDIUMDASHDOTDOT. 
|               PATTERN_ADJUST = A list of line patterns separted by spaces to set the line types
|                                in the adjusted survival plots.  Default is AUTO (picks numbers if only one color,
|                                does solid if multiple colors).  Options are to do numbers between
|                                1 and 46, or: SOLID, SHORTDASH, MEDIUMDASH, LONGDASH,
|                                MEDIUMDASHSHORTDASH, DASHDASHDOT, DASHDOTDOT, DASH, LONGDASHSHORTDASH,
|                                DOT, THINDOT, SHORTDASHDOT, and MEDIUMDASHDOTDOT. If AUTO is used and PLOT_UNADJUST=1
|                                the macro will pick different patterns than the unadjusted curves.
|               LINESIZE = Size of the lines in the Kaplan-Meier curve.  Default = 1pt.          
|           2.5.2.2: Censor Symbols
|               CENSORMARKERS = A flag variable to turn on the display of censor marks on the plot. 
|                               Options are 1 (on),0 (off), and 2 (on but not in legend).  Default=1.   
|               SYMBOLSIZE = Size of the censor markers.  Default = 2pt.   
|               SYMBOLWEIGHT = Weight of the censor markers.  Options are NORMAL and BOLD.  Default = NORMAL. 
|       2.5.3: Font Options 
|           2.5.3.1: Plot specific titles/footnotes (Can be assigned to individual plots in multi-cell graphs)
|               TITLE = Sets the title for the plot.
|               TITLEALIGN = Sets the horizontal alignment for the title.  Options are left, right and center.  
|                            Default=center
|               TFAMILY = Sets the font family for the title. Default is Albany AMT.
|               TSIZE = Sets the font size for the text in the plot title. Default=12pt. Must be followed by pt.
|               TWEIGHT = Sets the weight of the text in the plot title. Default=bold. Options = normal or bold.
|               FOOTNOTE = Creates footnotes in the individual plot panes.  Multiple footnotes can be created by 
|                          separating them with the ` delimiter.
|               FOOTNOTEALIGN = Sets the horizontal alignment for the footnotes. Options are left, right and center.  
|                               Default=center
|               FNFAMILY = Sets the font family for the title. Default is Albany AMT.
|               FNSIZE = Sets the font size for the text in the footnotes. Default=12pt. Must be followed by pt.
|               FNWEIGHT = Sets the weight of the text in the plot title. Default=bold. Options = medium or bold.
|           2.5.3.2: Overall image titles/footnotes (Only one Title/Footnote can be assigned to the whole graph)
|               OVTITLE = Sets the overall title for lattice plots
|               OVTFAMILY = Sets the font family for the overall title. Default is Albany AMT.
|               OVTSIZE = Sets the font size for the text in the overall plot title. Default=12pt. Must be followed by pt.
|               OVTWEIGHT = Sets the weight of the text in the overall plot title. Default=bold. Options = medium or bold.
|               OVFOOTNOTE = Sets the overall footnote for lattice plots
|               OVFNFAMILY = Sets the font family for the overall footnote. Default is Albany AMT.
|               OVFNSIZE = Sets the font size for the text in the overall plot footnote. Default=12pt. Must be followed by pt.
|               OVFNWEIGHT = Sets the weight of the text in the overall plot footnote. Default=bold. Options = medium or bold.
|       2.5.4: Plot Summary Table Options
|           2.5.4.1: Summary Table Location Option
|               AUTOALIGN = Set a list separated by spaces of locations for the summary
|                           table in the plot.  Default is topright bottomleft.  Options
|                           are topleft, top, topright, left, center, right, bottomleft,
|                           bottom, bottomright. The plot will attempt to fit the table
|                           into the alignment that least interferes with the plot lines
|                           starting from left to right.
|                           NOTE: There is a bug with listing multiple locations that can
|                                 cause the table to be displayed in multiple locations if 
|                                 there is collision with the graph.  If this happens
|                                 change the list to one location.
|               LOCATION = Sets the summary table in the plot to be inside the plot window
|                          (inside) or outside the plot window (outside).
|           2.5.4.1: Statistical Output Options
|               DISPLAY = A list of metrics to display in the plot window. Any combination
|                         of the following may be entered separated by spaces: 
|                             LEGEND: legend for the Kaplan-Meier curves
|                             TOTAL: Total sample size
|                             EVENT: Total number of events
|                             EV_N: Combined events and total in the format EVENTS/Total
|                             N_EV: Combined events and total in the format TOTAL (EVENTS)
|                             MEDIAN: Median time-to-event and confidence interval
|                             TIMELIST: Time-point event-free rates, confidence intervals, and time-points
|                             HR: Hazard ratio and confidence interval
|                             PVAL: P-value specified by the PLOTPVAL parameter
|                             COVPVAL: P-value from Cox model parameters table
|                             CINDEX: C-index from the Cox model (not available for CIF methods)
|                             RMST: Restricted Means Survival Time (Only available when METHOD=KM)
|                             RMTL: Restricted Means Time Lost (Only available when METHOD=KM)
|                             RMST_PVAL: P-value comparing RMST or RMTL values (Only available when METHOD=KM)
|                             PVAL_INTER: Interaction P-value of BY and CLASS variables (Only available when BY is used)
|                             TABLECOMMENTS: User added annotation text specified by TABLECOMMENTS parameter
|                         When CLASSCOV or CONTCOV are specified, the following items can have the MV suffix
|                         added for the multivariate versions:
|                             TOTALMV, EVENTMV, EV_NMV, N_EVMV, HRMV, PVALMV, PVAL_INTERMV, COVPVALMV, and CINDEXMV
|                         The order that items are written in this parameter determines the order they 
|                         appear on the plot with one exception, Legend will always come first if listed.
|                         The default is set to STANDARD, which will select default options depending on
|                         plot method and SAS version.
|               MODEL_STATS_DISPLAY = Determines if model level statistics (p-value, c-index, etc) are printed
|                                     as rows or columns in the plot summary table.  Options are 1 (rows),
|                                     2 (columns), or 3 (columns except for interaction p-value).  Default is 1.
|               TABLECOMMENTS = A text string or series of text strings that will show up in
|                               the bottom of the plot summary table.  These text strings are
|                               manually entered into the plot and show up as typed.  These
|                               comments can be broken into multiple lines by splitting them
|                               with the ` delimiter (e.g. Comment 1`Comment2).
|               CLASSVALALIGN = Sets the horizontal alignment for the listed class values.
|                               Options are left,right, and center.  Default=center.
|               CLASSDESC = Functions differently depending on whether CLASS is missing or not.
|                           No CLASS variable: Gives a description to be used in the legend and
|                           patients-at-risk table label CLASS Variable present: Gives a column header 
|                           for the class levels in the plot summary table.  Using %STR( ) will make this 
|                           header blank.  Using a ` delimiter will cause a line break
|               LEGENDLINELENGTH = Sets the length of the legend lines within the summary table.
|                           Default is null.  Example is 0.5in.  Only available in SAS 9.4M1 or later.
|               LISTTIMEPOINTS = Flag variable to determine if the column showing the time-points 
|                                for the Kaplan-Meier event-free rates is displayed.  
|                                Options are 1 (on) and 0 (off).  Default is 1.
|           2.5.4.2: Font Options
|               PTABFAMILY = Sets the font family for the text in the plot.  Default is Albany AMT.
|               PTABSIZE = Sets the font size for the text in the plot.  Default is 8pt. Must be followed by pt.
|               STATCOLOR = A flag variable to enable the text for the statistical values
|                           to be colored to match the color of the corresponding plot lines.
|                           Options are 1 (on) and 0 (off).  Default is 0.
|           2.5.4.3: Column Header Options
|               prefixHEADER = The header for each item in DISPLAY can be modified by entering the 
|                              DISPLAY item followed by HEADER.  Options are:
|                                  LEGENDHEADER, TOTALHEADER, EVENTHEADER, EV_NHEADER, N_EVHEADER, MEDIANHEADER,
|                                  TIMELISTHEADER, KMESTHEADER, HRHEADER, PVALHEADER, COVPVALHEADER,
|                                  CINDEXHEADER, TOTALMVHEADER, EVENTMVHEADER, EV_NMVHEADER, N_EVMVHEADER,
|                                  HRMVHEADER, PVALMVHEADER, COVPVALMVHEADER, CINDEXMVHEADER, RMSTHEADER, RMTLHEADER,
|                                  RMST_PVALHEADER
|                              NOTE: KMESTHEADER is for the Kaplan-Meier estimate and confidence interval 
|                                    portion of the TIMELIST display item.  TIMELISTHEADER is for the 
|                                    time-points portion of the TIMELIST display item.
|                              NOTE: The ` symbol acts as a line break for the headers.
|       2.5.5: Patients-at-Risk Options
|           2.5.5.1: Time and location options
|               RISKLIST = Numeric time-points to collect number of patients at risk.
|                          Can be entered as numeric values separated by spaces, or in
|                          a list format (example: 0 to 60 by 10). 
|                          These must be in the TRANSFORMED time units if XDIVISOR is used.
|               RISKLOCATION = Location for the number at risk to show on the plot.  The
|                              The default is bottom, the options are: bottom (below
|                              the x-axis) and INSIDE (above x-axis below the plot).
|           2.5.5.1: Display options
|               PARDISPLAY = Determines which values are displayed in the patients-at-risk tables.
|                            Options are:
|                                PAR: Number of Patients-at-risk 
|                                NCENS: Cumulative number of censors
|                                NEVENTS: Cumulative number of events
|                                PAR_NCENS: Patients-at-Risk (Cumulative number of censors)
|                                PAR_NEVENTS: Patients-at-Risk (Cumulative number of events)
|                            Listing multiple items will display in the order they are listed.
|               PARALIGN = Determines the alignment of the patients-at-risk subtitle.
|                          Default is CENTER.  Options are LEFT|CENTER|RIGHT|LABELS.  LABELS
|                          will place subtitle above class values when RISKLABELLOCATION=LEFT.
|               RISKCOLOR = A flag variable that causes the patients-at-risk numbers to 
|                           match the colors of the plot lines.  Options are 1 or 0.
|                           Default is 0.
|               RISKROWWEIGHTS = When the RISKLOCATION=BOTTOM option is selected this sets
|                                the weights for the plot window and patients-at-risk table.
|                                Gives the weight value uniformly  to each class level.  
|                                Default is 0.025, range is (0,1)
|               UNIFORMHEIGHT = When any of the requested plots have a patients-at-risk table with
|                               RISKLOCATION=BOTTOM, the space below each plot's x-axis will be
|                               uniformly set based on the maximum number of rows and maximum value
|                               of RISKROWWEIGHTS.  Default is 0 (off).  Values are 1 (On) or 0 (Off).
|           2.5.5.2: Label options
|               RISKLABELLOCATION = Sets the location of the class level labels for the
|                                   patients-at-risk tables.  Default is LEFT.  Options are LEFT, ABOVE, and nothing.
|               RISKLABELALIGN = Sets the alignment for the patients-at-risk table when RISKLABELLOCATION=ABOVE.  
|                                Default is LEFT. Options are LEFT, CENTER, or RIGHT.
|               RISKLABELDLM = Determines the delimiter between the Risk label and the patients-at-risk
|                              RISKLABELLOCATION is set to LEFT. Default is -.
|               RISKLABELWEIGHT = Sets the font weights of the class level labels for the
|                               patients-at-risk tables.  Default is NORMAL.  Options are NORMAL and BOLD.
|           2.5.5.3: Header options
|               PARHEADER / NCENSHEADER / NEVENTSHEADER = Creates a sub-title above the corresponding table of values.
|                                                         Making this blank will remove the header.
|               PARFAMILY = Determines the font of the text in the patients-at-risk subtitle. Default is Albany AMT.  
|               PARSIZE = Determines the size of the text in the patients-at-risk subtitle. Default is 10pt.  Must contain a number and pt.
|               PARWEIGHT = Determines the font weight of the text in the patients-at-risk subtitle.  Default is normal. Options are bold and normal.
|           2.5.5.4: Line Divider options (When RISKLOCATION=INSIDE)
|               RISKDIVIDER = A flag variable to turn off the dividor line between the plot window and the patients at risk table when
|                             RISKLIST=INSIDE.  Default=1, options are 1 or 0.
|               RISKDIVCOLOR = Sets the color of the RISKDIVIDER line.  Default=black.
|               RISKDIVSTYLE = Sets the line style of the RISKDIVIDER line.  Default=solid. Values can be numbers from 1-46 or values such as: shortdash,
|                              mediumdash, longdash, dashdashdot, dash, dot, thindot.
|       2.5.6: Reference Line Options
|           REFLINES = Indicates which time-point referencelines are requested at. The options are TIMEPOINTS for
|                      the times listed in the TIMELIST parameter and MEDIANS for the median time-to-event times.
|                      The default is the null value, which will prevent any reference lines from showing.
|           REFLINEMETHOD = Sets the type of reference line.  Options are FULL for lines that go from one of the axis
|                           to the other and DROP for lines that go from the Kaplan-Meier curves to the axis.  Default is DROP.
|           REFLINEAXIS = Sets the axis that the reference lines are based off of.  The options are X, Y and Both.  Default is X.
|           REFLINESIZE = Sets the line thickness for the reference lines.  Default is 1pt.
|           REFLINEPATTERN = Sets the pattern for the reference lines.  The default is 2 (for dashed lines).
|                            Options are to do numbers between 1 and 46, or: SOLID, SHORTDASH, MEDIUMDASH, LONGDASH,
|                            MEDIUMDASHSHORTDASH, DASHDASHDOT, DASHDOTDOT, DASH, LONGDASHSHORTDASH, DOT, THINDOT, SHORTDASHDOT, and MEDIUMDASHDOTDOT.
|           REFLINECOLOR = Sets the color of the reference lines.  Default is grey.
|       2.5.7: Confidence Interval Options
|           2.5.7.1: Enable Options
|               PLOTCI = Determines if confidence bounds will be drawn.  Options are 0 (No), 1 (Yes), and 2 (Auto).  Setting 2 will
|                        enable confidence bounds on plots without a CLASS variable but not on plots with a CLASS variable.
|                        Default is 2.
|               PLOTCIFILL = Determines if a band plot will be drawn to fill in the space between the confidence intervals.  Options
|                            are 0 (No) and 1 (Yes).  Default is 1.
|           2.5.7.2: Line Options
|               PLOTCILINESIZE = Sets the thickness of the confidence interval lines.  Setting this to zero will cause the lines to not show,
|                                but the fill can still be enabled.  Default is 0pt.
|               PLOTCILINECOLOR = Sets the colors for the confidence lines.  Default is null.  Leaving this option null will cause the colors
|                                 to match the plot lines from the COLOR parameter.  Otherwise 1 color must be specified per CLASS
|                                 group (e.g. Black Red Blue).
|               PLOTCILINEPATTERN = Sets the pattern of the confidence interval lines. Default is null. Leaving this option null will cause the patterns
|                                   to match the plot lines from the PATTERN parameter.  Otherwise 1 pattern must be specified per CLASS
|                                   group (e.g. 1 2 3), or one pattern must be specified to be applied to all. Default is 2 (for dashed lines).  
|                                   Options are to do numbers between 1 and 46, or: SOLID, SHORTDASH, MEDIUMDASH, LONGDASH,
|                                   MEDIUMDASHSHORTDASH, DASHDASHDOT, DASHDOTDOT, DASH, LONGDASHSHORTDASH,
|                                   DOT, THINDOT, SHORTDASHDOT, and MEDIUMDASHDOTDOT.
|           2.5.7.3: Fill Options
|               PLOTCIFILLCOLOR = Sets the colors for the band plots.  Default is null.  Leaving this option null will cause the colors
|                                 to match the plot lines from the COLOR parameter.  Otherwise 1 color must be specified per CLASS
|                                 group (e.g. Black Red Blue).
|               PLOTCIFILLTRANSPARENCY = Sets the transparency of the band plots.  This will cause the fill to be more see-through and the graphs
|                                        to be less cluttered.  Must be a number between 0 and 1, larger numbers are more transparent.
|                                        Default is 0.95.
|       2.5.8: Grid Line Options
|               GRIDLINES = Determines if gridlines will be drawn on the x and/or y axes.  Will be drawn at tick values.  
|                           1 is yes and 0 is no.  Default is 0.
|                           NOTE: This will make the background of the STAT table to be opaque.
|               GRIDLINE_AXIS = Determines which axes have gridlines drawn. Options are X, Y, and BOTH.  Default is BOTH
|               GRIDLINE_COLOR = Determines the color of the gridlines.  Default is GREY.
|               GRIDLINE_SIZE = Determines the size/thickness of the gridlines.  Default is 1pt.
|               GRIDLINE_PATTERN = Determines the pattern of the gridlines.  Default is SOLID.
|                                  Options are to numbers between 1 and 46, or: SOLID, SHORTDASH, MEDIUMDASH, LONGDASH,
|                                  MEDIUMDASHSHORTDASH, DASHDASHDOT, DASHDOTDOT, DASH, LONGDASHSHORTDASH,
|                                  DOT, THINDOT, SHORTDASHDOT, and MEDIUMDASHDOTDOT.
|       2.5.9: Underline Options
|         DRAW_UNDERLINES = Determines if underlines are drawn under the headers and above the model summary section.  If a BY variable is present
|                           then a line is drawn for each level of the BY variable.  Options are AUTO, 0 (off), and 1 (on).  When AUTO is used then
|                           lines are drawn if a BY variable is used but not when a BY variable is absent.
|         UL_SIZE = Determines the thickness of the underlines.  Default is 1pt.
|         UL_COLOR = Determines the color of the underlines.  Default is BLACK.
|         UL_PATTERN = Determines the pattern of the underlines.  Default is SOLID.
|                      Options are to numbers between 1 and 46, or: SOLID, SHORTDASH, MEDIUMDASH, LONGDASH,
|                                  MEDIUMDASHSHORTDASH, DASHDASHDOT, DASHDOTDOT, DASH, LONGDASHSHORTDASH,
|                                  DOT, THINDOT, SHORTDASHDOT, and MEDIUMDASHDOTDOT.
|       2.5.9: Image Options
|         AXISCOLOR = Sets the color for the axes and border.  Default is black.
|         BACKGROUND = Sets the color for the background.  Ignored if TRANSPARENT=1.  Default is white.
|         BORDER = Turns the black border in the plot image on (1) or off (0).  Options are
|                  1 or 0, default is 1. 
|         FONTCOLOR = Determines the color for the graph text.  Default is black
|         SHOWWALLS = Flag indicator to turn the top and right walls of the plot border
|                     on (1) or off (0).  Default is 1, options are 0 and 1.
|         TRANSPARENT = Determines if the background will be transparent.  Only available in 9.4M3+.
|                       Default is 0 (no).  Options are 1 (yes) and 0 (no).
|       2.5.10: Lattice (Multiple graph) Options
|           COLUMNS = Sets the number of columns in a plot lattice
|           ROWS = Sets the number of rows in a plot lattice
|                  NOTE: ROWS and COLUMNS should be set to be able to contain the number of models
|                        specified by NMODELS.
|           ORDER = sets the order that plots are placed into a lattice.  Options are
|                   columnmajor and rowmajor.  Rowmajor fills rows first.  Columnmajor fills 
|                   columns first.  Default is rowmajor.
|   2.6: Table Summary (Not related to plot) Options
|       2.6.1: Display Options
|           TABLEDISPLAY = Controls the columns displayed in the table summary.  
|                Options are:
|                    TITLE: Shows the titles from the TITLE parameter (See plot section)
|                    FOOTNOTE: Shows the titles from the FOOTNOTE parameter (See plot section)
|                    TOTAL: Total sample size
|                    EVENT: Total number of events
|                    EV_N: Combined events and total in the format EVENTS/Total
|                    N_EV: Combined events and total in the format TOTAL (EVENTS)
|                    MEDIAN: Median time-to-event and confidence interval
|                    TIMELIST: Time-point event-free rates, confidence intervals, and time-points
|                    HR: Hazard ratio and confidence interval
|                    PVAL: P-value specified by the PLOTPVAL parameter
|                    COVPVAL: P-value from Cox model parameters table
|                    CINDEX: C-index from the Cox model (not available for CIF methods)
|                    RMST: Restricted Means Survival Time (Only available when METHOD=KM)
|                    RMTL: Restricted Means Time Lost (Only available when METHOD=KM)
|                    RMST_PVAL: P-value comparing RMST or RMTL values (Only available when METHOD=KM)
|                The following statistics can have the MV suffix added to get the multivariate versions
|                when CLASSCOV or CONTCOV is specified:
|                EVENTMV, EV_NMV, N_EVMV, HRMV, PVALMV, COVPVALMV, and CINDEXMV
|                The desired list of outputted columns should be listed with each parameter separated by a space.  
|                Statistics are printed in the table in the order they are listed.
|                Default=title footnote ev_n median hr timelist pval.
|           TABLESHADING = Flag variable to turn alternating shading on or off within the table
|                          summary. 1 is on and 0 is off.  Default is 1.
|           TABLETITLE = Creates a title at the top of the summary table.
|           TABLEFOOTNOTE = Creates a footnote at the end of the summary table.  Use a ` to
|                           designate a line break to create multiple footnotes.
|       2.6.2: Statistical Column Options
|           2.6.2.1: Font Options
|               TABLEDATAFAMILY = Determines the font for the area between the headers and footnotes.
|                                 Default is Arial.
|               TABLEDATASIZE = Determines the font size for the area between the headers and footnotes.
|                               Default is 9pt.
|               TABLEDATAWEIGHT = Determines the font weight for the area between the headers and footnotes.
|                                 Default is medium. Options are medium and bold.
|               TABLEFOOTNOTEFAMILY = Determines the font for the footnotes. Default is Arial.
|               TABLEFOOTNOTESIZE = Determines the font size for the footnotes. Default is 10pt.
|               TABLEFOOTNOTEWEIGHT = Determines the font weight for the footnotes.
|                                     Default is medium. Options are medium and bold.
|               TABLEHEADERFAMILY = Determines the font for the headers. Default is Arial.
|               TABLEHEADERSIZE = Determines the font size for the headers. Default is 10pt.
|               TABLEHEADERWEIGHT = Determines the font weight for the headers.
|                                   Default is bold. Options are medium and bold. 
|           2.6.2.2: Header Options
|               TprefixHEADER = The header for each item in DISPLAY can be modified by entering T, the 
|                               DISPLAY item, and then HEADER.  Options are:
|                               TTOTALHEADER, TEVENTHEADER, TEV_NHEADER, TN_EVHEADER, TMEDIANHEADER,
|                               TTIMELISTHEADER, THRHEADER, TPVALHEADER, TCOVPVALHEADER, TCINDEXHEADER,
|                               TTOTALMVHEADER, TEVENTMVHEADER, TEV_NMVHEADER, TN_EVMVHEADER, 
|                               THRMVHEADER, TPVALMVHEADER, TCOVPVALMVHEADER, TCINDEXMVHEADER,
|                               RMSTHEADER, RMTLHEADER, RMST_PVALHEADER
|                               NOTE: The ^n character acts as a line break.
|           2.6.2.3: Width Options
|               TprefixWIDTH = The column width for each item in DISPLAY can be modified by entering T, the 
|                               DISPLAY item, and then WIDTH.  Options are:
|                               TTOTALWIDTH, TEVENTWIDTH, TEV_NWIDTH, TN_EVWIDTH, TMEDIANWIDTH,
|                               TTIMELISTWIDTH, THRWIDTH, TPVALWIDTH, TCOVPVALWIDTH, TCINDEXWIDTH,
|                               TTOTALMVWIDTH, TEVENTMVWIDTH, TEV_NMVWIDTH, TN_EVMVWIDTH, 
|                               THRMVWIDTH, TPVALMVWIDTH, TCOVPVALMVWIDTH, TCINDEXMVWIDTH,
|                               RMSTWIDTH, RMTLWIDTH, RMST_PVALWIDTH
|                               TSUBTITLEWIDTH controls the width of the first column that contains
|                               the CLASS variable levels.
|                               NOTE: Widths should be entered with a unit of measure such as in (inches).
|   2.6: Macro Debugging Options
|       DEBUG = Determines if temporary datasets are deleted and notes are turned off.  Options are:
|               0: No notes are displayed and temporary datasets are deleted (Default)
|               1: Notes are displayed and temporary datasets are left behind.  MPRINT is also turned on.
*------------------------------------------------------------------*
| OPERATING SYSTEM COMPATIBILITY
| UNIX SAS v9.4M3   :   YES
| PC SAS v9.4M3     :   YES
*------------------------------------------------------------------*
| MACRO CALL
|
| %newsurv (
|            DATA=,
|            TIME=,
|            CENS=,
|            CEN_VL=
|          );
*------------------------------------------------------------------*
| REQUIRED PARAMETERS
|
| Name      : DATA
| Default   : 
| Type      : Dataset Name
| Purpose   : REFER TO REFERENCE SECTION
|
| Name      : TIME
| Default   :
| Type      : Variable Name
| Purpose   : REFER TO REFERENCE SECTION
|
| Name      : CENS
| Default   :
| Type      : Variable Name
| Purpose   : REFER TO REFERENCE SECTION
|
| Name      : CEN_VL
| Default   :
| Type      : Variable Name
| Purpose   : REFER TO REFERENCE SECTION
|
| Name      : METHOD
| Default   :
| Type      : Variable Name
| Purpose   : REFER TO REFERENCE SECTION
|
| Name      : EV_VL
| Default   :
| Type      : Variable Name
| Purpose   : REFER TO REFERENCE SECTION
|
*------------------------------------------------------------------*
| EXAMPLES (MUST BE RUN IN 9.3+ AS SASHELP.BMT DOES NOT EXIST IN 9.2)
|
| Example 1: Basic Example Call:
| %newsurv(data=sashelp.bmt, time=T, cens=STATUS, cen_vl=0);
| Example 2: Call with Class Variable:
| %newsurv(data=sashelp.bmt, time=T, cens=STATUS, cen_vl=0,class=Group);
| Example 3: Call with Class Variable and more options:
| %newsurv(data=sashelp.bmt, time=T, cens=STATUS, cen_vl=0,class=group,
    classref=ALL,summary=0, plot=1,xmin=0, xmax=2500,ptabsize=7pt,
    xincrement=500,symbolsize=5pt,
    xlabel=Time (Days), ylabel=Proportion Alive,autoalign=topright,
    title=Example 3,ytype=ppt,color=BLACK BLUE RED);
| Example 4: Call with Patients-At-Risk table within Plot Window:
| %newsurv(data=sashelp.bmt, time=T, cens=STATUS, cen_vl=0,class=group,
    summary=0, risklist=0 to 2500 by 500, risklocation=INSIDE,
    xmin=0, xmax=2500, xincrement=500,outp=test);
| Example 5: Call with Patients-At-Risk table below Plot Window:
| %newsurv(data=sashelp.bmt, time=T, cens=STATUS, cen_vl=0,class=group,
    summary=0, risklist=0 to 2500 by 500, risklocation=BOTTOM,
    xmin=0, xmax=2500, xincrement=500,risklabellocation=above);
| Example 6: Call with Kaplan-Meier Time-Point Estimates:
| %newsurv(data=sashelp.bmt, time=T, cens=STATUS, cen_vl=0,class=group,
    summary=0, timelist=1000 2000, timedx=Days);
| Example 7: Call with a basic lattice set-up:
| %newsurv(data=sashelp.bmt, time=T, cens=STATUS, cen_vl=0,class=|group,
    nmodels=2,summary=0,height=8in,rows=2,xdivisor=365.25,
    title=Overall Survival|Overall Survival by Disease Group);
| Example 8: Call with confidence intervals:
| %newsurv(data=sashelp.bmt, time=T, cens=STATUS, cen_vl=0,class=group,
    summary=0,height=4in,xdivisor=365.25,plotci=1,color=black red blue);
| Example 9: Call with reference lines:
| %newsurv(data=sashelp.bmt, time=T, cens=STATUS, cen_vl=0,class=group,
    summary=0,height=4in,xdivisor=365.25,reflines=medians,reflineaxis=both);

|**Example dataset for CIF Examples:
proc format;
   value grpLabel 1='ALL' 2='AML low risk' 3='AML high risk';
run;

data BMT;
        input DIAGNOSIS Ftime Status Gender@@;
        label Ftime="Days";
        format Diagnosis grpLabel.;
datalines;
1       2081       0       1       1       1602    0       1
1       1496       0       1       1       1462    0       0
1       1433       0       1       1       1377    0       1
1       1330       0       1       1       996     0       1
1       226        0       0       1       1199    0       1
1       1111       0       1       1       530     0       1
1       1182       0       0       1       1167    0       0
1       418        2       1       1       383     1       1
1       276        2       0       1       104     1       1
1       609        1       1       1       172     2       0
1       487        2       1       1       662     1       1
1       194        2       0       1       230     1       0
1       526        2       1       1       122     2       1
1       129        1       0       1       74      1       1
1       122        1       0       1       86      2       1
1       466        2       1       1       192     1       1
1       109        1       1       1       55      1       0
1       1          2       1       1       107     2       1
1       110        1       0       1       332     2       1
2       2569       0       1       2       2506    0       1
2       2409       0       1       2       2218    0       1
2       1857       0       0       2       1829    0       1
2       1562       0       1       2       1470    0       1
2       1363       0       1       2       1030    0       0
2       860        0       0       2       1258    0       0
2       2246       0       0       2       1870    0       0
2       1799       0       1       2       1709    0       0
2       1674       0       1       2       1568    0       1
2       1527       0       0       2       1324    0       1
2       957        0       1       2       932     0       0
2       847        0       1       2       848     0       1
2       1850       0       0       2       1843    0       0
2       1535       0       0       2       1447    0       0
2       1384       0       0       2       414     2       1
2       2204       2       0       2       1063    2       1
2       481        2       1       2       105     2       1
2       641        2       1       2       390     2       1
2       288        2       1       2       421     1       1
2       79         2       0       2       748     1       1
2       486        1       0       2       48      2       0
2       272        1       0       2       1074    2       1
2       381        1       0       2       10      2       1
2       53         2       0       2       80      2       0
2       35         2       0       2       248     1       1
2       704        2       0       2       211     1       1
2       219        1       1       2       606     1       1
3       2640       0       1       3       2430    0       1
3       2252       0       1       3       2140    0       1
3       2133       0       0       3       1238    0       1
3       1631       0       1       3       2024    0       0
3       1345       0       1       3       1136    0       1
3       845        0       0       3       422     1       0
3       162        2       1       3       84      1       0
3       100        1       1       3       2       2       1
3       47         1       1       3       242     1       1
3       456        1       1       3       268     1       0
3       318        2       0       3       32      1       1
3       467        1       0       3       47      1       1
3       390        1       1       3       183     2       0
3       105        2       1       3       115     1       0
3       164        2       0       3       93      1       0
3       120        1       0       3       80      2       1
3       677        2       1       3       64      1       0
3       168        2       0       3       74      2       0
3       16         2       0       3       157     1       0
3       625        1       0       3       48      1       0
3       273        1       1       3       63      2       1
3       76         1       1       3       113     1       0
3       363        2       1
;
run;
| Example 10: Example using the CIF method (run above data step first)
| %newsurv(data=bmt, time=ftime, cens=STATUS, cen_vl=0,ev_vl=1,
    method=cif,class=diagnosis,height=4in,
    summary=0, risklist=0 to 2500 by 500, risklocation=BOTTOM,
    xmin=0, xmax=2500, xincrement=500,risklabellocation=left,
    timelist=250 500, timedx=Days,parheader=);
*------------------------------------------------------------------*
| REFERENCES
| The code for the CIF method is transcribed from the SAS autocall
|   macro %CIF.  The code was originally written primarily in IML,
|   so within this macro it was rewritten in data step language
|   instead.  The references within the %CIF macro are as follows:
|   1. Marubini, E. and Valsecchi, M.G. (1995), Analysing survival data from clinical trials
|      and observational studies, John Wiley.
|   2. Gray, R.J. (1988). "A class of K-sample tests for comparing the cumulative incidence of
|      a competing risk," Annals of statistics, 16(3), 1141--1154.
|   3. Klein, J.P. and Moeschberger, M.L., (2003), Survival analysis: techniques for censored
|      and truncated data, Springer Verlag.
| Methods for calculating concordance index and standard error in Cox 
| proportional hazards regression taken from:
| Therneau T (2014). _A Package for Survival Analysis in S_. R package
| version 2.37-7, <URL: http://CRAN.R-project.org/package=survival>.
|
| Terry M. Therneau and Patricia M. Grambsch (2000). _Modeling Survival
| Data: Extending the Cox Model_. Springer, New York. ISBN
| 0-387-98784-3.
|
| Therneau TM, Crowson CS, Atkinson EJ. 2015. Adjusted Survival Curves. 
| https://cran.r-project.org/web/packages/survival/vignettes/adjcurve.pdf. 
| Accessed October 31, 2015.
*------------------------------------------------------------------*
| This program is free software; you can redistribute it and/or
| modify it under the terms of the GNU General Public License as
| published by the Free Software Foundation; either version 2 of
| the License, or (at your option) any later version.
|
| This program is distributed in the hope that it will be useful,
| but WITHOUT ANY WARRANTY; without even the implied warranty of
| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
| General Public License for more details.
*------------------------------------------------------------------*/


%macro newsurv(
    /*** 1.0: Required Variables ***/
    cens=,cen_vl=0,data=,ev_vl=,method=KM, time=,
    
    /*** 2.0: Optional Variables ***/
    /** 2.1: Global Options **/
    /* 2.1.1: Output Controlling Options*/
    border=1,destination=rtf,newtable=1,out=,outdoc=,
    outp=,plot=1,summary=1,plot_unadjust=1,
    /* 2.1.2: Image Controlling Options*/
    antialiasmax=1000,axiscolor=black,background=white,dpi=200,
    fontcolor=black,gpath=,height=6in,plotname=_surv,plottype=png,
    svg=0,tiffdevice=TIFFP,transparent=0,width=8in, 
    /* 2.1.3: Lattice Controlling Options*/
    columns=1,nmodels=1,order=columnmajor,rows=1,rowgutter=0,columngutter=0,
    
    /** 2.2: Class Variables **/
    class=,classcov=,classdesc=,classorder=,
    classref=,contcov=,desc=0,hrties=BRESLOW,plotpval=,
    refhrtext=Reference,refptext=--,strata=,plotpvalmv=,
    /** 2.2.1: By Variables (Subset CLASS Variables) **/
    by=,byorder=,bylabel=,pval_inter=wald,pval_intermv=wald,
    /** 2.3: Dataset Modifiers **/
    landmark=,xdivisor=, where=,
    
    /** 2.4: Patient Kaplan-Meier Options **/
        /**2.4.1: Time-point Estimate Options **/
        alpha=0.05,cifvar=COUNT,conftype=LOG,timelist=,    
        /**2.4.2: Restricted Means Survival Time Options **/
        tau=,rmst_ci=1,
    /** 2.5: Plot Options **/ 
    /* 2.5.1: Axis Options*/
    lfamily=Arial,lsize=10pt,lweight=bold,showwalls=1,
    xincrement=,xlabel=,xmax=,xmaxoffset=,xmin=0,xminoffset=,
    xtickvalfamily=Arial,xtickvalsize=8pt,xtickvalweight=normal,
    yincrement=,ylabel=,ymax=,ymaxoffset=,ymin=0,yminoffset=0,
    ytickvalfamily=Arial,ytickvalsize=8pt,ytickvalweight=normal,ytype=pct,
    /* 2.5.2: Plot Statistal Table Options*/
    autoalign=topright bottomleft,classvalalign=center,
    cindexheader=,cindexmvheader=,
    covpvalheader=Wald P-value,covpvalmvheader=Adj Wald P-value,
    display=standard, model_stats_display=1,
    
    legendlinelength=,
    ev_nheader=Events/Total,eventheader=Event,hrheader=/*HR (95% CI)*/,KMEstheader=/*KM Est (95% CI)*/,
    legendheader=%str( ),location=inside,medianheader=/*Median (95% CI)*/,
    rmstheader=/*RMST (95% CI)*/,rmtlheader=/*RMTL (95% CI)*/,rmst_pvalheader=,n_evheader=Total (Events),
    totalmvheader=MV Total,eventmvheader=MV Event,hrmvheader=/*Adj HR (95% CI)*/,
    ev_nmvheader=MV Events/Total,n_evmvheader=MV Total (Events),pvalmvheader=,pval_interheader=,pval_intermvheader=,
    ptabsize=8pt, ptabfamily=Albany AMT,pvalheader=,risktableheader=N at Risk,statcolor=0,tablecomments=,
    timedx=,listtimepoints=1,timelistheader=Time-Point,totalheader=Total,       
    /* 2.5.3: Patients-at-Risk Options*/
    paralign=CENTER,parfamily=Albany AMT,parheader=Patients-at-Risk,parsize=10pt,parweight=normal,
    ncensheader=No. Cumulative Censors,neventsheader=No. Cumulative Events,pardisplay=par,
    riskcolor=0,riskdivcolor=black,riskdivider=1,riskdivstyle=solid,
    risklabelalign=LEFT,risklabeldlm=,risklabellocation=LEFT,risklabelweight=normal,
    risklist=,risklocation=bottom,riskrowweights=0.025,uniformheight=0,
    /* 2.5.4: Plot Lines/Symbols Options*/
    censormarkers=1,color=black,linesize=1pt,pattern=AUTO,pattern_adjust=auto,
    sreverse=0, symbolsize=3pt,symbolweight=normal,
    /*2.5.5: Title/Footnote Options*/
    fnfamily=Albany AMT,fnsize=8pt,fnweight=normal,
    footnote=, footnotealign=left,
    ovfnfamily=Albany AMT,ovfnsize=8pt,ovfnweight=normal,
    ovfootnote=,ovfootnotealign=left,
    ovtfamily=Albany AMT,ovtitle=,ovtitlealign=center,ovtsize=12pt,ovtweight=bold, 
    tfamily=Albany AMT,title=, titlealign=center,
    tsize=12pt,tweight=bold,
    /* 2.5.4: Reference Line Options*/
    reflines=,reflinesize=1pt,reflinepattern=2,reflinecolor=grey,reflinemethod=drop,reflineaxis=X,
    /* 2.5.5: Confidence Interval Options*/
    plotci=2,plotcifill=1,plotcifillcolor=,plotcifilltransparency=0.95,
    plotcilinecolor=,plotcilinesize=1pt,plotcilinepattern=2,
    /* Grid line Options*/
    gridlines=0, gridline_axis=both, gridline_color=grey, gridline_size=1pt,gridline_pattern=solid,
    /* Underlines*/
    draw_underlines=AUTO,ul_color=black,ul_pattern=solid,ul_size=1pt,
    /** 2.6: Optional Table Summary Options **/ 
    /* 2.6.1: Display Options */
    tablemergepval=0,tablefootnote=,
    tabledisplay=standard, 
    tableshading=1,tabletitle=,
    /* 2.6.2: Font Options */
    tabledatafamily=Arial,tabledatasize=9pt,tabledataweight=medium,
    tablefootnotefamily=Arial,tablefootnotesize=10pt,tablefootnoteweight=medium,
    tableheaderfamily=Arial,tableheadersize=10pt,tableheaderweight=bold,
    /* 2.6.3: Column Heading Options*/
    tcovpvalheader=Covariate Level~P-values,tcovpvalmvheader=Adjusted~Covariate Level~P-values,
    teventheader=Event,teventmvheader=MV Event,tev_nheader=Event/Total,tev_nmvheader=MV Event/Total,
    tn_evheader=Total (Events),tn_evmvheader=MV Total (Events),
    thrheader=,thrmvheader=,
    ttimelistheader=,
    ttimelistmvheader=,
    tmedianheader=,tmedianmvheader=,trmstheader=,trmtlheader=,trmst_pvalheader=RMST~P-value,
    tpvalheader=P-value,tpvalmvheader=Adjusted~P-value,
    tpval_interheader=Interaction~P-value,tpval_intermvheader=Adjusted~Interaction~P-value,
    ttotalheader=Total,ttotalmvheader=MV Total, 
    tcindexheader=,tcindexmvheader=,
    /* 2.6.4: Column Width Options */
    tsubtitlewidth=2in,
    tcovpvalwidth=0.7in,tcovpvalmvwidth=0.7in,
    ttotalwidth=0.5in,teventwidth=0.5in,tev_nwidth=1in,tn_evwidth=1in,
    ttotalmvwidth=0.5in,teventmvwidth=0.5in,tev_nmvwidth=1in,tn_evmvwidth=1in,
    thrwidth=1.1in,thrmvwidth=1.1in,tmedianwidth=1.3in,tmedianmvwidth=1.3in,
    trmstwidth=1.3in,trmtlwidth=1.3in,trmst_pvalwidth=0.7in,
    tpvalwidth=0.7in,tpvalmvwidth=0.7in,tpval_interwidth=0.7in,tpval_intermvwidth=0.7in,
    ttimelistwidth=1.6in,ttimelistmvwidth=1.6in,
    tcindexwidth=1.1in,tcindexmvwidth=1.1in,
    /* 2.7: Statistic Significant Digits Options */
    pvaldigits=4,hrdigits=2,mediandigits=1,kmestdigits=auto,cindexdigits=2,rmstdigits=1,
    /* 2.8: Debugging Options */
    debug=0);
    
    


    /**Save current options to reset after macro runs**/
    %local _mergenoby _notes _qlm _odspath _starttime _device _gsfname _mprint
        _xmax _ymax _xpixels _ypixels _imagestyle _iback _listing _linesize _center _gborder _msglevel;
    %let _starttime=%sysfunc(time());
    %let _notes=%sysfunc(getoption(notes));
    %let _mergenoby=%sysfunc(getoption(mergenoby));
    %let _qlm=%sysfunc(getoption(quotelenmax)); 
    %let _device=%sysfunc(getoption(device));
    %let _gsfname=%sysfunc(getoption(gsfname));
    %let _xmax=%sysfunc(getoption(xmax));
    %let _ymax=%sysfunc(getoption(ymax));
    %let _xpixels=%sysfunc(getoption(xpixels));
    %let _ypixels=%sysfunc(getoption(ypixels));
    %let _imagestyle=%sysfunc(getoption(imagestyle));
    %let _gborder=%sysfunc(getoption(border));
    %let _iback=%sysfunc(getoption(iback));
    %let _linesize=%sysfunc(getoption(linesize));
    %let _center=%sysfunc(getoption(center));
    %let _msglevel=%sysfunc(getoption(msglevel));
    %let _mprint=%sysfunc(getoption(mprint));
    %let _odspath=&sysodspath;
    %if %sysevalf(%superq(_odspath)=,boolean) %then %let _odspath=WORK.TEMPLAT(UPDATE) SASHELP.TMPLMST (READ);
    /**Turn off warnings for merging without a by and long quote lengths**/
    /**Turn off notes**/
    options mergenoby=NOWARN nonotes noquotelenmax msglevel=N;
    ods path WORK.TEMPLAT(UPDATE) SASHELP.TMPLMST (READ);
    
    /*Don't send anything to output window, results window, and set escape character*/
    ods select none;
    ods noresults escapechar='^';
    
    /**Process Error Handling**/
    %if &sysvlong < 9.04.01M3P062415 %then %do;
        %put ERROR: SAS must be version 9.4M3 or later;
        %goto errhandl;
    %end;       
    %else %if %sysfunc(exist(&data))=0 %then %do;
        %put ERROR: Dataset &data does not exist;
        %put ERROR: Please enter a valid dataset;
        %goto errhandl;
    %end;
    %else %if %sysevalf(%superq(data)=,boolean)=1 %then %do;
        %put ERROR: DATA parameter is required;
        %put ERROR: Please enter a valid dataset;
        %goto errhandl;
    %end;
    
    /**Pull dataset information**/
    /*proc contents data=&data out=_temp noprint;
    run;

    /**See if the listing output is turned on**/
    proc sql noprint;
        select 1 into :_listing separated by '' from sashelp.vdest where upcase(destination)='LISTING';
    quit;
    /**Create list of macro variables that can vary across different models called**/
    /**Sets up for lattice plots**/
    %local _mvarlist;
    %let _mvarlist=%sysfunc(compress(alpha|autoalign|
        cens|cen_vl|censormarkers|class|classcov|classdesc|
        classorder|classref|classvalalign|cifvar|color|conftype|contcov|covpvalheader|covpvalmvheader|
        desc|display|ev_nheader|ev_nmvheader|ev_vl|eventheader|eventmvheader|fnfamily|fnsize|fnweight|footnote|
        footnotealign|hrheader|hrmvheader|hrties|KMEstheader|RMSTheader|RMTLheader|RMST_PVALheader|
        landmark|legendheader|legendlinelength|lfamily|linesize|listtimepoints|location|lsize|lweight|medianheader|
        method|n_evheader|n_evmvheader|paralign|parfamily|parheader|parsize|parweight|
        pattern|plotci|plotcifill|plotcifillcolor|plotcifilltransparency|
        plotcilinecolor|plotcilinesize|plotcilinepattern|    
        plotpval|plotpvalmv|pvalheader|pvalmvheader|ptabfamily|ptabsize|refhrtext|refptext|
        riskcolor|riskdivcolor|riskdivider|riskdivstyle|risklabelalign|risklabeldlm|risklabellocation|
        risklabelweight|risklist|risklocation|riskrowweights|risktableheader|sreverse|statcolor|strata|symbolsize|symbolweight|tablecomments|
        tfamily|time|timedx|timelist|timelistheader|title|titlealign|totalheader|totalmvheader|tsize|tweight|
        where|xdivisor|xincrement|xlabel|xmax|xmaxoffset|xmin|xminoffset|
        reflinecolor|reflinemethod|reflinepattern|reflines|reflinesize|reflineaxis|
        xtickvalfamily|xtickvalsize|xtickvalweight|
        ylabel|yincrement|ymax|ymaxoffset|ymin|yminoffset|ytickvalfamily|ytickvalsize|ytickvalweight|ytype|
        pvaldigits|hrdigits|mediandigits|kmestdigits|rmstdigits|ncensheader|neventsheader|pardisplay|
        cindexheader|cindexmvheader|cindexdigits|gridlines|gridline_axis|gridline_color|gridline_size|gridline_pattern|plot_unadjust|pattern_adjust|
        tau|rmst_ci|
        by|byorder|bylabel|pval_inter|pval_intermv|pval_interheader|pval_intermvheader|model_stats_display|draw_underlines|ul_color|ul_pattern|ul_size));
        
    %local i j;
    %do i = 1 %to &nmodels;
        /**Cycle through each macro parameter**/
        %do j = 1 %to %sysfunc(countw(&_mvarlist,|));
            %local v&j;
            %let v&j=%scan(%superq(_mvarlist),&j,|);
            %local &&v&j..&i;
            /**If the | delimiter is detected, assign the different values between | to numbered parameters**/
            /**Else Assign the same value to all numbered parameters**/
            %if %index(%superq(&&v&j),|)>0 %then %let &&v&j..&i=%qscan(%superq(&&v&j),&i,|,m);
            %else %let &&v&j..&i=%qscan(%superq(&&v&j),1,|,m); 
        %end;   
    %end;                       
    %local z nerror;
    %let nerror=0;
    /**Error Handling on Individual Model Variables**/
    %macro _varcheck(var,require,numeric);
        %local _z _numcheck;
        %do z = 1 %to &nmodels;
            /**Check if variable parameter is missing**/
            %if %sysevalf(%superq(&var.&z)=,boolean)=0 %then %do;
                %if %sysfunc(notdigit(%superq(&var.&z))) > 0 %then
                    %do _z = 1 %to %sysfunc(countw(%superq(&var.&z),%str( )));
                    /**Check to make sure variable names are not just numbers**/    
                    %local datid;
                    /**Open up dataset to check for variables**/
                    %let datid = %sysfunc(open(&data));
                    /**Check if variable exists in dataset**/
                    %if %sysfunc(varnum(&datid,%scan(%superq(&var.&z),&_z,%str( )))) = 0 %then %do;
                        %put ERROR: (Model &z: %qupcase(&var)) Variable %qupcase(%scan(%superq(&var.&z),&_z,%str( ))) does not exist in dataset &data;
                        %local closedatid;
                        /**Close dataset**/
                        %let closedatid=%sysfunc(close(&datid));
                        %let nerror=%eval(&nerror+1);
                    %end;
                    %else %do;
                        %local closedatid;
                        %let closedatid=%sysfunc(close(&datid));
                        %if &numeric=1 %then %do;
                            data _null_;
                                set &data (obs=1);
                                call symput('_numcheck',strip(vtype(%superq(&var.&z))));
                            run;
                            %if %sysevalf(%superq(_numcheck)^=N,boolean) %then %do;
                                %put ERROR: (Model &z: %qupcase(&var)) Variable must be numeric;
                                %let nerror=%eval(&nerror+1);
                            %end;   
                        %end;                         
                    %end;
                %end;
                %else %do;
                    /**Give error message if variable name is number**/
                    %put ERROR: (Model &z: %qupcase(&var)) Variable is not a valid SAS variable name (%superq(&var.&z));
                    %let nerror=%eval(&nerror+1);
                %end;
            %end;
            %else %if &require=1 %then %do;
                /**Give error if required variable is missing**/
                %put ERROR: (Model &z: %qupcase(&var)) Variable is a required variable but has no value;
                %let nerror=%eval(&nerror+1);
            %end;
        %end;
    %mend;
    /**Check time variables**/
    %_varcheck(time,1,1)
    /**Censor Variables**/
    %_varcheck(cens,1,1)
    /**Class Variables**/
    %_varcheck(class,0)
    /**BY Variables**/
    %_varcheck(by,0)
    %local z;
    %do z = 1 %to &nmodels;
        %if (%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT) and 
            (%sysevalf(%superq(class&z)=,boolean) or (%sysevalf(%superq(classcov&z)=,boolean) and %sysevalf(%superq(contcov&z)=,boolean))) %then %do;
            /**Give error message**/
            %if %sysevalf(%superq(class&z)=,boolean) %then %put ERROR: (Model &z: %qupcase(METHOD)) CLASS cannot be missing when METHOD=%qupcase(&&method&z);
            %if %sysevalf(%superq(classcov&z)=,boolean) and %sysevalf(%superq(contcov&z)=,boolean) %then 
                %put ERROR: (Model &z: %qupcase(METHOD)) CLASSCOV and CONTCOV cannot both be missing when METHOD=%qupcase(&&method&z);
            %let nerror=%eval(&nerror+1);
        %end;
    %end;        
    /**Strata Variables**/
    %_varcheck(strata,0)
    /**Class Type Covariate Variables**/
    %_varcheck(classcov,0)
    /**Continuous Type Covariate Variables**/
    %_varcheck(contcov,0)
    /**Landmark Variables**/
    %do z = 1 %to &nmodels;
        /**Check if variable parameter is missing**/
        %if %sysevalf(%superq(landmark&z)=,boolean)=0 %then %do;
            %if %sysfunc(notdigit(%sysfunc(compress(%superq(landmark&z),.-)))) > 0 %then %do;
                /*Check if number is first digit*/
                %if %sysfunc(anydigit(%sysfunc(compress(%superq(landmark&z),.-)))) = 1 %then %do;
                    %put ERROR: (Model &z: %qupcase(landmark)) Must be a valid variable name or a number greater than 0. %qupcase(%superq(landmark&z)) is not valid.;
                    %let nerror=%eval(&nerror+1);
                %end;
                %else %do;
                    /**Check to make sure variable names are not just numbers**/
                    %local datid;
                    /**Open up dataset to check for variables**/
                    %let datid = %sysfunc(open(&data));
                    /**Check if variable exists in dataset**/
                    %if %sysfunc(varnum(&datid,%superq(landmark&z))) = 0 %then %do;
                        %put ERROR: (Model &z: %qupcase(landmark)) Variable %qupcase(%superq(landmark&z)) does not exist in dataset &data;
                        %local closedatid;
                        /**Close dataset**/
                        %let closedatid=%sysfunc(close(&datid));
                        %let nerror=%eval(&nerror+1);
                    %end;
                    %else %do;
                        %local closedatid;
                        %let closedatid=%sysfunc(close(&datid));
                    %end;
                %end;
            %end;
            %else %if %superq(landmark&z) lt 0 %then %do;
                /**Check if value is below minimum threshold**/
                %put ERROR: (Model &z: %qupcase(landmark)) Must be greater than 0. %qupcase(%superq(landmark&z)) is not valid.;
                %let nerror=%eval(&nerror+1);
            %end;  
        %end;
    %end;
    
             
    /**Error Handling on Individual Model Parameters Involving units**/
    %macro _unitcheck(parm,allowmissing);
        %do z = 1 %to &nmodels;
            %if %sysevalf(%superq(&parm.&z)=,boolean)=1 %then %do;
                %if %sysevalf(&allowmissing^=1,boolean) %then %do;
                    /**Check for missingness**/
                    %put ERROR: (Model &z: %qupcase(&parm)) Cannot be set to missing;
                    %let nerror=%eval(&nerror+1);
                 %end;
            %end;
            %else %if %sysfunc(compress(%superq(&parm.&z),ABCDEFGHIJKLMNOPQRSTUVWXYZ,i)) lt 0 %then %do;
                /**Check if value is less than zero**/
                %put ERROR: (Model &z: %qupcase(&parm)) Cannot be less than zero (%qupcase(%superq(&parm.&z)));
                %let nerror=%eval(&nerror+1);
            %end;
        %end;
    %mend;
    /**Font Sizes**/
    /*Y Tick Value Font Size**/
    %_unitcheck(ytickvalsize)
    /**X Tick Value Font Size**/
    %_unitcheck(xtickvalsize)
    /***Label Font Size**/
    %_unitcheck(lsize)
    /**Plot Table Font Size**/
    %_unitcheck(ptabsize)
    /**Title Font Size**/
    %_unitcheck(tsize)
    /**Footnote Font Size**/
    %_unitcheck(fnsize)
    /***Patients-at-Risk Header Font Size**/
    %_unitcheck(parsize)
    /**Plot Line Size**/
    %_unitcheck(linesize)
    /**Grid Line Size**/
    %_unitcheck(gridline_size)
    /**Plot Symbol Size**/
    %_unitcheck(symbolsize)
    /**Plot Confidence Bounds Line Size**/
    %_unitcheck(plotcilinesize)
    /**Plot Reference Line Size**/
    %_unitcheck(reflinesize)
    /**Plot Reference Line Size**/
    %_unitcheck(legendlinelength,1)
    /**Underline Size**/
    %_unitcheck(ul_size)
    /**Error Handling on Individual Model Numeric Variables**/
    %macro _numcheck(parm,min,contain,default,max);
        %do z = 1 %to &nmodels;
            /**Check if missing**/
            %if %sysevalf(%superq(&parm.&z)=,boolean)=0 %then %do;
                %if %sysfunc(notdigit(%sysfunc(compress(%superq(&parm.&z),-.)))) > 0 %then %do;
                    /**Check if character values are present**/
                    %put ERROR: (Model &z: %qupcase(&parm)) Must be numeric.  %qupcase(%superq(&parm.&z)) is not valid.;
                    %let nerror=%eval(&nerror+1);
                %end;  
                %else %if %superq(&parm.&z) le &min and &contain=0 %then %do;
                    /**Check if value is below minimum threshold**/
                    %put ERROR: (Model &z: %qupcase(&parm)) Must be greater than &min..  %qupcase(%superq(&parm.&z)) is not valid.;
                    %let nerror=%eval(&nerror+1);
                %end;  
                %else %if %superq(&parm.&z) lt &min and &contain=1 %then %do;
                    /**Check if value is below minimum threshold**/
                    %put ERROR: (Model &z: %qupcase(&parm)) Must be greater than or equal to &min..  %qupcase(%superq(&parm.&z)) is not valid.;
                    %let nerror=%eval(&nerror+1);
                %end; 
                %else %if %sysevalf(%superq(max)^=,boolean) %then %do;
                    %if %superq(&parm.&z) gt &max and &contain=1 %then %do;
                        /**Check if value is above maximum threshold**/
                        %put ERROR: (Model &z: %qupcase(&parm)) Must be less than or equal to &max..  %qupcase(%superq(&parm.&z)) is not valid.;
                        %let nerror=%eval(&nerror+1);
                    %end; 
                    %else %if %superq(&parm.&z) ge &max %then %do;
                        /**Check if value is above maximum threshold**/
                        %put ERROR: (Model &z: %qupcase(&parm)) Must be less than &max..  %qupcase(%superq(&parm.&z)) is not valid.;
                        %let nerror=%eval(&nerror+1);
                    %end; 
                %end;
            %end;   
            %else %let &parm.&z=&default;        
        %end;
    %mend;
    /**X Axis Minimum Value**/
    %_numcheck(xmin,0,1,0)
    /**Y Axis Minimum Value**/
    %_numcheck(ymin,0,1,0) 
    /**Alpha Value**/
    %_numcheck(alpha,0,0,0.05,1) 
    /**RMST Tau**/
    %_numcheck(tau,0,0,) 
    %local _test j;
    %let _test=1;
    %if &nerror=0 %then %do i = 1 %to &nmodels;
        /**Check alpha titles**/
        %if %sysevalf(%superq(hrheader&i)=,boolean) %then %let hrheader&i=HR (%sysevalf(100-&&alpha&i*100)% CI);
        %if %sysevalf(%superq(kmestheader&i)=,boolean) %then %let kmestheader&i=KM Est (%sysevalf(100-&&alpha&i*100)% CI);
        %if %sysevalf(%superq(medianheader&i)=,boolean) %then %let medianheader&i=Median (%sysevalf(100-&&alpha&i*100)% CI);
        %if %sysevalf(%superq(hrmvheader&i)=,boolean) %then %let hrmvheader&i=Adj HR (%sysevalf(100-&&alpha&i*100)% CI);
        %if %sysevalf(%superq(rmstheader&i)=,boolean) %then %let rmstheader&i=RMST (%sysevalf(100-&&alpha&i*100)% CI);
        %if %sysevalf(%superq(rmtlheader&i)=,boolean) %then %let rmtlheader&i=RMTL (%sysevalf(100-&&alpha&i*100)% CI);
        %if &i>1 %then %do;
            %let j=%sysevalf(&i-1);
            %if &&alpha&i^=&&alpha&j %then %let _test=0;
        %end;
        %if &_test =0 and &summary=1 %then %do;
            %put WARNING: ALPHA is not consistent between models.  Default column headers containing confidence intervals may not have accurate range;
            %if %sysevalf(%superq(thrheader)=,boolean) %then %let thrheader=Hazard Ratio~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(thrmvheader)=,boolean) %then %let thrmvheader=Adjusted~Hazard Ratio~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(ttimelistheader)=,boolean) %then %let ttimelistheader=Survival Estimates~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(ttimelistmvheader)=,boolean) %then %let ttimelistmvheader=Adjusted~Survival Estimates~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(tmedianheader)=,boolean) %then %let tmedianheader=Median~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(tmedianmvheader)=,boolean) %then %let tmedianmvheader=Adjusted~Median~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(tcindexheader)=,boolean) %then %let tcindexheader=C-index~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(tcindexmvheader)=,boolean) %then %let tcindexmvheader=Adjusted~C-index~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(trmstheader)=,boolean) %then %let trmstheader=RMST (%sysevalf(100-&&alpha&i*100)% CI);
            %if %sysevalf(%superq(trmtlheader)=,boolean) %then %let trmtlheader=RMTL (%sysevalf(100-&&alpha&i*100)% CI);
        %end;
        %else %do;
            %if %sysevalf(%superq(thrheader)=,boolean) %then %let thrheader=Hazard Ratio~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(thrmvheader)=,boolean) %then %let thrmvheader=Adjusted~Hazard Ratio~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(ttimelistheader)=,boolean) %then %let ttimelistheader=Survival Estimates~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(ttimelistmvheader)=,boolean) %then %let ttimelistmvheader=Adjusted~Survival Estimates~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(tmedianheader)=,boolean) %then %let tmedianheader=Median~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(tmedianmvheader)=,boolean) %then %let tmedianmvheader=Adjusted~Median~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(tcindexheader)=,boolean) %then %let tcindexheader=C-index~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(tcindexmvheader)=,boolean) %then %let tcindexmvheader=Multivariate~C-index~(%sysevalf(100-&alpha1*100)% CI);
            %if %sysevalf(%superq(trmstheader)=,boolean) %then %let trmstheader=RMST (%sysevalf(100-&&alpha&i*100)% CI);
            %if %sysevalf(%superq(trmtlheader)=,boolean) %then %let trmtlheader=RMTL (%sysevalf(100-&&alpha&i*100)% CI);
        %end;
    %end;
    /**Censor Value**/
    %do z = 1 %to &nmodels;
        /**Check if missing**/
        %if %sysevalf(%superq(cen_vl&z)=,boolean)=0 %then %do;
            %if %sysfunc(notdigit(%sysfunc(compress(%superq(CEN_VL&z),.)))) > 0 %then %do;
                /**Check if character values are present**/
                %put ERROR: (Model &z: %qupcase(CEN_VL)) Must be numeric.  %qupcase(%superq(CEN_VL&z)) is not valid.;
                %let nerror=%eval(&nerror+1);
            %end;  
        %end;   
        %else %do;
            /**Check if character values are present**/
            %put ERROR: (Model &z: %qupcase(CEN_VL)) Is Required and cannot be missing;
            %let nerror=%eval(&nerror+1);
        %end; 
    %end;   
    /**Event Code for CIF Method**/
    %do z = 1 %to &nmodels;
        %if %sysevalf(%qupcase(%superq(method&z))=CIF,boolean) %then %do;
            /**Check if missing**/
            %if %sysevalf(%superq(ev_vl&z)=,boolean)=0 %then %do;
                %if %sysfunc(notdigit(%sysfunc(compress(%superq(EV_VL&z),.)))) > 0 %then %do;
                    /**Check if character values are present**/
                    %put ERROR: (Model &z: %qupcase(EV_VL)) Must be numeric.  %qupcase(%superq(EV_VL&z)) is not valid.;
                    %let nerror=%eval(&nerror+1);
                %end;  
            %end;   
            %else %do;
                /**Check if character values are present**/
                %put ERROR: (Model &z: %qupcase(EV_VL)) Is Required and cannot be missing when METHOD=CIF;
                %let nerror=%eval(&nerror+1);
            %end; 
        %end;
    %end;             
    /***Error checking for RISKROWWEIGHTS***/
    %do z = 1 %to &nmodels;
        /**Check if RISKROWWEIGHTS is missing when RISKLOCATION is set to BOTTOM**/
        %if %sysevalf(%superq(riskrowweights&z)=,boolean)=0 and %superq(risklocation&z)=BOTTOM %then %do;
            %if %sysfunc(notdigit(%sysfunc(compress(%superq(riskrowweights&z),.-)))) %then %do;
                %put ERROR: (Model &z: RISKROWWEIGHTS) Must be a numeric value (%superq(riskrowweights&z));
                %put ERROR: Macro NEWSURV will cease;
                %let nerror=%eval(&nerror+1);
            %end;
            %else %if %sysevalf(%superq(riskrowweights&z)<0,boolean)=1 or %sysevalf(%superq(riskrowweights&z)>=1.0,boolean)=1 %then %do;
                %put ERROR: (Model &z: RISKROWWEIGHTS) Is not between 0 and 1 (%superq(riskrowweights&z));
                %put ERROR: Macro NEWSURV will cease;
                %let nerror=%eval(&nerror+1);
            %end;
        %end;
        %else %if %superq(risklocation&z)=BOTTOM %then %do;
            %put ERROR: (Model &z: RISKROWWEIGHTS) No risk row weights specified when risklocation=BOTTOM;
            %put ERROR: Macro NEWSURV will cease;
            %let nerror=%eval(&nerror+1);
        %end;
    %end;
        
    /**Error Handling on Individual Model Parameters**/
    %macro _parmcheck(parm, parmlist,kmlist,ciflist,invwtslist,dirlist);
        %do z = 1 %to &nmodels;  
            %if %sysevalf(%superq(&parm.&z)=,boolean)=0 %then %let &parm.&z=%sysfunc(compress(%qupcase(%superq(&parm.&z)),'""'));
            %local _test _z;
            %let _test=;
            %if %sysevalf(%superq(parmlist)^=,boolean) %then %do;
                %do _z=1 %to %sysfunc(countw(&parmlist,|,m));
                    %if %superq(&parm.&z)=%scan(&parmlist,&_z,|,m) %then %let _test=1;
                %end;
                %if &_test ^= 1 %then %do;
                    %put ERROR: (Model &z: %qupcase(&parm)): %superq(&parm.&z) is not a valid value;
                    %put ERROR: (Model &z: %qupcase(&parm)): Possible values are &parmlist;
                    %let nerror=%eval(&nerror+1);
                %end;
            %end;
            %else %do;
                %if %sysevalf(%qupcase(%superq(method&z))=KM,boolean) %then %do;
                    %do _z=1 %to %sysfunc(countw(&kmlist,|,m));
                        %if %superq(&parm.&z)=%scan(&kmlist,&_z,|,m) %then %let _test=1;
                    %end;
                    %if &_test ^= 1 %then %do;
                        %put ERROR: (Model &z: %qupcase(&parm)): %superq(&parm.&z) is not a valid value;
                        %put ERROR: (Model &z: %qupcase(&parm)): Possible values for Kaplan-Meier method are &kmlist;
                        %let nerror=%eval(&nerror+1);
                    %end;
                %end;
                %else %if %sysevalf(%qupcase(%superq(method&z))=CIF,boolean) %then %do;
                    %do _z=1 %to %sysfunc(countw(&ciflist,|,m));
                        %if %superq(&parm.&z)=%scan(&ciflist,&_z,|,m) %then %let _test=1;
                    %end;
                    %if &_test ^= 1 %then %do;
                        %put ERROR: (Model &z: %qupcase(&parm)): %superq(&parm.&z) is not a valid value;
                        %put ERROR: (Model &z: %qupcase(&parm)): Possible values for Competing Risks method are &ciflist;
                        %let nerror=%eval(&nerror+1);
                    %end;
                %end;
                %else %if %sysevalf(%qupcase(%superq(method&z))=INVWTS,boolean) %then %do;
                    %do _z=1 %to %sysfunc(countw(&invwtslist,|,m));
                        %if %superq(&parm.&z)=%scan(&invwtslist,&_z,|,m) %then %let _test=1;
                    %end;
                    %if &_test ^= 1 %then %do;
                        %put ERROR: (Model &z: %qupcase(&parm)): %superq(&parm.&z) is not a valid value;
                        %put ERROR: (Model &z: %qupcase(&parm)): Possible values for Inverse Weights Adjusting method are &invwtslist;
                        %let nerror=%eval(&nerror+1);
                    %end;
                %end;
                %else %if %sysevalf(%qupcase(%superq(method&z))=DIRECT,boolean) %then %do;
                    %do _z=1 %to %sysfunc(countw(&dirlist,|,m));
                        %if %superq(&parm.&z)=%scan(&dirlist,&_z,|,m) %then %let _test=1;
                    %end;
                    %if &_test ^= 1 %then %do;
                        %put ERROR: (Model &z: %qupcase(&parm)): %superq(&parm.&z) is not a valid value;
                        %put ERROR: (Model &z: %qupcase(&parm)): Possible values for Direct Adjusting method are &dirlist;
                        %let nerror=%eval(&nerror+1);
                    %end;
                %end;
            %end;
        %end;
    %mend;
    /**Method**/
    %_parmcheck(method,KM|CIF|INVWTS|DIRECT) 
    /**Y axis Type**/
    %_parmcheck(ytype,PPT|PCT)    
    /**Risk List Location**/
    %_parmcheck(risklocation,|BOTTOM|INSIDE)
    /**Inside Risk List Dividor Line On/Off Option**/
    %_parmcheck(riskdivider,0|1)
    /**Class Descending Order**/
    %_parmcheck(desc,0|1)
    /**Grid lines On/Off**/
    %_parmcheck(gridlines,0|1)
    /**Grid lines Axis**/
    %_parmcheck(gridline_axis,X|Y|BOTH)
    /**Set possible default differences between KM and CIF methods**/
    /**Plot Unadjusted Curves**/
    %_parmcheck(plot_unadjust,,1,1,1|0,1|0);
    %local z;
    %do z = 1 %to &nmodels;
        /**Special PLOT_UNADJUST error catching**/    
        %if %sysevalf(%qupcase(%superq(method&z))=KM) %then %do;
            %if %sysevalf(%superq(plotpval&z)=,boolean) %then %let plotpval&z=LOGRANK;
            %if %sysevalf(%superq(plotpvalmv&z)=,boolean) %then %let plotpvalmv&z=LR;
        %end;
        %else %if %sysevalf(%qupcase(%superq(method&z))=INVWTS) %then %do;
            %if %sysevalf(%superq(plotpval&z)=,boolean) %then %let plotpval&z=LOGRANK;
            %if %sysevalf(%superq(plotpvalmv&z)=,boolean) %then %let plotpvalmv&z=LOGRANK;
            %if %sysevalf(%superq(by&z)^=,boolean) and %sysevalf(%superq(plot_unadjust&z)=1,boolean) %then %do;
                %put ERROR: (Model &z: PLOT_UNADJUST): PLOT_UNADJUST must be 0 when specifying a BY variable and METHOD=INVWTS;
                %let nerror=%eval(&nerror+1);
            %end;
        %end;
        %else %if %sysevalf(%qupcase(%superq(method&z))=DIRECT) %then %do;
            %if %sysevalf(%superq(plotpval&z)=,boolean) %then %let plotpval&z=LOGRANK;
            %if %sysevalf(%superq(plotpvalmv&z)=,boolean) %then %let plotpvalmv&z=LR;
            %if %sysevalf(%superq(by&z)^=,boolean) and %sysevalf(%superq(plot_unadjust&z)=1,boolean) %then %do;
                %put ERROR: (Model &z: PLOT_UNADJUST): PLOT_UNADJUST must be 0 when specifying a BY variable and METHOD=DIRECT;
                %let nerror=%eval(&nerror+1);
            %end;
        %end;
        %else %if %sysevalf(%qupcase(%superq(method&z))=CIF) %then %do;
            %if %sysevalf(%superq(plotpval&z)=,boolean) %then %let plotpval&z=GRAY;
            %if %sysevalf(%superq(plotpvalmv&z)=,boolean) %then %let plotpvalmv&z=WALD;
        %end; 
    %end; 
    /**Confidence Interval Type**/
    %_parmcheck(conftype,,LOG|ASINSQRT|LOGLOG|LINEAR|LOGIT,LOG|ASINSQRT|LOGLOG|LINEAR|LOGIT,LOG|ASINSQRT|LOGLOG|LINEAR|LOGIT,LOG|ASINSQRT|LOGLOG|LINEAR|LOGIT)
    /**CIF Variance Calculation Method**/
    %_parmcheck(cifvar,COUNT|DELTA,COUNT|DELTA)
    /**Plot P-Value**/
    %_parmcheck(plotpval,,SCORE|LR|LOGRANK|LOGRANK_ONESIDED|WILCOXON|WALD|%str( ),WALD|GRAY|%str( ),SCORE|LR|LOGRANK|WILCOXON|WALD|%str( ),SCORE|LR|LOGRANK|WILCOXON|WALD|%str( ))
    /**Adjusted Plot P-Value**/
    %_parmcheck(plotpvalmv,,SCORE|LR|WALD|%str( ),WALD|%str( ),SCORE|LR|LOGRANK|WILCOXON|WALD|%str( ),SCORE|LR|WALD|%str( ))
    /**Interaction Plot P-Value**/
    %_parmcheck(pval_inter,,SCORE|LR|WALD|%str( ),WALD|%str( ),SCORE|LR|LOGRANK|WILCOXON|WALD|%str( ),SCORE|LR|WALD|%str( ))
    /**Interaction Plot P-Value**/
    %_parmcheck(pval_intermv,,SCORE|LR|WALD|%str( ),WALD|%str( ),SCORE|LR|LOGRANK|WILCOXON|WALD|%str( ),SCORE|LR|WALD|%str( ))
    /**S-Reverse Options**/  
    %_parmcheck(sreverse,,0|1,0,0,0)
    /**Censor Values On/Off Option**/
    %_parmcheck(censormarkers,0|1|2)
    /**Class Value Align Option**/
    %_parmcheck(classvalalign,LEFT|CENTER|RIGHT)
    /**Title Align Option**/
    %_parmcheck(titlealign,LEFT|CENTER|RIGHT)
    /**Foot Note Align Option**/
    %_parmcheck(footnotealign,LEFT|CENTER|RIGHT)
    /**X-axis tick value weight Option**/
    %_parmcheck(xtickvalweight,NORMAL|BOLD)
    /**Y-axis tick value weight Option**/
    %_parmcheck(ytickvalweight,NORMAL|BOLD)
    /**Location Options**/
    %_parmcheck(location,INSIDE|OUTSIDE)
    /**Label weight Option**/
    %_parmcheck(lweight,NORMAL|BOLD)
    /**Title weight Option**/
    %_parmcheck(tweight,NORMAL|BOLD)
    /**Footnote weight Option**/
    %_parmcheck(fnweight,NORMAL|BOLD)
    /**Patients-at-Risk weight Option**/
    %_parmcheck(parweight,NORMAL|BOLD)
    /**Hazard Ratio Ties Method Option**/
    %_parmcheck(hrties,BRESLOW|DISCRETE|EFRON|EXACT)
    /**Risk Table Label Location Option**/
    %_parmcheck(risklabellocation,LEFT|ABOVE|)
    /**Risk Table Label Alignment Option**/
    %_parmcheck(risklabelalign,LEFT|CENTER|RIGHT)
    /**Risk Table Patients-at-Risk Subheader Alignment Option**/
    %_parmcheck(paralign,LEFT|CENTER|RIGHT|LABELS)
    /**Risk Table Label weight Option**/
    %_parmcheck(risklabelweight,NORMAL|BOLD)
    /**Risk Numbers and Colors Option**/
    %_parmcheck(riskcolor,0|1)
    /**List Time-points in Plot Summary Table Option**/
    %_parmcheck(listtimepoints,0|1)
    /**Color Statistics in Plot Summary Table Option**/
    %_parmcheck(statcolor,0|1)
    /**Plot Confidence Intervals Enabled**/
    %_parmcheck(plotci,0|1|2)
    /**Plot Confidence Intervals Background Fill Enabled**/
    %_parmcheck(plotcifill,0|1)
    /**Plot X Reference Line Location**/
    %_parmcheck(reflines,|TIMEPOINTS|MEDIANS)
    /**Plot Reference Line Method**/
    %_parmcheck(reflinemethod,FULL|DROP)
    /**Plot Reference Line Axis**/
    %_parmcheck(reflineaxis,X|Y|BOTH)
    /**Symbol Weight**/
    %_parmcheck(symbolweight,NORMAL|BOLD)
    /**RMST Confidence Intervals**/
    %_parmcheck(rmst_ci,0|1)
    /**Model Stat Display Mode**/
    %_parmcheck(model_stats_display,1|2|3)
    /**Draw Underlines**/
    %_parmcheck(draw_underlines,AUTO|0|1)
    
    /**Auto Align Options**/
    %local _z _z2 _test;
    %do z = 1 %to &nmodels;
        /**Check for missing values**/
        %if %sysevalf(%superq(autoalign&z)=,boolean)=0 %then %do _z2=1 %to %sysfunc(countw(%superq(autoalign&z),%str( )));
            /**Check all given values against the possible allowed values**/
            %let _test=;
            %do _z = 1 %to %sysfunc(countw(TOPLEFT|TOP|TOPRIGHT|LEFT|CENTER|RIGHT|BOTTOMLEFT|BOTTOM|BOTTOMRIGHT,|));
                %if %qupcase(%scan(%superq(autoalign&z),&_z2,%str( )))=%scan(TOPLEFT|TOP|TOPRIGHT|LEFT|CENTER|RIGHT|BOTTOMLEFT|BOTTOM|BOTTOMRIGHT,&_z,|,m) %then %let _test=1;
            %end;
            /**If any values are not in the possible list then throw an error**/
            %if &_test ^= 1 %then %do;
                %put ERROR: (Model &z: %qupcase(autoalign)): %qupcase(%scan(%superq(autoalign&z),&_z2,%str( ))) is not in the list of valid values;
                %put ERROR: (Model &z: %qupcase(autoalign)): Possible values are TOPLEFT|TOP|TOPRIGHT|LEFT|CENTER|RIGHT|BOTTOMLEFT|BOTTOM|BOTTOMRIGHT;
                %let nerror=%eval(&nerror+1);
            %end;
        %end;
        %else %do;
            /**If missing then show error**/
            %put ERROR: (Model &z: %qupcase(autoalign)): Cannot be missing;
            %put ERROR: (Model &z: %qupcase(autoalign)): Possible values are TOPLEFT|TOP|TOPRIGHT|LEFT|CENTER|RIGHT|BOTTOMLEFT|BOTTOM|BOTTOMRIGHT;
            %let nerror=%eval(&nerror+1);
        %end;           
    %end;
    /**Plot Display Variables**/
    %local _z _z2 _test _displaylist;
    %do z = 1 %to &nmodels;
        %if %sysevalf(%qupcase(%superq(display&z))=STANDARD,boolean) %then %do;
            %if %sysevalf(%superq(class&z)^=,boolean) %then %let display&z=legend ev_n median hr timelist pval tablecomments;
            %else %if %sysevalf(%superq(class&z)=,boolean) %then %let display&z=ev_n median timelist tablecomments;
        %end;
        %if %sysevalf(%qupcase(%superq(method&z))=KM,boolean) %then 
            %let _displaylist=legend|hr|median|total|event|ev_n|n_ev|timelist|pval|cindex|tablecomments|totalmv|eventmv|ev_nmv|n_evmv|hrmv|pvalmv|cindexmv|covpval|covpvalmv|pval_inter|pval_intermv;
			/*|RMST|RMTL|RMST_PVAL;*/
        %else %if %sysevalf(%qupcase(%superq(method&z))^=CIF,boolean) %then 
            %let _displaylist=legend|hr|median|total|event|ev_n|n_ev|timelist|covpval|pval|pval_inter|cindex|tablecomments;
        %else %if %sysevalf(%qupcase(%superq(method&z))=CIF,boolean) %then %let _displaylist=legend|hr|median|total|event|ev_n|n_ev|timelist|pval|tablecomments|totalmv|eventmv|ev_nmv|n_evmv|hrmv|pvalmv|covpval|covpvalmv|pval_inter|pval_intermv;
        /**Check for missing values**/
        %if %sysevalf(%superq(display&z)=,boolean)=0 %then %do _z2=1 %to %sysfunc(countw(%superq(display&z),%str( )));
            /**Check all given values against the possible allowed values**/
            %let _test=;
            %do _z = 1 %to %sysfunc(countw(&_displaylist,|));
                %if %qupcase(%scan(%superq(display&z),&_z2,%str( )))=%scan(%qupcase(&_displaylist),&_z,|,m) %then %let _test=1;
            %end;
            /**If any values are not in the possible list then throw an error**/
            %if &_test ^= 1 %then %do;
                %put ERROR: (Model &z: %qupcase(display)): %qupcase(%scan(%superq(display&z),&_z2,%str( ))) is not in the list of valid values;
                %put ERROR: (Model &z: %qupcase(display)): Possible values are %qupcase(&_displaylist);
                %let nerror=%eval(&nerror+1);
            %end;
        %end;
    %end;
    /**Plot Patients-at-Risk Display Variables**/
    %local _z _z2 _test _displaylist;
    %do z = 1 %to &nmodels;
        %let _displaylist=par|ncens|nevents|par_ncens|par_nevents;
        %if %sysevalf(%superq(risklist&z)^=,boolean) %then %do;
            /**Check for missing values**/
            %if %sysevalf(%superq(pardisplay&z)=,boolean)=0 %then %do _z2=1 %to %sysfunc(countw(%superq(pardisplay&z),%str( )));
                /**Check all given values against the possible allowed values**/
                %let _test=;
                %do _z = 1 %to %sysfunc(countw(&_displaylist,|));
                    %if %qupcase(%scan(%superq(pardisplay&z),&_z2,%str( )))=%scan(%qupcase(&_displaylist),&_z,|,m) %then %let _test=1;
                %end;
                /**If any values are not in the possible list then throw an error**/
                %if &_test ^= 1 %then %do;
                    %put ERROR: (Model &z: %qupcase(pardisplay)): %qupcase(%scan(%superq(pardisplay&z),&_z2,%str( ))) is not in the list of valid values;
                    %put ERROR: (Model &z: %qupcase(pardisplay)): Possible values are %qupcase(&_displaylist);
                    %let nerror=%eval(&nerror+1);
                %end;
            %end;
        %end;
    %end;
    /**Line Pattern Variables**/
    %macro _linepattern(parm=,_patternlist=AUTO|SOLID|SHORTDASH|MEDIUMDASH|LONGDASH|MEDIUMDASHSHORTDASH|
        DASHDASHDOT|DASH|LONGDASHSHORTDASH|DOT|THINDOT|SHORTDASHDOT|MEDIUMDASHDOTDOT);
        %local _z _z2 _test;
        %do z = 1 %to &nmodels;
            /**Check for missing values**/
            %if %sysevalf(%superq(&parm.&z)=,boolean)=0 %then %do _z2=1 %to %sysfunc(countw(%superq(&parm.&z),%str( )));
                %let _test=;
                /**Check if values are either in the approved list, or are between 1 and 46**/
                %if %sysfunc(notdigit(%scan(%superq(&parm.&z),&_z2,%str( ))))>0 %then %do _z = 1 %to %sysfunc(countw(&_patternlist,|));
                    %if %qupcase(%scan(%superq(&parm.&z),&_z2,%str( )))=%scan(%qupcase(%sysfunc(compress(&_patternlist))),&_z,|,m) %then %let _test=1;
                %end;
                %else %if %scan(%superq(&parm.&z),&_z2,%str( )) ge 1 and %scan(%superq(&parm.&z),&_z2,%str( )) le 46 %then %let _test=1;
                %if &_test ^= 1 %then %do;
                    /**Throw error**/
                    %put ERROR: (Model &z: %qupcase(&parm.)): %qupcase(%scan(%superq(&parm.&z),&_z2,%str( ))) is not in the list of valid values;
                    %put ERROR: (Model &z: %qupcase(&parm.)): Possible values are %qupcase(&_patternlist) or Numbers Between 1 and 46;
                    %let nerror=%eval(&nerror+1);
                %end;
            %end;
            %else %do;
                /**Throw error**/
                %put ERROR: (Model &z: %qupcase(&parm.)): %qupcase(%superq(&parm.&z)) is not in the list of valid values;         
                %put ERROR: (Model &z: %qupcase(&parm.)): Possible values are %qupcase(&_patternlist) or Numbers Between 1 and 46;
                %let nerror=%eval(&nerror+1);       
            %end;
        %end;
    %mend;
    /**Plot Line Patterns**/
    %_linepattern(parm=pattern)
    %_linepattern(parm=pattern_adjust)
    /**Patients-at-Risk INSIDE option Dividor Line Style**/
    %_linepattern(parm=riskdivstyle)
    /**Plot Confidence Bounds Confidence Interval**/
    %_linepattern(parm=plotcilinepattern)
    /**Plot X Referenceline line style**/
    %_linepattern(parm=reflinepattern)
    /**Grid Line Pattern**/
    %_linepattern(parm=gridline_pattern)
    /**Underline Pattern**/
    %_linepattern(parm=ul_pattern)
        
    /**Range Value Check**/
    %macro _rangecheck(parm=,min=,max=,incmax=,incmin=);
        %do z = 1 %to &nmodels;
            /**Check for missing values**/
            %if %sysevalf(%superq(&parm.&z)=,boolean)=0 %then %do;
                %if %sysfunc(notdigit(%sysfunc(compress(%superq(&parm.&z),-.)))) > 0 %then %do;
                    /**Checks for character values**/
                    %put ERROR: (Model &z: %qupcase(&parm.)) Must be numeric. %qupcase(%superq(&parm.&z)) is not valid.;
                    %let nerror=%eval(&nerror+1);
                %end;
                %else %if %superq(&parm.&z) le &min and &incmin=0 %then %do;
                    /**Checks if less than or equal to min**/
                    %put ERROR: (Model &z: %qupcase(&parm.)) Cannot be less than or equal to &min (%superq(&parm.&z));
                    %let nerror=%eval(&nerror+1);
                %end;
                %else %if %superq(&parm.&z) lt &min %then %do;
                    /**Checks if less than min**/
                    %put ERROR: (Model &z: %qupcase(&parm.)) Cannot be less than &min (%superq(&parm.&z));
                    %let nerror=%eval(&nerror+1);
                %end;
                %else %if %superq(&parm.&z) ge &max and &incmax=0 %then %do;
                    /**Checks if greater than or equal to max**/
                    %put ERROR: (Model &z: %qupcase(&parm.)) Cannot be Greater Than or Equal to &max;
                    %let nerror=%eval(&nerror+1);
                %end;
                %else %if %superq(&parm.&z) gt &max %then %do;
                    /**Checks if greater than max**/
                    %put ERROR: (Model &z: %qupcase(&parm.)) Cannot be Greater Than &max;
                    %let nerror=%eval(&nerror+1);
                %end;
            %end;
        %end; 
    %mend; 
    /**X Axis Minimum offset Value**/
    %_rangecheck(parm=xminoffset,min=0,max=1,incmax=0,incmin=1)
    /**Y Axis Minimum offset Value**/
    %_rangecheck(parm=yminoffset,min=0,max=1,incmax=0,incmin=1)
    /**X Axis Maximum offset Value**/
    %_rangecheck(parm=xmaxoffset,min=0,max=1,incmax=0,incmin=1)
    /**Y Axis Maximum offset Value**/
    %_rangecheck(parm=ymaxoffset,min=0,max=1,incmax=0,incmin=1)
    /**Plot Confidence Intervals transparency**/
    %_rangecheck(parm=plotcifilltransparency,min=0,max=1,incmax=1,incmin=1)
    
    /**Y Axis Maximum Value**/
    %do z = 1 %to &nmodels;
        /**Check for missing values**/
        %if %sysevalf(%superq(ymax&z)=,boolean)=0 %then %do;
            %if %sysfunc(notdigit(%sysfunc(compress(%superq(ymax&z),-.)))) > 0 %then %do;
                /**Checks for character values**/
                %put ERROR: (Model &z: %qupcase(ymax)) Must be numeric.  %qupcase(%superq(ymax&z)) is not valid.;
                %let nerror=%eval(&nerror+1);
            %end;
            %else %if %superq(ymax&z) le %superq(ymin&z) %then %do;
                /**Makes sure the maximum is not less than the minimum**/
                %put ERROR: (Model &z: %qupcase(ymax)) Cannot be less than or equal to YMIN (%superq(ymax&z) vs. %superq(ymin&z));
                %let nerror=%eval(&nerror+1);
            %end;
            %else %if (%sysevalf(&&ymax&z>100,boolean) and %qupcase(%superq(ytype&z))=PCT) or 
                (%sysevalf(&&ymax&z>1,boolean) and %qupcase(%superq(ytype&z))=PPT) %then %do;
                /**Makes sure the maximum cannot be greater than the maximum survival estimate**/
                %if %qupcase(%superq(ytype))=PPT %then %do;
                    %put ERROR: (Model &z: %qupcase(ymax)) Cannot be Greater Than 1 When YTYPE=PPT;
                    %let nerror=%eval(&nerror+1);
                %end;
                %else %if %qupcase(%superq(ytype))=PCT %then %do;
                    %put ERROR: (Model &z: %qupcase(ymax)) Cannot be Greater Than 100 When YTYPE=PCT;
                    %let nerror=%eval(&nerror+1);
                %end;
            %end;
        %end;
        %else %if %qupcase(%superq(ytype&z))=PPT %then %let ymax&z=1;
        %else %if %qupcase(%superq(ytype&z))=PCT %then %let ymax&z=100;
    %end;
    
    /**Y Axis Increment Value**/
    %do z = 1 %to &nmodels;
        /**Check for missing values**/
        %if %sysevalf(%superq(yincrement&z)=,boolean)=0 %then %do;
            %if %sysfunc(notdigit(%sysfunc(compress(%superq(yincrement&z),-.)))) > 0 %then %do;
                /**Checks for character values**/
                %put ERROR: (Model &z: %qupcase(yincrement)) Must be numeric.  %qupcase(%superq(yincrement&z)) is not valid.;
                %let nerror=%eval(&nerror+1);
            %end;
            %else %if %sysevalf(&&yincrement&z>(&&ymax&z-&&ymin&z),boolean) %then %do;                    
                /**Makes sure the increment is not greater than the distance between max and min**/
                %put ERROR: (Model &z: %qupcase(yincrement)) Cannot be less than or equal to difference between YMAX and YMIN (%superq(yincrement&z) vs. %sysfunc(sum(%superq(ymax&z),-%superq(ymin&z))));
                %let nerror=%eval(&nerror+1);
            %end;
            %else %if %superq(yincrement&z) le 0 %then %do;
                /**Makes sure the increment is greater than zero**/
                %put ERROR: (Model &z: %qupcase(yincrement)) Cannot be less than or equal to 0;
                %let nerror=%eval(&nerror+1);
            %end;
        %end;
        %else %if %qupcase(%superq(ytype&z))=PPT %then %let yincrement&z=0.1;
        %else %if %qupcase(%superq(ytype&z))=PCT %then %let yincrement&z=10;
    %end;
            
    /**Error Handling on Global Parameters**/
    %macro _gparmcheck(parm, parmlist);          
        %local _test _z;
        /**Check if values are in approved list**/
        %do _z=1 %to %sysfunc(countw(&parmlist,|,m));
            %if %qupcase(%superq(&parm))=%qupcase(%scan(&parmlist,&_z,|,m)) %then %let _test=1;
        %end;
        %if &_test ^= 1 %then %do;
            /**If not then throw error**/
            %put ERROR: (Global: %qupcase(&parm)): %superq(&parm) is not a valid value;
            %put ERROR: (Global: %qupcase(&parm)): Possible values are &parmlist;
            %let nerror=%eval(&nerror+1);
        %end;
    %mend;
    /**Plot On/Off Options**/
    %_gparmcheck(plot,0|1)
    /**Summary On/Off Options**/
    %_gparmcheck(summary,0|1)
    /**Plot Wall On/Off Options**/
    %_gparmcheck(showwalls,0|1)
    /**New Table On/Off Options**/
    %_gparmcheck(newtable,0|1)
    /**Lattice Order Options**/
    %_gparmcheck(order,COLUMNMAJOR|ROWMAJOR)
    /**Destination Options**/
    %_gparmcheck(destination,RTF|PDF|HTML|EXCEL|POWERPOINT)
    /**Overall Title weight Option**/
    %_gparmcheck(ovtweight,NORMAL|BOLD)
    /**Overall Footnote weight Option**/
    %_gparmcheck(ovfnweight,NORMAL|BOLD)
    /**Table Header weight Option**/
    %_gparmcheck(tableheaderweight,MEDIUM|BOLD)
    /**Table Footnote weight Option**/
    %_gparmcheck(tablefootnoteweight,MEDIUM|BOLD)
    /**Table Data Columns weight Option**/
    %_gparmcheck(tabledataweight,MEDIUM|BOLD)
    /**Overall Title Align Option**/
    %_gparmcheck(ovtitlealign,LEFT|CENTER|RIGHT)
    /**Overall Foot Note Align Option**/
    %_gparmcheck(ovfootnotealign,LEFT|CENTER|RIGHT) 
    /**Border around plot image Option**/
    %_gparmcheck(border,0|1)
    /**Transparent Background Option**/
    %_gparmcheck(transparent,0|1)  
    /**Table Background Shading**/
    %_gparmcheck(tableshading,0|1)
    /**Uniform Height below X-axis**/
    %_gparmcheck(uniformheight,0|1)
    /**Merge covariate p-values into overall p-value column**/
    %_gparmcheck(tablemergepval,0|1)
    /**Turn debugging options on/off**/
    %_gparmcheck(debug,0|1)
    /*Tiff Device Check*/
    %if %sysevalf(%qupcase(&plottype)=TIFF,boolean) or  %sysevalf(%qupcase(&plottype)=TIF,boolean) %then %do;
        ods output gdevice=_gdevice;
        proc gdevice catalog=sashelp.devices nofs;
            list _all_;
        run;
        quit;
        %global _tifflist _tiffcheck;
        proc sql noprint;
            select 1 into :_tiffcheck from _gdevice where upcase(name)=upcase("&tiffdevice");
            select distinct upcase(name) into :_tifflist separated by '|' from _gdevice
                where substr(upcase(name),1,3)='TIF';
            %if %sysevalf(%superq(_tiffcheck)^=1,boolean) %then %do;
                /**If not then throw error**/
                %put ERROR: (Global: TIFFDEVICE): %qupcase(%superq(tiffdevice)) is not on the installed list of devices;
                %put ERROR: (Global: TIFFDEVICE): Please select from the following list &_tifflist;
                %let nerror=%eval(&nerror+1);
            %end;
            drop table _gdevice;
        quit;
    %end;
    /**Scalable Vector Graphics On/Off Options**/
    %_gparmcheck(svg,0|1);
    
    /**Error Handling on Global Parameters Involving units**/
    %macro _gunitcheck(parm);
        %if %sysevalf(%superq(&parm)=,boolean)=1 %then %do;
            /**Check if missing**/
            %put ERROR: (Global: %qupcase(&parm)) Cannot be set to missing;
            %let nerror=%eval(&nerror+1);
        %end;
        %else %if %sysfunc(compress(%superq(&parm),ABCDEFGHIJKLMNOPQRSTUVWXYZ,i)) lt 0 %then %do;
            /**Throw error**/
            %put ERROR: (Global: %qupcase(&parm)) Cannot be less than zero (%qupcase(%superq(&parm)));
            %let nerror=%eval(&nerror+1);
        %end;
    %mend;
    /**Overall Title Font Size**/
    %_gunitcheck(ovtsize)
    /**Overall Footnote Font Size**/
    %_gunitcheck(ovfnsize)
    /**Plot Width**/
    %_gunitcheck(width)
    %if %sysevalf(%qupcase(%superq(plottype))=TIFF,boolean) or
        %sysevalf(%qupcase(%superq(plottype))=TIF,boolean) %then %do;
        %if %sysfunc(find(%superq(width),px,i))>0 %then %do;
            /**Throw error**/
            %put ERROR: (Global: WIDTH) Must use units of IN when PLOTTYPE=%qupcase(&plottype);
            %let nerror=%eval(&nerror+1);
        %end;
    %end;
    /**Plot Height**/
    %_gunitcheck(height)
    %if %sysevalf(%qupcase(%superq(plottype))=TIFF,boolean) or
        %sysevalf(%qupcase(%superq(plottype))=TIF,boolean) %then %do;
        %if %sysfunc(find(%superq(height),px,i))>0 %then %do;
            /**Throw error**/
            %put ERROR: (Global: HEIGHT) Must use units of IN when PLOTTYPE=%qupcase(&plottype);
            %let nerror=%eval(&nerror+1);
        %end;
    %end;
    /**Table Header Font Size**/
    %_gunitcheck(tableheadersize)
    /**Table Footnote Font Size**/
    %_gunitcheck(tablefootnotesize)
    /**Table Data Columns Font Size**/
    %_gunitcheck(tabledatasize)
    /**Table Total Count Column Width**/
    %_gunitcheck(ttotalwidth)
    /**Table Events Count Column Width**/
    %_gunitcheck(teventwidth)
    /**Table Combined Total Counts and Events Count Column Width**/
    %_gunitcheck(tev_nwidth)
    /**Table Median Column Width**/
    %_gunitcheck(tmedianwidth)
    /**Table Hazard Ratio Column Width**/
    %_gunitcheck(thrwidth)
    /**Table Time point estimates Column Width**/
    %_gunitcheck(ttimelistwidth)
    /**Table P-value Column Width**/
    %_gunitcheck(tpvalwidth)
    
    /**Error Handling on Individual Model Numeric Variables**/
    %macro _gnumcheck(parm, min);
        /**Check if missing**/
        %if %sysevalf(%superq(&parm)^=,boolean) %then %do;
            %if %sysfunc(notdigit(%sysfunc(compress(%superq(&parm),-.)))) > 0 %then %do;
                /**Check if character value**/
                %put ERROR: (Global: %qupcase(&parm)) Must be numeric.  %qupcase(%superq(&parm)) is not valid.;
                %let nerror=%eval(&nerror+1);
            %end;
            %else %if %superq(&parm) < &min %then %do;
                /**Makes sure number is not less than the minimum**/
                %put ERROR: (Global: %qupcase(&parm)) Must be greater than or equal to %superq(min). %qupcase(%superq(&parm)) is not valid.;
                %let nerror=%eval(&nerror+1);
            %end;
        %end;
        %else %do;
            /**Throw Error**/
            %put ERROR: (Global: %qupcase(&parm)) Cannot be missing;
            %put ERROR: (Global: %qupcase(&parm)) Possible values are numeric values greater than or equal to %superq(min);
            %let nerror=%eval(&nerror+1);           
        %end;       
    %mend;
    /**Digital Pixels per Inch Value**/
    %_gnumcheck(dpi,100)
    /**Anti-alias Maximum Value**/
    %_gnumcheck(antialiasmax,1)
    /**Number of Lattice Columns**/
    %_gnumcheck(columns,1)
    /**Number of Lattice Rows**/
    %_gnumcheck(rows,1)
    /**Number of Models**/
    %_gnumcheck(nmodels,1)
    %if &plot=1 and &nerror=0 %then %do;
        %if %sysevalf(&rows*&columns < &nmodels,boolean) %then %do;
            /**Throw Error**/
            %put ERROR: (Global: NMODELS) if PLOT=1 then NMODELS (&nmodels) must be less than or equal to ROWS*COLUMNS (%sysevalf(&rows*&columns));
            %let nerror=%eval(&nerror+1);
        %end;
    %end;
    
    /**Summary Table Display Variables**/
    %local _z2 _adjtest _unadjtest _list;
    /**Check if missing**/
    %let _adjtest=0;%let _unadjtest=0;
    %if %qupcase(&tabledisplay)=STANDARD %then %do;
        %do _z2 = 1 %to &nmodels;
            %if %qupcase(&&method&_z2)=INVWTS or %qupcase(&&method&_z2)=DIRECT %then %let _adjtest=1;
            %else %if (%sysevalf(%superq(classcov&_z2)^=,boolean) or %sysevalf(%superq(contcov&_z2)^=,boolean)) %then %let _unadjtest=3;
            %else %if %sysevalf(%superq(class&_z2)^=,boolean) and &_unadjtest<3 %then %let _unadjtest=2;
            %else %if &_unadjtest<2 %then %let _unadjtest=1;
        %end;
        %if &_adjtest=1 and &_unadjtest=0 %then %let tabledisplay=title footnote ev_nmv medianmv timelistmv hrmv pvalmv;
        %else %if &_adjtest=1 and &_unadjtest>0 %then %let tabledisplay=title footnote ev_n median hr timelist pval ev_nmv medianmv timelistmv hrmv pvalmv;
        %else %if &_unadjtest=3 %then %let tabledisplay=title footnote ev_n median hr timelist pval ev_nmv hrmv pvalmv;
        %else %if &_unadjtest=2 %then %let tabledisplay=title footnote ev_n median hr timelist covpval pval;
        %else %let tabledisplay=title footnote ev_n median timelist;
    %end;
    %if %sysevalf(%superq(tabledisplay)=,boolean)=0 %then %do _z2=1 %to %sysfunc(countw(%superq(tabledisplay),%str( )));
        %let _test=;
        %let _list=title|footnote|total|event|ev_n|n_ev|median|hr|timelist|pval|covpval|cindex|totalmv|eventmv|ev_nmv|n_evmv|hrmv|medianmv|timelistmv|pvalmv|covpvalmv|cindexmv|pval_inter|pval_intermv;
		/*|RMST|RMTL|RMST_PVAL;*/
        /**Check if submitted values are in approved list**/
        %do _z = 1 %to %sysfunc(countw(%qupcase(&_list),|));
            %if %qupcase(%scan(%superq(tabledisplay),&_z2,%str( )))=
                %scan(%qupcase(&_list),&_z,|,m) %then %let _test=1;
        %end;
        %if &_test ^= 1 %then %do;
            /**Throw errors**/
            %put ERROR: (Global: %qupcase(tabledisplay)): %qupcase(%scan(%superq(tabledisplay),&_z2,%str( ))) is not in the list of valid values;
            %put ERROR: (Global: %qupcase(tabledisplay)): Possible values are %qupcase(&_list);
            %let nerror=%eval(&nerror+1);
        %end;
    %end;
    /**If confidence interval colors are missing, set them to line colors**/
    %local z;
    %do z = 1 %to &nmodels;
        %if %sysevalf(%superq(plotcifillcolor&z)=,boolean) %then %let plotcifillcolor&z=%superq(color&z);
        %if %sysevalf(%superq(plotcilinecolor&z)=,boolean) %then %let plotcilinecolor&z=%superq(color&z);
    %end;
    /*** If any errors exist, stop macro and send to end ***/
    %if &nerror > 0 %then %do;
        %put ERROR: &nerror pre-run errors listed;
        %put ERROR: Macro NEWSURV will cease;
        %goto errhandl;
    %end;    
    
    %if &debug=1 %then %do;
        options notes mprint;
    %end;
    
    %do z = 1 %to &nmodels;/**Start of Analysis Section**/
        %local nerror_run;
        %let nerror_run=0;
        data _null_;
            set &data (obs=1);
            /**If missing x label, grab label or name of time variable**/
            %if %sysevalf(%superq(xlabel&z)=,boolean) %then %do;
                call symput("xlabel&z",vlabel(%superq(time&z)));
            %end;
            %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do;
                %local format&z label&z;            
                /**Selects format for class variable**/
                call symput("label&z",vlabel(%superq(class&z)));
                call symput("format&z",vformat(%superq(class&z)));
                %if %sysevalf(%superq(by&z)^=,boolean) and %sysevalf(%superq(bylabel&z)=,boolean) %then %do;
                    call symput("bylabel&z",vlabel(%superq(by&z)));
                %end; 
            %end;
        run;
        /**If missing y label, make generic label depending on y axis type**/
        %if %sysevalf(%superq(ylabel&z)=,boolean)=1 %then %do;
            %if %qupcase(%superq(ytype&z))=PCT %then %do;
                %if %qupcase(%superq(method&z))=CIF or %superq(sreverse&z)=1 %then %let ylabel&z=Percent With Event;
                %else %let ylabel&z=Percent Without Event;
            %end;
            %else %do;
                %if %qupcase(%superq(method&z))=CIF or %superq(sreverse&z)=1 %then %let ylabel&z=Proportion With Event;
                %else %let ylabel&z=Proportion Without Event; 
            %end;               
        %end;
        /**Set Underline option if currently set to AUTO**/
        %if %sysevalf(%superq(draw_underlines&z)=AUTO,boolean) %then %do;
            %if %sysevalf(%superq(by&z)^=,boolean) %then %let draw_underlines&z=1;
            %else %let draw_underlines&z=0; 
        %end;
        /**Apply where clause**/
        proc sql;
            create table _tempdsn&z as
                select * from &data (
                    %if %sysevalf(%superq(where&z)=,boolean)=0 %then %do;
                    where=(&&where&z)
                    %end;);
            %local datcheck;
            
            %if %sysfunc(exist(_tempdsn&z))>0 %then %do;
                select count(*) into :datcheck from _tempdsn&z;
            %end;
            %else %let datcheck=0;
        quit;
        %if &datcheck =0 %then %do;
            %put ERROR: (Model &z: WHERE): Issue parsing the WHERE clause;
            %let nerror_run=%eval(&nerror_run+1);
            %goto errhandl2; 
        %end;
        /**Create temporary dataset for analysis**/
        /**Changes class variable to character variable by applying format**/
        /**Makes sure there are no missing values in the key variables**/
        data _tempdsn&z;
            merge %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do; 
                     _tempdsn&z (keep=%superq(class&z) rename=(%superq(class&z)=_tempvar_))
                     %if %sysevalf(%superq(by&z)=,boolean)=0 %then %do; 
                         _tempdsn&z (keep=%superq(by&z) rename=(%superq(by&z)=_tempvarb_))
                      %end;
                  %end;
                  _tempdsn&z (keep=%superq(time&z) rename=(%superq(time&z)=_time_))
                  _tempdsn&z (keep=%superq(cens&z) rename=(%superq(cens&z)=_cens_))
                  %do i=1 %to %sysfunc(countw(%superq(classcov&z),%str( ))); 
                      _tempdsn&z (keep=%scan(%superq(classcov&z),&i,%str( )) rename=(%scan(%superq(classcov&z),&i,%str( ))=_classcov_&i))
                  %end;
                  %do i=1 %to %sysfunc(countw(%superq(contcov&z),%str( ))); 
                      _tempdsn&z (keep=%scan(%superq(contcov&z),&i,%str( )) rename=(%scan(%superq(contcov&z),&i,%str( ))=_contcov_&i))
                  %end;                                             
                  %do i=1 %to %sysfunc(countw(%superq(strata&z),%str( ))); 
                      _tempdsn&z (keep=%scan(%superq(strata&z),&i,%str( )) rename=(%scan(%superq(strata&z),&i,%str( ))=_strata_&i))
                  %end;
                  %if %sysfunc(notdigit(0%sysfunc(compress(%superq(landmark&z),.-)))) > 0 %then %do;
                      _tempdsn&z (keep=%superq(landmark&z) rename=(%superq(landmark&z)=_landmark_))
                  %end;;
            %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do; 
                length _class_ $300. _class_order _class_ref 8.;
                _class_=strip(vvalue(_tempvar_));
                _class_order=.;
                _class_ref=.;
                drop _tempvar_;
                if ^missing(_class_);
                %if %sysevalf(%superq(by&z)=,boolean)=0 %then %do; 
                    length _by_ $300. _by_order 8.;
                    _by_=strip(vvalue(_tempvarb_));
                    _by_order=.;
                    drop _tempvarb_;
                    if ^missing(_by_);
                %end;
            %end;
            %if %sysevalf(%superq(by&z)=,boolean) %then %do;
                length _by_ $300. _by_order 8.;
                _by_='';
                _by_order=1;
            %end;
            if ^missing(_time_) and ^missing(_cens_) 
                %if %sysfunc(notdigit(0%sysfunc(compress(%superq(landmark&z),.-)))) > 0 %then %do; and ^missing(_landmark_) %end;
                %if %sysevalf(%superq(strata&z)=,boolean)=0 %then %do i = 1 %to %sysfunc(countw(%superq(strata&z),%str( ))); and ^missing(_strata_&i) %end;;        
            %if %sysevalf(%superq(xdivisor&z)=,boolean)=0 %then %do;
                _time_=_time_/%superq(xdivisor&z);
            %end;
            /**Apply landmark if not missing**/
            %if %sysevalf(%superq(landmark&z)=,boolean)=0 %then %do;
                _time_=_time_-
                    %if %sysfunc(notdigit(0%sysfunc(compress(%superq(landmark&z),.-)))) > 0 %then %do;
                        _landmark_
                    %end;
                    %else %do;
                        %superq(landmark&z)
                    %end;;
                if _time_ gt 0;
            %end;
            
            %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;
                if missing(_class_)=0
                    %if %sysevalf(%superq(classcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(classcov&z),%str( )));
                        and missing(_classcov_&i)=0
                    %end;
                    %if %sysevalf(%superq(contcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(contcov&z),%str( )));
                        and missing(_contcov_&i)=0
                    %end; then do;
                    _multi_=1;output;
                end;
                %if %sysevalf(%superq(classcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(classcov&z),%str( ))); 
                    if vtype(_classcov_&i)='N' then _classcov_&i=1;
                    else _classcov_&i='1';
                %end;
                %if %sysevalf(%superq(contcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(contcov&z),%str( ))); _contcov_&i=1; %end;
                _multi_=0;output;
            %end;
            %else %do;
                _multi_=0;output;
            %end;
        run;
        proc sort data=_tempdsn&z;
            by _multi_;
        run;
        proc sql noprint;
            %local nby_&z b;
            %if %sysevalf(%superq(by&z)=,boolean)=0 %then %do;
                /**Saves number of by levels into macro variable**/
                select count(distinct _by_) into :nby_&z from _tempdsn&z
                    where _multi_=0;

                /***Save by values into macro variables***/
                %do i = 1 %to %superq(nby_&z);
                    %local by_&z._&i.;
                %end;   
                /**If order set to auto then select by values in order sorted by by variable**/
                %if %sysevalf(%superq(byorder&z)=,boolean) %then %do;
                    select distinct _by_ into :by_&z._1-
                        from _tempdsn&z
                        where _multi_=0 
                        order by _by_
                        %if %superq(desc&z)=1 %then %do; DESC %end; /**If desc=DESC then reverse order of by variable**/;
                %end;
                %else %do;
                    /**Select by values per requested order**/
                    %local ___bylevels_&z;
                    select distinct _by_ into :___bylevels_&z separated by '|'
                        from _tempdsn&z
                        where _multi_=0
                        order by _by_;
                    %local i2;
                    %do i = 1 %to %superq(nby_&z);
                        %if %superq(desc&z)=1 %then %let i2=%sysfunc(abs(%sysfunc(sum(&i,-%superq(nby_&z),-1))));
                        %else %let i2=&i;
                        %let by_&z._&i = %scan(%superq(___bylevels_&z),%scan(%superq(byorder&z),&i2),|,m);
                    %end;
                %end;
                update _tempdsn&z
                    set _by_order=case (_by_)
                        %do i=1 %to &&nby_&z;
                            when "%superq(by_&z._&i)" then &i
                        %end;
                        else . end
                    ;
            %end;
            %else %let nby_&z=1;
            %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do;
                %local nclass_&z;    
                /**Saves number of total class/by combination levels into macro variable**/
                select count(*) into :nclass_&z 
                    from (select distinct _by_,_class_ from _tempdsn&z
                            where _multi_=0); 
                %do b=1 %to &&nby_&z;
                    %local nclass_&b._&z;    
                    /**Saves number of total classcombination levels within by level into macro variable**/
                    select count(distinct _class_) into :nclass_&b._&z 
                        from _tempdsn&z where _multi_=0 and _by_order=&b; 
                    %local ___classlevels_&b._&z;
                    select distinct _class_ into :___classlevels_&b._&z separated by '|'
                        from _tempdsn&z
                        where _multi_=0 and _by_order=&b
                        order by _class_;
                    /***Save class values into macro variables***/
                    %do i = 1 %to %superq(nclass_&b._&z);
                        %local class_&b._&z._&i nclass_&b._&z._&i;
                    %end;   
                
                    /**If order set to auto then select class values in order sorted by class variable**/
                    %if %sysevalf(%qscan(%superq(classorder&z),&b,`,m)=,boolean) %then %do;
                        select distinct _class_ into :class_&b._&z._1-
                            from _tempdsn&z
                            where _multi_=0 and _by_order=&b
                            order by _class_
                            %if %superq(desc&z)=1 %then %do; DESC %end; /**If desc=DESC then reverse order of class variable**/;
                    %end;
                    %else %do;
                        /**Select class values per requested order**/
                        %local i2;
                        %do i = 1 %to %superq(nclass_&b._&z);
                            %if %superq(desc&z)=1 %then %let i2=%sysfunc(abs(%sysfunc(sum(&i,-%superq(nclass_&b._&z),-1))));
                            %else %let i2=&i;
                            %let class_&b._&z._&i = %scan(%superq(___classlevels_&b._&z),%scan(%qscan(%superq(classorder&z),&b,`,m),&i2),|,m);
                        %end;
                        
                        /***Error Check if Enough Class Order Levels were Specified if Not Missing***/
                        %if %superq(nclass_&b._&z) ^= %sysfunc(countw(%qscan(%superq(classorder&z),&b,`,m),%str( ))) %then %do;                            
                            %if %sysevalf(%superq(by&z)=,boolean) %then %do;
                                %put ERROR: (Model &z: CLASSORDER) Different number of classorders specified than number of class levels (%qscan(%superq(classorder&z),&b,`,m) vs. %sysfunc(strip(%superq(nclass_&b._&z))) Class Levels);
                            %end;
                            %else %do;
                                %put ERROR: (Model &z BY group &b: CLASSORDER) Different number of classorders specified than number of class levels (%qscan(%superq(classorder&z),&b,`,m) vs. %sysfunc(strip(%superq(nclass_&b._&z))) Class Levels);
                            %end;                            
                            %let nerror_run=%eval(&nerror_run+1);
                        %end; 
                        /***Error Check if each level of class was specified in the class order list***/    
                        %if &nerror_run=0 %then %do;
                            %do i = 1 %to %superq(nclass_&b._&z);
                                %local _test;
                                %let _test=;
                                %do j = 1 %to %sysfunc(countw(%qscan(%superq(classorder&z),&b,`,m),%str( )));
                                    %if &i = %scan(%qscan(%superq(classorder&z),&b,`,m),&j,%str( )) %then %let _test=1;
                                %end;
                                %if &_test ^=1 %then %do;                        
                                    %if %sysevalf(%superq(by&z)=,boolean) %then %do;
                                        %put ERROR: (Model &z: CLASSORDER) Number &i was not found in the CLASSORDER list;
                                        %put ERROR: (Model &z: CLASSORDER): Each number from 1 to maximum number of levels in CLASS variable %qupcase(%superq(class&z)) (%sysfunc(strip(%superq(nclass_&b._&z)))) must be represented;
                                    %end;
                                    %else %do;
                                        %put ERROR: (Model &z BY group &b: CLASSORDER) Number &i was not found in the CLASSORDER list;
                                        %put ERROR: (Model &z BY group &b: CLASSORDER): Each number from 1 to maximum number of levels of CLASS variable %qupcase(%superq(class&z)) (%sysfunc(strip(%superq(nclass_&b._&z)))) within BY group (%superq(by_&z._&b)) must be represented;
                                    %end;
                                    %let nerror_run=%eval(&nerror_run+1);                            
                                %end;
                            %end;
                        %end;     
                    %end;
                    
                    
                    
                    
                    %local j nby_last;
                    %let j=%sysevalf(&b-1);
                    %if &b>1 %then %let nby_last=%eval(&nby_last+&&nclass_&j._&z);
                    %else %let nby_last=0;
                    update _tempdsn&z
                        set _class_order=case (_class_)
                            %do i=1 %to &&nclass_&b._&z;
                                when "%superq(class_&b._&z._&i)" then &i + &nby_last
                            %end;
                            else . end
                        where _by_order=&b;
                    %if %qupcase(&&method&z)=INVWTS %then %do;
                            select _by_order,min(_class_order) into :null,:minclass_&z separated by ','
                                from _tempdsn&z group by _by_order;

                        /**Get counts for each Class Level**/            
                        select %do i = 1 %to %superq(nclass_&b._&z);
                                    %if &i>1 %then %do; , %end;
                                    sum(ifn(_class_order=&i+%scan(%superq(minclass_&z),&b,%str(,))-1,1,0))
                                %end;
                                into %do i = 1 %to %superq(nclass_&b._&z);
                                        %if &i>1 %then %do; , %end;
                                        :nclass_&b._&z._&i separated by ' '
                                     %end;
                                from _tempdsn&z
                                where _multi_=1 and _by_order=&b;
                    %end;
                    /**Check if the provided class reference value is in the data**/
                    %local _classrefcheck _classreforder&b._&z;
                    %let _classrefcheck=;
                    %if %sysevalf(%qscan(%superq(classref&z),&b,`,m)^=,boolean) %then %do;
                        select distinct 1,_class_order into :_classrefcheck,:_classreforder&b._&z separated by ''
                            from _tempdsn&z
                            where strip(_class_)=strip("%qscan(%superq(classref&z),&b,`,m)") and _multi_=0 and _by_order=&b;
                        %if %superq(_classrefcheck) ^=1 %then %do;
                            %if %sysevalf(%superq(by&z)^=,boolean) %then %do;
                                /**If not in the dataset then throw error**/
                                %put ERROR: (Model &z BY group &b: CLASSREF): Indicated class reference value (%qscan(%superq(classref&z),&b,`,m)) does not exist in dataset (%superq(data));
                                %put ERROR: (Model &z BY group &b: CLASSREF): Class reference value must be an existing value with same case as formatted value in dataset;
                                %put ERROR: (Model &z BY group &b: CLASSREF): Options are: %superq(___classlevels_&b._&z);  
                                %let nerror_run=%eval(&nerror_run+1);
                            %end;
                            %else %do;
                                /**If not in the dataset then throw error**/
                                %put ERROR: (Model &z: CLASSREF): Indicated class reference value (%qscan(%superq(classref&z),&b,`,m)) does not exist in dataset (%superq(data));
                                %put ERROR: (Model &z: CLASSREF): Class reference value must be an existing value with same case as formatted value in dataset;
                                %put ERROR: (Model &z BY group &b: CLASSREF): Options are: %superq(___classlevels_&b._&z);  
                                %let nerror_run=%eval(&nerror_run+1);
                            %end;
                        %end;                                     
                    %end;  
                    %else %do;
                        select distinct _class_,_class_order into :classref&z separated by '',:_classreforder&b._&z separated by '' from _tempdsn&z 
                            where ^missing(_class_) and _by_order=&b and _multi_=0
                            having _class_=max(_class_);
                    %end;
                    %if &nerror_run=0 %then %do;
                        update _tempdsn&z
                            set _class_ref=case(_class_order)
                                when %superq(_classreforder&b._&z) then 0
                                else _class_order end
                                where _by_order=&b;
                    %end;
                %end;
                %do i = 1 %to %superq(nclass_&z);
                    %local class_&z._&i;
                %end;
                select distinct _class_order,_class_ into :null,:class_&z._1-
                    from _tempdsn&z;
            %end;
            %else %do;
                %let nclass_1_&z=1;
                %let nclass_&z=1;
                %let class_&z._1=%str(&&classdesc&z);
            %end;
            /**Check that CLASS has at least two levels in every group**/
            %if %sysevalf(%superq(class&z)^=,boolean) %then %do;
                %local _min_class_lvls&z;
                select min(n) into :_min_class_lvls&z separated by ''
                    from  (select _by_order,count(distinct _class_order) as n from 
                             _tempdsn&z group by _by_order);
                 %if %sysevalf(%superq(_min_class_lvls&z) lt 2,boolean) %then %do;
                    %if %sysevalf(%superq(by&z)^=,boolean) %then 
                        %put ERROR: (Model &z: CLASS) Cannot have less than 2 levels in any BY variable subgroup;
                    %else %put ERROR: (Model &z: CLASS) Cannot have less than 2 levels;
                    %let nerror_run=%eval(&nerror_run+1);       
                 %end;
            %end;
            /*** X-axis Maximum ***/
            %local _maxt;
            /**Select greatest time value**/
            select max(_time_) into :_maxt
                from _tempdsn&z;
            /**Check if missing**/
            %if %sysevalf(%superq(xmax&z)=,boolean)=1 %then %do;
                /**Set x-maximum to maximum time rounded to next number divisible by 5**/
                %let xmax&z=%sysfunc(sum(%sysfunc(ceil(%superq(_maxt))),5,-%sysfunc(mod(%sysfunc(ceil(%superq(_maxt)-%superq(xmin&z))),5))));
            %end;
            %else %if %sysfunc(notdigit(%sysfunc(compress(%superq(xmax&z),.-)))) > 0 %then %do;
                /**If character value then throw error**/
                %put ERROR: (Model &z: XMAX) Must be a numeric value (%qupcase(%superq(xmax&z)));
                %let nerror_run=%eval(&nerror_run+1);
            %end;               
            %else %if %sysevalf(%superq(xmax&z) le %superq(xmin&z),boolean) %then %do;
                /**Make sure maximum is not less or equal to than minimum**/
                %put ERROR: (Model &z: XMAX) Cannot be less than or equal to XMIN (%superq(xmax&z) vs. %superq(xmin&z));
                %let nerror_run=%eval(&nerror_run+1);                   
            %end;
            /*** X-axis Incrementation ***/
            %if %sysevalf(%superq(xincrement&z)=,boolean)=1 %then %do;
                /**If missing then set increment to be range cut into 5 pieces**/
                %let xincrement&z=%sysfunc(sum((%superq(xmax&z)-%superq(xmin&z))/5));
            %end;
            %else %if %sysfunc(notdigit(%sysfunc(compress(%superq(xincrement&z),.-)))) > 0 %then %do;
                /**If character value then throw error**/
                %put ERROR: (Model &z: XINCREMENT) Must be a numeric value (%qupcase(%superq(xincrement&z)));
                %let nerror_run=%eval(&nerror_run+1);
            %end;
            %else %if %sysevalf(%superq(xincrement&z) gt (%sysevalf(%superq(xmax&z)-%superq(xmin&z))),boolean) %then %do;
                /**Make sure increment cannot be greater than range**/
                %put ERROR: (Model &z: XINCREMENT) Cannot be greater than the difference between XMAX and XMIN (%superq(xincrement&z) vs. %sysfunc(sum(%superq(xmax&z),-%superq(xmin&z))));                 
                %let nerror_run=%eval(&nerror_run+1);
            %end;
        quit;  
        %if &nerror_run > 0 %then %goto errhandl2; 
        %if %sysevalf(&&nby_&z>1,boolean) %then %do;
            proc sort data=_tempdsn&z;
                by _multi_ _by_order;
            run;
        %end;
        
        /**Assign multiplicative constant for y-axis**/
        %local xmult_&z tfmt_&z;
        %if %qupcase(%superq(ytype&z))=PCT %then %let xmult_&z=100;
        %else %if %qupcase(%superq(ytype&z))=PPT %then %let xmult_&z=1;
        /**Assign formats for y-axis**/
        %if %sysevalf(%qupcase(%superq(kmestdigits&z))=AUTO,boolean) %then %do;
            %if %qupcase(%superq(ytype&z))=PCT %then %let tfmt_&z=5.1;
            %else %if %qupcase(%superq(ytype&z))=PPT %then %let tfmt_&z=4.2;
        %end;
        %else %let tfmt_&z=%sysevalf(12.&&kmestdigits&z);
        
        
        %local risklist_&z;
        %if %qupcase(%superq(risklocation&z))=TIMELIST %then %let risklist_&z=%superq(timelist&z);
        %else %let risklist_&z = %superq(risklist&z);
             
        /**Get survival times**/
        /**median, logrank test, time estimates**/
        /*Unadjusted Curves*/
        ods graphics on;
        proc lifetest data=_tempdsn&z alpha=&&alpha&z alphaqt=&&alpha&z            
            /*%if %sysevalf(%qupcase(&&method&z)=KM,boolean) %then %do;
                RMST (%if %sysevalf(%superq(tau&z)^=,boolean) %then %do; tau=&&tau&z %end; cl)
                RMTL (%if %sysevalf(%superq(tau&z)^=,boolean) %then %do; tau=&&tau&z %end; cl)
            %end;*/
            /**Set up dataset with time-point estimate numbers**/
            %if %sysevalf(%superq(timelist&z)=,boolean)=0 %then %do;
                timelist=%sysfunc(compress(%superq(timelist&z), "'"))
                %if %sysevalf(%qupcase(&&method&z)^=CIF,boolean) %then %do;
                    reduceout outs=_timelist (where=(_multi_=0))
                %end;
            %end;  
            %if %sysevalf(%qupcase(&&method&z)=CIF,boolean) %then %do;
                /**Set up dataset with survival numbers**/
                cifvar 
                %if %sysevalf(%qupcase(%superq(cifvar&z))=COUNT,boolean) %then %do;
                    error=aalen
                %end;
                %else %do;
                    error=delta
                %end;
                outcif=_surv (where=(_multi_=0 and (event>. or t1=0 or c1=1))
                    rename=(%if %sysevalf(%superq(class&z)=,boolean)=0 %then %do; _class_order=stratumnum %end;
                        _time_=t1 cif=s1 censored=c1 cif_lcl=lcl1 cif_ucl=ucl1 variance=vcif)
                    keep=_multi_ _by_order cif _time_ censored event atrisk cif_lcl cif_ucl variance
                        %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do; _class_order %end;)
            %end;     
            /***Set up dataset with survival estimates***/ 
            %if %sysevalf(%qupcase(&&method&z)^=CIF,boolean) %then %do;
                plot=(survival(cl
                    %if %sysevalf(%superq(risklist&z)=,boolean)=0  %then %do; 
                        atrisk= %sysfunc(compress(%superq(risklist_&z), "'"))
                    %end;))
            %end;   
            /**Set confidence interval type**/ 
            conftype=%superq(conftype&z);
            /**If class variable then split into groups**/
            %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do;
                strata _class_order
                %if (%qupcase(%superq(plotpval&z))=LOGRANK or %qupcase(%superq(plotpval&z))=LOGRANK_ONESIDED) and %sysevalf(%superq(strata&z)=,boolean) %then %do;
                    / logrank
                %end;
                %else %if %qupcase(%superq(plotpval&z))=WILCOXON and %sysevalf(%superq(strata&z)=,boolean) %then %do;
                    / wilcoxon
                %end;;
            %end;
            by _multi_ _by_order;
            /**Run Model**/
            time _time_ * _cens_(%superq(cen_vl&z)) %if %qupcase(&&method&z)=CIF %then %do; / failcode=&&ev_vl&z %end;;            
            
            ods output 
                /*%if %qupcase(&&method&z)=KM %then %do; 
                    rmst=_rmst rmtl=_rmtl 
                    %if %sysevalf(%superq(strata&z)=,boolean) and %sysevalf(%superq(class&z)^=,boolean) %then %do; 
                        rmsttest=_rmstp 
                    %end; 
                %end; *//**Restricted Means Survival**/
                %if %qupcase(&&method&z)^=CIF %then %do; censoredsummary=_sum /**Contains events/totals**/%end;
                %else %do; failuresummary=_sum (rename=(event=failed)) /**Contains events/totals**/%end; 
                
                %if %qupcase(&&method&z)^=CIF %then %do; quartiles=_quart (where=(percent=50 and _multi_=0)) /**Containts Medians**/%end;
                
                %if %sysevalf(%superq(class&z)=,boolean)=0 and %sysevalf(%superq(strata&z)=,boolean) %then %do;
                    %if %qupcase(&&method&z)^=CIF %then %do;
                        homtests=_ltest /**Contains Logrank and/or Wilcoxon test**/
                    %end;
                    %else %do;
                        graytest=_ltest /**Contains Logrank and/or Wilcoxon test**/
                    %end;
                %end;
                
                %if %qupcase(&&method&z)^=CIF %then %do;
                    survivalplot=_surv (where=(event>. and _multi_=0)
                        rename=(survival=s1 time=t1 censored=c1 sdf_lcl=lcl1 sdf_ucl=ucl1)
                        keep=_multi_ _by_order survival time censored event sdf_lcl sdf_ucl
                            %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do; stratum stratumnum %end;
                            )/**Set up dataset with survival numbers**/
                    %if %qupcase(&&method&z)=DIRECT %then %do;
                        survivalplot=_surv_adj2 (where=(event>. and _multi_=1)
                            rename=(survival=s1 time=t1 censored=c1 sdf_lcl=lcl1 sdf_ucl=ucl1)
                            keep=_by_order _multi_ survival time censored event sdf_lcl sdf_ucl
                                %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do; stratum stratumnum %end;)/**Set up dataset with survival numbers**/                    
                    %end;
                    %if %sysevalf(%superq(risklist&z)=,boolean)=0 %then %do;
                        survivalplot=_splot
                        (keep=_multi_ _by_order time tatrisk 
                              %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do; stratum stratumnum %end; 
                              atrisk event censored) /**Outputs dataset with patients-at-risk numbers**/
                    %end;
                %end;
                %else %if %sysevalf(%superq(timelist&z)=,boolean)=0 %then %do;
                    cif=_timelist (rename=(cif=survival cif_lcl=sdf_lcl cif_ucl=sdf_ucl))                
                %end;;
        run;  
        %if %qupcase(&&method&z)=CIF %then %do;
            %if %sysevalf(%superq(class&z)=,boolean) %then %do;
                data _surv;
                    set _surv;
                    stratumnum=1;
                run;
            %end;
            %if %sysevalf(%superq(risklist&z)^=,boolean) %then %do;       
                %local _cif_riskpoints_&z;
                /**Created from XX to XX by XX format**/
                %if %sysfunc(find(%superq(risklist_&z),to,i))>0 %then %do;
                    data _list_;            
                        do i = %superq(risklist_&z);
                            riskpoint=i;
                            output;
                        end;
                        drop i;
                    run;
                    proc sort data=_list_;
                        by riskpoint;
                    run;
                    proc sql noprint;
                        select distinct riskpoint into :_cif_riskpoints_&z separated by ',' from _list_;
                        drop table _list_;
                    quit;
                %end;
                /**Created from numeric list**/
                %else %do i = 1 %to %sysfunc(countw(%superq(risklist&z),%str( )));
                    %if &i=1 %then %let _cif_riskpoints_&z=%scan(%superq(risklist&z),&i,%str( ));
                    %else %let _cif_riskpoints_&z=&&_cif_riskpoints_&z,%scan(%superq(risklist&z),&i,%str( ));
                %end; 
                data _splot;
                    set _surv;
                    by _by_order stratumnum;
                    where _multi_=0;
                    
                    array _risktimes_ {%sysfunc(countw(%superq(_cif_riskpoints_&z),%str(,)))} (&&_cif_riskpoints_&z);
                    retain _ncens _nevents _lagncens _lagnevents;
                    
                    if first.stratumnum then  do;
                        _ncens=0;_nevents=0;
                        _lagncens=0;_lagnevents=0;
                        _count_=1;
                    end;
                    if _count_ le dim(_risktimes_) then do;
                        if t1 <= _risktimes_{_count_} then do;
                            if event>. then _nevents=_nevents+event;
                            if c1>. then _ncens=_ncens+c1;
                        end;
                        else do;
                            if event>. then _lagnevents=_lagnevents+event;
                            if c1>. then _lagncens=_lagncens+c1;
                        end;
                    end;
                    if _count_ le dim(_risktimes_) then do;
                        if t1 >=_risktimes_(_count_) then do until(x=1);
                            tatrisk=_risktimes_(_count_);
                            time=tatrisk;
                            call missing(event,c1);
                            output;
                            _nevents=_nevents+_lagnevents;
                            _ncens=_lagncens+_ncens;
                            _lagnevents=0;
                            _lagncens=0;
                            _count_+1;
                            if _count_>dim(_risktimes_) then x=1;
                            else if t1 < _risktimes_(_count_) then x=1;
                            else x=0;
                        end;
                        if last.stratumnum then do i=_count_ to dim(_risktimes_);
                            tatrisk=_risktimes_(_count_);
                            time=tatrisk;
                            atrisk=0;
                            call missing(event,c1);
                            output;
                            _nevents=_nevents+_lagnevents;
                            _ncens=_lagncens+_ncens;
                            _lagnevents=0;
                            _lagncens=0;
                        end;
                    end;
                    keep _multi_ tatrisk time stratumnum atrisk event c1 t1 _count_ _nevents _ncens _lagnevents _lagncens;
                    rename c1=censored;
                run;
            %end;
            %if %sysevalf(%superq(class&z)^=,boolean) %then %do;
                data _sum;
                    set _sum;
                    if stratum<1 then control_var='-';
                    else control_var='';
                    stratumnum=stratum;
                run;
            %end;
            proc sql;
            /**Update censor values for plot**/
                update _surv
                    set c1=ifn(c1=1,s1,.);
            /**Generate Median quartiles table**/
                create table _quart as  
                    select _multi_,_by_order,%if %sysevalf(%superq(class)^=,boolean) %then %do; stratumnum as _class_order,%end; 50 as percent,
                    min(ifn(s1>=0.5,t1,.)) as estimate,
                    min(ifn(ucl1>=0.5,t1,.)) as lowerlimit,
                    min(ifn(lcl1>=0.5,t1,.)) as upperlimit
                    from _surv
                    where _multi_=0
                    group by _multi_,_by_order
                    %if %sysevalf(%superq(class)^=,boolean) %then %do; ,_class_order,_class_order %end;;
            quit;
        %end;
        %else %if %sysevalf(%superq(class&z)^=,boolean) %then %do;
            proc sql;
                update _surv
                    set stratumnum=input(stratum,12.);
                alter table _surv   
                    drop stratum;
                %if %qupcase(&&method&z)=DIRECT %then %do;
                    update _surv_adj2
                        set stratumnum=input(stratum,12.);
                    alter table _surv_adj2  
                        drop stratum;                
                %end;
                %if %sysevalf(%superq(risklist&z)=,boolean)=0 %then %do;
                    update _splot
                        set stratumnum=input(stratum,12.);
                    alter table _splot  
                        drop stratum;
                %end;
            quit;
        %end;
        %if %sysevalf(%superq(risklist&z)=,boolean)=0 and %qupcase(&&method&z)^=CIF %then %do;
            /*proc sort data=_splot;
                by _multi_ _by_order %if %sysevalf(%superq(class&z)^=,boolean) %then %do; stratumnum %end;;
            run;*/
            data _splot;
                set _splot;
                by _multi_ _by_order %if %sysevalf(%superq(class&z)^=,boolean) %then %do; stratumnum %end;;
                %if %sysevalf(%superq(class&z)^=,boolean) %then %do;  if first.stratumnum then do; %end;
                %else %do; if _n_=1 then do; %end;
                    _ncens=0;_nevents=0;
                end;
                if event>. then _nevents=_nevents+event;
                if censored>. then _ncens=_ncens+1;
                retain _ncens _nevents;
                if ^missing(tatrisk) then output;
                drop event censored;
            run;
        %end;
        ods graphics off; 
        /**calculate stratified p-values**/
        %if %sysevalf(%superq(strata&z)^=,boolean) and %sysevalf(%superq(class&z)^=,boolean) %then %do;
            %local rm_check;
            %if %sysevalf(%qupcase(&&method&z)=KM,boolean) and %sysfunc(find(%superq(display&z),rmst_pval,i))>0 %then %do;
                proc sql noprint;
                    %local rm_max rm_min;
                    select max(n),min(n) into :rm_max separated by '',:rm_min separated by '' from
                        (select %do i = 1 %to %sysfunc(countw(%superq(strata&z),%str( ))); _strata_&i, %end; count(distinct _class_order) as n
                            from _tempdsn&z where _multi_=0 group by _multi_,_by_order %do i = 1 %to %sysfunc(countw(%superq(strata&z),%str( ))); ,_strata_&i %end;);
                quit;
                %let rm_check=%sysevalf(&rm_max=&rm_min,boolean);
                %if &rm_check^=1 %then %do;
                    %put WARNING: CLASS values must be the same across STRATA to calculate RMST p-value;
                    %put WARNING: RMST p-value will not be calculated;
                    %let display&z=%compbl(%sysfunc(tranwrd(%qupcase(&&display&z),RMST_PVAL,%str())));
                %end;
            %end;
            proc lifetest data=_tempdsn&z         
                /*%if %sysevalf(%qupcase(&&method&z)=KM,boolean) and %sysfunc(find(%superq(display&z),rmst_pval,i))>0 and &rm_check=1 %then %do;
                    RMST (%if %sysevalf(%superq(tau&z)^=,boolean) %then %do; tau=&&tau&z %end; cl)
                    RMTL (%if %sysevalf(%superq(tau&z)^=,boolean) %then %do; tau=&&tau&z %end; cl)
                %end;*/;
                by _multi_ _by_order;
                strata %do i = 1 %to %sysfunc(countw(%superq(strata&z),%str( ))); _strata_&i %end; / group=_class_order
                    %if %qupcase(%superq(plotpval&z))=LOGRANK or %qupcase(%superq(plotpval&z))=LOGRANK_ONESIDED %then %do;
                        test=logrank
                    %end;
                    %else %if %qupcase(%superq(plotpval&z))=WILCOXON %then %do;
                        test=wilcoxon
                    %end;;
                /**Run Model**/
                time _time_ * _cens_(%superq(cen_vl&z)) %if %qupcase(&&method&z)=CIF %then %do; / failcode=&&ev_vl&z %end;;                
                ods output 
                    /*%if %sysevalf(%qupcase(&&method&z)=KM,boolean) and %sysfunc(find(%superq(display&z),rmst_pval,i))>0 and &rm_check=1 %then %do;
                        rmstStratifiedtest=_rmstp
                    %end;*/
                    %if %qupcase(&&method&z)^=CIF %then %do;
                        homtests=_ltest /**Contains Stratified Logrank and/or Wilcoxon test**/
                    %end;
                    %else %do;
                        graytest=_ltest /**Contains Stratified Gray K-sample test**/
                    %end;; 
            run; 
        %end;
        /**Inverse Weights Adjusted Survival Curves**/
        %if %qupcase(&&method&z)=INVWTS %then %do;
            /**Create weights for adjustment**/
            data _adjwts_prep;
                set _tempdsn&z; 
                where _multi_=1;
                array _bys_ {&&nby_&z} (%do b = 1 %to &&nby_&z; %if &b>1 %then %do; , %end; &&nclass_&b._&z %end;);
                array _min_ {&&nby_&z} (&&minclass_&z );
                _order_=_n_;
                do i = 1 to _bys_(_by_order);
                    _set_=i;
                    if _class_order=i+_min_(_by_order)-1 then _adjevent_=1;
                    else _adjevent_=0;
                    output;
                end;
                drop i _bys_:;
            run;
            
            proc sort data=_adjwts_prep;    
                by _by_order _set_;
            run;
            
            proc logistic data=_adjwts_prep;
                by _by_order _set_;
                /***Set up class variable and class covariates***/
                class 
                    /**Class covariates**/
                    %if %sysevalf(%superq(classcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(classcov&z),%str( ))); 
                        _classcov_&i 
                    %end;
                    %if %sysevalf(%superq(strata&z)=,boolean)=0 %then %do i = 1 %to %sysfunc(countw(%superq(strata&z),%str( ))); _strata_&i %end;;
                        
                /**Run model**/
                model _adjevent_ (event='1') = 
                    %if %sysevalf(%superq(classcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(classcov&z),%str( ))); _classcov_&i %end; /**Class Covariates**/
                    %if %sysevalf(%superq(contcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(contcov&z),%str( ))); _contcov_&i %end; /**Numeric Covariates**/
                    /**Apply strata**/
                    %if %sysevalf(%superq(strata&z)=,boolean)=0 %then %do i = 1 %to %sysfunc(countw(%superq(strata&z),%str( ))); _strata_&i %end; ; 
                output out=_adjwts prob=_adjwts_;
            run;
            
            proc sort data=_adjwts;
                by _by_order _order_ _set_;
            run;

            data _adjwts;
                merge %do b = 1 %to &&nby_&z;
                        %do i = 1 %to &&nclass_&b._&z; 
                            _adjwts (where=(_set_=&i and _by_order=&b) rename=(_adjwts_=_adjwts&b._&i) 
                            %if &i>1 %then %do; keep=_set_ _by_order _order_ _adjwts_ %end;) 
                        %end;
                      %end;;
                by _by_order _order_;   
                %do b=1 %to &&nby_&z;
                    array _class_n_&b._ {&&nclass_&b._&z} (%do i = 1 %to &&nclass_&b._&z; 
                                                            %if &i>1 %then %do; , %end; 
                                                            %superq(nclass_&b._&z._&i) 
                                                        %end;);
                %end;
                array _min_ {&&nby_&z} (&&minclass_&z );
                %do b=1 %to &&nby_&z;
                    array _adjwts&b._ {&&nclass_&b._&z};
                %end;
                %do b = 1 %to &&nby_&z;
                    %if &b>1 %then %do; else %end;
                    if _by_order=&b then do;
                        do i = 1 to dim(_class_n_&b._);
                            if _class_order=i+_min_(&b)-1 then _adjwt_=(_class_n_&b._(i)/sum(of _class_n_&b._(*)))/_adjwts&b._(i);
                        end; 
                    end;
                %end;
                drop  _set_ _order_ _adjevent_ _level_ i _class_n_: _adjwts:;
            run;
            /**Weighted Survival Curves**/ 
            ods graphics on;
            proc lifetest data=_adjwts alpha=&&alpha&z alphaqt=&&alpha&z
                /**Set up dataset with time-point estimate numbers**/
                %if %sysevalf(%superq(timelist&z)=,boolean)=0 %then %do;
                    reduceout timelist=%sysfunc(compress(%superq(timelist&z), "'"))
                    outs=_timelistmv
                %end;     
                /***Set up dataset with patients-at-risk numbers***/ 
                plot=(survival(cl
                    %if %sysevalf(%superq(risklist&z)=,boolean)=0  %then %do; 
                        atrisk= %sysfunc(compress(%superq(risklist_&z), "'"))
                    %end;))     
                    /**Set confidence interval type**/ 
                    conftype=%superq(conftype&z);
                by _by_order;
                /**If class variable then split into groups**/
                strata _class_order
                    %if %qupcase(%superq(plotpvalmv&z))=LOGRANK and %sysevalf(%superq(strata&z)=,boolean) %then %do;
                        / logrank
                    %end;
                    %else %if %qupcase(%superq(plotpvalmv&z))=WILCOXON and %sysevalf(%superq(strata&z)=,boolean) %then %do;
                        / wilcoxon
                    %end;;
                /**Run Model**/
                time _time_ * _cens_(%superq(cen_vl&z));
                weight _adjwt_;
                ods output censoredsummary=_summv /**Contains events/totals**/
                    quartiles=_quartmv (where=(percent=50)) /**Containts Medians**/
                    %if %sysevalf(%superq(class&z)=,boolean)=0 and %sysevalf(%superq(strata&z)=,boolean) %then %do;
                        homtests=_ltestmv /**Contains Logrank and/or Wilcoxon test**/
                    %end;
                    survivalplot=_surv_adj (where=(event>.)
                        rename=(survival=s1 time=t1 censored=c1 sdf_lcl=lcl1 sdf_ucl=ucl1)
                        keep=_by_order survival time censored event sdf_lcl sdf_ucl
                            %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do; stratum stratumnum %end;);/**Set up dataset with survival numbers**/
            run;
            %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do; 
                proc sql;
                    update _surv_adj
                        set stratumnum=input(stratum,12.);
                quit;
            %end;
            ods graphics off; 
            /**calculate stratified p-values**/
            %if %sysevalf(%superq(strata&z)^=,boolean) %then %do;
                /**Weighted Survival Curves**/ 
                proc lifetest data=_adjwts alpha=&&alpha&z alphaqt=&&alpha&z;
                    strata %do i = 1 %to %sysfunc(countw(%superq(strata&z),%str( )));_strata_&i %end; / group=_class_order
                        %if %qupcase(%superq(plotpvalmv&z))=LOGRANK %then %do;
                            test=logrank
                        %end;
                        %else %if %qupcase(%superq(plotpvalmv&z))=WILCOXON %then %do;
                            test=wilcoxon
                        %end;;
                by _by_order ;
                /**Run Model**/
                time _time_ * _cens_(%superq(cen_vl&z));                
                weight _adjwt_;
                ods output homtests=_ltestmv; /**Contains Stratified Logrank and/or Wilcoxon test**/
                run; 
            %end;
            
                
            /**Causes the survival estimates to continue**/
            /**for each time value in the dataset**/
            data _surv_adj;
                set _surv_adj;
                if s1 = . then s1 = c1;
                drop event;
            run;  
        %end;
        %if %qupcase(&&method&z)=DIRECT %then %do;
            /**Adjusted Survival Curves -- Direct Adjustment Method**/
            /*Setup population to be used as the covariates dataset*/
            data _covs;
                set %do b=1 %to &&nby_&z;
                        %do j=1 %to %superq(nclass_&z);
                            _tempdsn&z (in=b&b.a&j where=(_multi_=1 and _by_order=&b))
                        %end;
                    %end;;
                %do b=1 %to &&nby_&z;
                    %do j=1 %to %superq(nclass_&z);       
                        %if &j > 1 %then %do; else %end; 
                        if b&b.a&j then _class_order=&j;
                    %end;
                %end;
            run; 
            ods graphics on;
            /*Run Direct Adjusted Model*/        
            proc phreg data=_tempdsn&z plots(cl) =survival ;
                by _by_order;
                strata _class_order;
                class %if %sysevalf(%superq(strata&z)=,boolean)=0 %then %do i = 1 %to %sysfunc(countw(%superq(strata&z),%str( ))); _strata_&i %end;
                    %if %sysevalf(%superq(classcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(classcov&z),%str( ))); _classcov_&i %end;;
                model _time_*_cens_(%superq(cen_vl&z))= 
                    %if %sysevalf(%superq(classcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(classcov&z),%str( ))); _classcov_&i %end; /**Class Covariates**/
                    %if %sysevalf(%superq(contcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(contcov&z),%str( ))); _contcov_&i %end; /**Numeric Covariates**/
                    %if %sysevalf(%superq(strata&z)=,boolean)=0 %then %do i = 1 %to %sysfunc(countw(%superq(strata&z),%str( ))); _strata_&i %end; /**Stratification Factors**/ / rl;
                baseline covariates=_covs out=_surv_adj (rename=(_time_=t1 _class_order=stratumnum)) survival=s1 lower=lcl1 upper=ucl1/ diradj cltype=%superq(conftype&z); 
                where missing(_class_order)=0 and _multi_=1
                    %if %sysevalf(%superq(classcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(classcov&z),%str( )));
                        and missing(_classcov_&i)=0
                    %end;
                    %if %sysevalf(%superq(contcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(contcov&z),%str( )));
                        and missing(_contcov_&i)=0
                    %end;;
            run; 
            /*Add blank censor variable to outputted dataset*/
            %local _n;
            data _surv_adj;
                set _surv_adj;
                c1=.;
            run;
            %if %sysevalf(%superq(risklist&z)=,boolean)=0 %then %do;
                data _splot_adj;
                    set _splot;
                    where _multi_=1;
                run;
            %end; 
                        
            proc sort data=_surv_adj nodup;
                by _by_order stratumnum t1;
            proc sort data=_surv_adj2 nodup;
                by _by_order stratumnum t1;
            run;
            data _surv_adj;
                merge _surv_adj (in=a drop=c1) _surv_adj2 (keep=_by_order stratumnum t1 c1 where=(^missing(c1)));
                by _by_order stratumnum t1;
                if a or ^missing(c1);
            run;
            data _surv_adj;
                set _surv_adj end=last;
                by _by_order stratumnum;
                if first.stratumnum then do;
                    lcl1=1;
                    ucl1=1;
                end;
                if missing(lcl1) and ^missing(_lcl) then lcl1=_lcl;
                if missing(ucl1) and ^missing(_ucl) then ucl1=_ucl;
                retain _s1  _lcl _ucl;
                if s1 = . and c1>. then do;
                    s1 = _s1;
                    c1=_s1;
                    lcl1=_lcl;
                    ucl1=_ucl;
                end;
                else if c1>. then c1=s1;
                _s1=s1;_lcl=lcl1;_ucl=ucl1;
                
                if last then do;
                    call symput('_n',strip(put(_n_+1,12.0)));
                end;
                drop _s1 _lcl _ucl;            
            run;  
            
            /*Create a quartiles and timelist dataset for adjusted curves*/   
            %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;            
                %local _timepoints_&z;
                /**Created from XX to XX by XX format**/
                %if %sysfunc(find(%superq(timelist&z),to,i))>0 %then %do;
                    data _list_;            
                        do i = %superq(timelist&z);
                            timepoint=i;
                            output;
                        end;
                        drop i;
                    run;
                    proc sort data=_list_;
                        by timepoint;
                    run;
                    proc sql noprint;
                        select distinct timepoint into :_timepoints_&z separated by ',' from _list_;
                        drop table _list_;
                    quit;
                %end;
                /**Created from numeric list**/
                %else %do i = 1 %to %sysfunc(countw(%superq(timelist&z),%str( )));
                    %if &i=1 %then %let _timepoints_&z=%scan(%superq(timelist&z),&i,%str( ));
                    %else %let _timepoints_&z=&&_timepoints_&z,%scan(%superq(timelist&z),&i,%str( ));
                %end; 
            %end;    

            data _quartmv (keep=_by_order stratumnum _class_order percent estimate lowerlimit upperlimit )
                %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                    _timelistmv (keep= stratumnum _by_order _class_order timelist _time_ survival sdf_lcl sdf_ucl)
                %end;;
                set _surv_adj (rename=(t1=_time_ s1=_surv_ lcl1=_lcl_ ucl1=_ucl_)) end=last;
                by _by_order stratumnum _time_;
                array times {&_n} _temporary_;
                array surv {&_n} _temporary_;
                array lcl {&_n} _temporary_;
                array ucl {&_n} _temporary_;
                %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                    array _timepoints {%sysfunc(countw(%superq(_timepoints_&z),%str(,)))} (&&_timepoints_&z);
                %end;
                if first._by_order then do;
                    _mult_=1;
                    if _surv_=1*_mult_ then _order_=1;
                    else if _surv_=0 then _order_=2;
                end;
                retain _mult_ _order_;
                drop _mult_ _order_;
                
                if first.stratumnum then do;
                    _count_=1;
                    call missing(of times(*),of surv(*), of lcl(*), of ucl(*));
                end;
                else _count_+1;
                
                times(_count_)=_time_;
                surv(_count_)=_surv_;
                lcl(_count_)=_lcl_;
                ucl(_count_)=_ucl_;
                
                if last.stratumnum then do;
                    _class_order=stratumnum;
                    /**Median Time-to-Event**/
                    estimate=.;lowerlimit=.;upperlimit=.;percent=50;sdf_lcl=.;sdf_ucl=.;
                    do i = 1 to _count_;
                        if _order_=1 and ^missing(surv(i)) then do;
                            if surv(i) le 0.5*_mult_ and estimate=. then estimate=times(i);
                            if lcl(i) le 0.5*_mult_ and lowerlimit=. then lowerlimit=times(i);
                            if ucl(i) le 0.5*_mult_ and upperlimit=. then upperlimit=times(i);
                        end;
                        else if _order_=2 and ^missing(surv(i)) then do;
                            if surv(i) ge 0.5*_mult_ and estimate=. then estimate=times(i);
                            if lcl(i) ge 0.5*_mult_ and lowerlimit=. then lowerlimit=times(i);
                            if ucl(i) ge 0.5*_mult_ and upperlimit=. then upperlimit=times(i);
                        end;
                    end;
                    output _quartmv;
                    /**Time-point values**/
                    %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                        do j = 1 to dim(_timepoints);
                            _tout=0;
                            survival=.;sdf_lcl=.;sdf_ucl=.;
                            if _timepoints(j) > times(_count_) then do;
                                do i = j to dim(_timepoints);
                                    survival=.;sdf_lcl=.;sdf_ucl=.;
                                    timelist=_timepoints(i);
                                    output _timelistmv;
                                end;
                                j=dim(_timepoints);
                            end;
                            else do i = _tout+1 to _count_;
                                if times(i) <= _timepoints(j) and ^missing(surv(i)) then do;
                                    survival=surv(i);
                                    sdf_lcl=lcl(i);
                                    sdf_ucl=ucl(i);
                                    _tout=i;
                                end;
                                if times(i) = _timepoints(j) or
                                    (times(i) < _timepoints(j) and times(i+1)>_timepoints(j)) then do;
                                    _time_=times(_tout);
                                    timelist=_timepoints(j);
                                    output _timelistmv;
                                end;
                                if times(i) ge _timepoints(j) then i=_count_;
                            end;
                        end;
                    %end;
                end;
            run;
        %end;
        /**Causes the survival estimates to continue**/
        /**for each time value in the dataset**/
        data _surv;
            set _surv;
            if s1 = . then s1 = c1;
            drop event;
        run;  
        %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do;
            /**hazard ratio, p-value**/
            proc phreg data=_tempdsn&z (drop=_class_order rename=(_class_ref=_class_order)) alpha=&&alpha&z;
                by _multi_ _by_order;
                /***Set up class variable and class covariates***/
                class 
                    /***Class Variable and Reference Group***/
                    _class_order (ref="0") 
                        /**Class covariates**/
                        %if %sysevalf(%superq(classcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(classcov&z),%str( ))); 
                            _classcov_&i 
                        %end;;
                /**Apply strata**/
                %if %sysevalf(%superq(strata&z)=,boolean)=0 %then %do;
                    strata %do i = 1 %to %sysfunc(countw(%superq(strata&z),%str( ))); _strata_&i %end;;
                %end;
                /**Run model**/
                model _time_ * _cens_(%superq(cen_vl&z)) = _class_order /**Class Variable**/
                    %if %sysevalf(%superq(classcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(classcov&z),%str( ))); _classcov_&i %end; /**Class Covariates**/
                    %if %sysevalf(%superq(contcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(contcov&z),%str( ))); _contcov_&i %end; /**Numeric Covariates**/
                    / rl ties=%superq(hrties&z) type3 (SCORE LR WALD)
                    %if %sysevalf(%qupcase(%superq(method&z))=CIF) %then %do;
                        eventcode=%superq(ev_vl&z)
                    %end;; /**P-values, Ties, and Confidence Bounds, Competing Risks Event Level**/
                ods output parameterestimates=_parm (where=(upcase(strip(parameter))=upcase(strip("_class_order")))) /***Hazard ratios and confidence intervals***/
                    modelanova=_t3 (where=(upcase(strip(effect))=upcase(strip("_class_order")))) /**P-values**/;
                %if %sysevalf(%qupcase(%superq(method&z))^=CIF,boolean) %then %do;
                    /**Outputs the betas for C-index calculations**/
                    output out=_xbeta  xbeta=xbeta;
                %end; 
            run;
            /**Get Interaction P-values**/
            %if %sysevalf(%superq(by&z)^=,boolean) %then %do;                
                proc phreg data=_tempdsn&z (drop=_class_order rename=(_class_ref=_class_order)) alpha=&&alpha&z;
                    by _multi_;
                    /***Set up class variable and class covariates***/
                    class 
                        /***Class Variable and Reference Group***/
                        _class_ _by_order
                            /**Class covariates**/
                            %if %sysevalf(%superq(classcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(classcov&z),%str( ))); 
                                _classcov_&i 
                            %end;;
                    /**Apply strata**/
                    %if %sysevalf(%superq(strata&z)=,boolean)=0 %then %do;
                        strata %do i = 1 %to %sysfunc(countw(%superq(strata&z),%str( ))); _strata_&i %end;;
                    %end;
                    /**Run model**/
                    model _time_ * _cens_(%superq(cen_vl&z)) = _class_|_by_order /**Class/BY Interaction Variable**/
                        %if %sysevalf(%superq(classcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(classcov&z),%str( ))); _classcov_&i %end; /**Class Covariates**/
                        %if %sysevalf(%superq(contcov&z)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(contcov&z),%str( ))); _contcov_&i %end; /**Numeric Covariates**/
                        / rl ties=%superq(hrties&z) type3 (SCORE LR WALD)
                        %if %sysevalf(%qupcase(%superq(method&z))=CIF) %then %do;
                            eventcode=%superq(ev_vl&z)
                        %end;; /**P-values, Ties, and Confidence Bounds, Competing Risks Event Level**/
                    ods output modelanova=_t3i (where=(find(effect,'*')>0)) /**P-values**/;
                run;            
            %end;
            
            %if %sysevalf(%qupcase(%superq(method&z))^=CIF,boolean) %then %do;
                data _xbeta;
                    set _xbeta (in=a)
                        /*%if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;
                            _xbetamv (in=b)
                        %end;*/;
                    rename _multi_=_bylevel_;
                    *if a then _bylevel_=0;
                    *else _bylevel_=1;
                run;
                /**Sorts by present stratification**/
                %if %sysevalf(%superq(strata&z)=,boolean)=0 %then %do;
                    proc sort data=_xbeta out=_strata;
                        by _bylevel_ _by_order
                        %do j=1 %to %sysfunc(countw(%superq(strata&z),%str( ))); 
                            _strata_&j
                        %end;;
                    run;
                %end;
                /**Counts number of stratification levels**/
                data _stratlevels;
                    /**Pulls if strata are present**/
                    %if %sysevalf(%superq(strata&z)^=,boolean) %then %do;
                        set _strata end=__last;
                        by _bylevel_ _by_order
                        %do j=1 %to %sysfunc(countw(%superq(strata&z),%str( ))); 
                            _strata_&j
                        %end;;
                        if first._bylevel_ then _strata_=0;
                        if first._strata_%sysfunc(countw(%superq(strata&z),%str( ))) then _strata_+1;
                    %end;
                    /**Pulls if no strata are present**/
                    %else %do;
                        set _xbeta end=__last;
                        by _bylevel_ _by_order;
                        _strata_=1;
                    %end;
                run;
                /**Updates censor variable to always be a 0/1 type for this portion of code**/
                proc sql;
                    update _stratlevels
                        set _cens_=case(_cens_)
                                when . then . /**missings stay missing - shouldn't happen**/
                                when %superq(cen_vl&z) then 0 /**Censor levels become 0**/
                            else 1 end;/**Other levels are considered events**/
                quit;
                /**Performs analysis by stratification levels**/
                /**Sorts by time and descending status**/
                proc sort data=_stratlevels out=__step1;
                    by _bylevel_ _by_order _strata_ _time_ descending _cens_;
                    /**Sorts by the distinct betas to get an ordered list**/
                proc sort data=__step1 (keep=_bylevel_ _by_order _strata_ xbeta) nodupkey out=__ordered_xbeta;
                    by _bylevel_ _by_order _strata_ xbeta;
                run;
                /**Counts the number of distinct betas**/
                data __ordered_xbeta;
                    set __ordered_xbeta;
                    by _bylevel_ _by_order _strata_;
                    if first._strata_ then n = 0;
                    n+1;
                run;
                /**Create macro variables for later**/
                proc sql noprint;
                    %local _null_ _null2_ n _uv;
                    select _bylevel_, _by_order,_strata_,count(*),count(distinct xbeta)
                        into :_null_,:_null_,:_null_,:n separated by ', ',:_uv separated by ', '
                        from __step1 group by _bylevel_,_by_order,_strata_;/**Number of patients, Number of distinct xbetas**/
                    create table __step2 as
                        select a.*,
                            (select distinct N from __ordered_xbeta
                                where a._by_order=_by_order and a._bylevel_=_bylevel_ and a._strata_=_strata_ and a.xbeta=xbeta) as __rank
                        from __step1 a; /**Ranks patients by xbeta**/
                quit;
                /**Creates binary tree for indexing patients**/
                /**This copies BTREE from survConcordance.fit**/
                proc sort data=__step2 out=_ranks_prep (keep=_bylevel_ _by_order _strata_) nodupkey;
                    by _bylevel_ _by_order _strata_;
                run;
                data _ranks;
                    set _ranks_prep;
                    by _bylevel_ _by_order _strata_;
                     
                    /**Create arrays and variables to match survConcordance.fit**/
                    array _uvs_ {%sysfunc(countw(%superq(_uv),%str(,)))} (&_uv);
                    array yet_to_do {%sysfunc(max(0,&_uv))};
                    array indx {%sysfunc(max(0,&_uv))};
                    array ranks {%sysfunc(max(0,&_uv))};
                    if first._strata_ then m+1;
                    
                    call missing(of yet_to_do(*),of indx(*),of ranks(*));
                    depth=floor(log2(_uvs_(m)));
                    start=2**depth;
                    lastrow_length=1+_uvs_(m)-start;
                    do i = 1 to _uvs_(m);
                        yet_to_do(i)=i;
                    end;
                    do i = 1 to lastrow_length;
                        indx(i)=1+2*(i-1);
                    end;
                    do i = 1 to _uvs_(m);
                        if ^missing(indx(i)) then do;
                            ranks(yet_to_do(indx(i)))=start+(i-1);
                        end;
                    end;
                    do i = 1 to dim(indx);
                        if ^missing(indx(i)) then do;
                            yet_to_do(indx(i))=.;
                        end;
                    end;
                    count=1;
                    do i = 1 to dim(yet_to_do);
                        if ^missing(yet_to_do(i)) then do;
                            yet_to_do(count)=yet_to_do(i);
                            count=count+1;
                            yet_to_do(i)=.;
                        end;
                    end;
                    do while(start>1);
                        start=int(start/2);
                        do i = 1 to dim(indx);
                            indx(i)=.;
                        end;
                        do i = 1 to start;
                            indx(i)=1+2*(i-1);
                        end;
                        do i = 1 to dim(indx);
                            if ^missing(indx(i)) then do;
                                ranks(yet_to_do(indx(i)))=start+(i-1);
                            end;
                        end;
                        do i = 1 to dim(indx);
                            if ^missing(indx(i)) then do;
                                yet_to_do(indx(i))=.;
                            end;
                        end;
                        count=1;
                        do i = 1 to dim(yet_to_do);
                            if ^missing(yet_to_do(i)) then do;
                               yet_to_do(count)=yet_to_do(i);
                               count=count+1;
                               yet_to_do(i)=.;
                            end;
                        end;
                    end;
                    do i = 1 to dim(ranks);
                        /**Outputs indexes for the ranks**/
                        *call symput('_r'||strip(put(i,12.)),strip(put(ranks(i)-1,12.)));
                        ranks(i)=ranks(i)-1;
                    end;
                                            
                    keep _bylevel_ _by_order _strata_ ranks:;
                run;
                data __step2;
                    merge __step2 _ranks;
                    by _bylevel_ _by_order _strata_;
                run;
                /**Performs the C-index calculations**/
                /**Copies the function docount from survConcordance.fit**/
                data _step3;
                    set __step2 end=last;
                    by _bylevel_ _by_order _strata_;
                    /**Initializes arrays and variables for calculations**/
                    /**Method is to save values into temporary arrays and 
                    then perform calculations on the last observation**/
                    array ranks {*} ranks:;
                    array _ns_ {%sysfunc(countw(%superq(n),%str(,)))} (&n);
                    array _uvs_ {%sysfunc(countw(%superq(_uv),%str(,)))} (&_uv);
                    array y {2,%sysfunc(max(0,&n))} _temporary_;
                    array wt {%sysfunc(max(0,&n))} _temporary_;
                    array indx {%sysfunc(max(0,&n))} _temporary_;
                    array count {5};
                    array twt (%sysevalf(2*%sysfunc(max(0,&_uv))));
                    array nwt (%sysevalf(2*%sysfunc(max(0,&_uv))));
                    if first._strata_ then do;
                        call missing(of y(*),of wt(*),of indx(*));
                        _nc_=0;
                        _nstrata_+1;
                    end;
                    _nc_+1;                            
                    n=_ns_(_nstrata_);
                    y(1,_nc_)=_time_;
                    y(2,_nc_)=_cens_;
                    wt(_nc_)=1;
                    i=1;
                    do while (__rank^=i);
                        i=i+1;
                    end;
                    indx(_nc_)=ranks(i);
                    ntree=_uvs_(_nstrata_);
                    vss=0;/**Initializes for standard error calculation**/
                    if last._strata_ then do;/**Starts calculations**/
                        do i = 1 to dim(twt);
                            twt(i)=0;nwt(i)=0;
                        end;
                        do i = 1 to dim(count);
                            count(i)=0;
                        end;
                        i=n-1;
                        do while(i ge 0);/**First While loop**/
                            ndeath=0;
                            if y[2,i+1]=1 then do;/**Event Section**/
                                j=i;
                                if j > 0 then do while(j>0 and y(2,j+1)=1 and y(1,j+1)=y(1,i+1));/**Start J Loop**/
                                    ndeath=ndeath+wt(j+1);
                                    index=indx(j+1);
                                    k=i;
                                    if k > j then do while(k>j);/**Adds to ties**/
                                        count(4)=count(4)+wt(j+1)*wt(k+1);
                                        k=k-1;
                                    end;
                                    else k=k-1;
                                    count(3)=count(3)+wt(j+1)*nwt(index+1);/**Adds to time ties**/
                                    child=2*index+1;
                                    if child < ntree then
                                        count(1)=count(1)+wt(j+1)*twt(child+1);/**Adds to concordant pairs**/
                                    child=child+1;
                                    if child < ntree then
                                        count(2)=count(2)+wt(j+1)*twt(child+1);/**Adds to discordant pairs**/
                                    do while(index >0);
                                        parent=int((index-1)/2);
                                        if band(index,1) then
                                            count(2)=count(2)+wt(j+1)*(twt(parent+1)-twt(index+1));/**Adds to discordant pairs**/
                                        else count(1)=count(1)+wt(j+1)*(twt(parent+1)-twt(index+1));/**Adds to concordant pairs**/
                                        index=parent;
                                    end;
                                    j=j-1;
                                end;/**Ends J loop**/
                                else if j=0 then do;/**Finishes J=0 level of original loop**/
                                    if y(2,j+1)=1 and y(1,j+1)=y(1,i+1) then do;/**Event Loop**/
                                        ndeath=ndeath+wt(j+1);
                                        index=indx(j+1);
                                        k=i;
                                        if k > j then do while(k>j);
                                            count(4)=count(4)+wt(j+1)*wt(k+1);/**Adds to ties**/
                                            k=k-1;
                                        end;
                                        else k=k-1;
                                        count(3)=count(3)+wt(j+1)*nwt(index+1);/**Adds to time ties**/
                                        child=2*index+1;
                                        if child < ntree then
                                            count(1)=count(1)+wt(j+1)*twt(child+1);/**Adds to concordant pairs**/
                                        child=child+1;
                                        if child < ntree then
                                            count(2)=count(2)+wt(j+1)*twt(child+1);/**Adds to discordant pairs**/
                                        do while(index >0);
                                            parent=int((index-1)/2);
                                            if band(index,1) then
                                                count(2)=count(2)+wt(j+1)*(twt(parent+1)-twt(index+1));/**Adds to discordant pairs**/
                                            else count(1)=count(1)+wt(j+1)*(twt(parent+1)-twt(index+1));/**Adds to concordant pairs**/
                                            index=parent;
                                        end;
                                        j=j-1;
                                    end;/**Ends event loop**/
                                end;/**Ends J=0 loop**/
                            end;/**Ends Event section**/
                            else j = i-1;
                            if i>j then do while (i>j);/**Sum of squares section**/
                                wsum1=0;
                                oldmean=twt(1)/2;
                                index=indx(i+1);
                                nwt(index+1)=nwt(index+1)+wt(i+1);
                                twt(index+1)=twt(index+1)+wt(i+1);
                                wsum2=nwt(index+1);
                                child=2*index + 1;
                                if child < ntree then wsum1=wsum1+twt(child+1);
                                do while(index > 0);
                                    parent=int((index-1)/2);
                                    twt(parent+1)=twt(parent+1)+wt(i+1);
                                    if ^(band(index,1)) then wsum1=wsum1+(twt(parent+1)-twt(index+1));
                                    index=parent;
                                end;
                                wsum3=twt(1) - (wsum1+wsum2);
                                lmean=wsum1/2;
                                umean=wsum1+wsum2+wsum3/2;
                                newmean=twt(1)/2;
                                myrank=wsum1+wsum2/2;
                                /**Adds to sum of squares**/
                                vss=vss+wsum1*(newmean+oldmean-2*lmean)*(newmean-oldmean);
                                vss=vss+wsum3*(newmean+oldmean+wt(i+1)-2*umean)*(oldmean-newmean);
                                vss=vss+wt(i+1)*(myrank-newmean)*(myrank-newmean);
                                i=i-1;
                            end;/**Ends sum of squares section**/
                            else i=i-1;
                            count(5)=count(5)+ndeath*vss/twt(1);/**Adds to variance**/
                        end;/**Ends first loop**/
                        concordant=count(1);
                        discordant=count(2);
                        ties=count(3);
                        ties_time=count(4);
                        std=2*sqrt(count(5));/**Calculates standard deviation**/
                        if concordant+discordant+ties > 0 then do;
                            c=(concordant+ties/2)/(concordant+discordant+ties);/**Calculates concordance**/
                            se=std/(2*sum(concordant,discordant,ties));/**Calculates standard error**/
                        end;
                        else do;
                            c=0;se=0;
                        end;
                        output;
                    end;
                    keep _bylevel_ _by_order _strata_ concordant discordant ties ties_time std c se;
                run;
                /**Puts all calculations by strata level into same dataset**/
                /**Split by current BY level**/
                /**Calculates overall C-index within current by level**/
                data _cindex;
                    set _step3;
                    by _bylevel_ _by_order _strata_;
                    drop _strata_;
                    array temp {5}  _temporary_;
                    if first._by_order then do i = 1 to 5;
                        temp(i)=0;
                    end;
                    if ^(first._by_order and last._by_order) then do;
                        temp(1)=temp(1)+concordant;
                        temp(2)=temp(2)+discordant;
                        temp(3)=temp(3)+ties;
                        temp(4)=temp(4)+ties_time;
                        temp(5)=temp(5)+std;
                        c=.;
                        se=.;
                        strata=_strata_;
                        output;
                        /**Calculates overall c-index**/
                        if last._by_order then do;
                            concordant=temp(1);
                            discordant=temp(2);
                            ties=temp(3);
                            ties_time=temp(4);
                            std=temp(5);
                            c=(concordant+ties/2)/(concordant+discordant+ties);
                            se=std/(2*sum(concordant,discordant,ties));
                            strata=.;
                            output;
                        end;
                        drop i;
                    end;
                    else do;
                        strata=.;
                        output;
                    end;
                run;
                    
                proc datasets nolist nodetails;
                    %if &debug=0 %then %do;
                        delete  _stratlevels _strata _xbeta _xbetamv __ordered_xbeta _ranks _ranks_prep
                            __step1 __step2 _step3 _ranks;
                    %end;
                quit;
            %end;
            %else %do;
                data _cindex;
                    strata=.;_bylevel_=.; _by_order=.;c=.;se=.;output;
                run;
            %end;
        %end;
        proc sql noprint;
            /**Creates a fresh table to insert analysis into for summary table if flagged**/
            /**If not flagged, inserts new anlysis into current table**/
            %local _model;
            ***Check if output table has been determined;
            %if (%sysevalf(%superq(out)=,boolean) or (%sysevalf(%superq(out)^=,boolean) and (&newtable=1 or %sysfunc(exist(&out))=0))) and &z = 1 %then %do;
                create table %if %sysevalf(%superq(out)=,boolean)=0 %then %do; &out %end;
                             %else %do; _summary %end;
                    (modelnum num 'Model Number',/**Tracks models to distinguish later in PROC REPORT**/
                    modeltype num, /**Tracks whether the data is KM or CIF**/
                    by_order num, /**Used for Ordering Values**/
                    by_label char(300), /**Used for By labels**/
                    by_level char(300), /**Labeling the BY Order value**/
                    classlevel num, /**Used for ordering Values**/
                    subind num 'Indentation Indicator',/**Determines if a row is indented or not in PROC REPORT**/
                    subtitle char(100) 'Factor Label',/**Sub-title for each model's analysis.  Generally a description of the class variable**/
                    title char(2000) 'Model Title',/**Title for each model's analysis**/
                    footnote char(2000) 'Model Footnote',/**Footnotes for each model's analysis**/
                    total num 'Total',/**Number of patients**/
                    event num 'Events',/**Number of events**/
                    ev_n char(50) 'Events/N',/**Combination of events and patients.  Format events/patients**/
                    n_ev char(50) 'N (Events)',/**Combination of events and patients.  Format patients (events)**/
                    median char(50) 'Median',/**Median time-to-event with confidence interval**/
                    median_estimate num 'Median Estimate', /**Median time-to-event estimate**/
                    median_lcl num 'Median Lower Limit', /**Median time-to-event lower limit**/
                    median_ucl num 'Median Upper Limit', /**Median time-to-event upper limit**/
                    hr char(50) 'Hazard Ratio (CI)',/**Hazard ratio with confidence interval**/ 
                    hr_estimate num 'Hazard Ratio Estimate', /**Hazard Ratio estimate**/
                    hr_lcl num 'Hazard Ratio Lower Limit', /**Hazard Ratio lower limit**/
                    hr_ucl num 'Hazard Ratio Upper Limit', /**Hazard Ratio upper limit**/
                    covpval char(50) 'Covariate Level P-value',/**Covariate Level P-value**/
                    covpvaln num 'Covariate Level P-value',/**Covariate Level P-value**/
                    pval char(50) 'Displayed P-value',/**Chosen P-value to display for model**/
                    pvaln num 'Displayed P-value',/**Chosen P-value to display for model**/
                    pval_inter char(50) 'Interaction P-value',/**Interaction P-value between BY and CLASS**/
                    pvaln_inter num 'Interaction P-value',/**Interaction P-value between BY and CLASS**/
                    cindex char(50) 'C-index (CI)',/**C-index for univariate model**/
                    cindex_estimate num 'C-index Estimate', /**C-index estimate**/
                    cindex_lcl num 'C-index Lower Limit', /**C-index lower limit**/
                    cindex_ucl num 'C-index Upper Limit', /**C-index upper limit**/
                    tau num 'Tau (RMST)', /**Restricted Means Survival Time Tau**/
                    rmst char(50) 'RMST (CI)',/**Restricted Means Survival Time**/
                    rmst_estimate num 'RMST Estimate', /**RMST estimate**/
                    rmst_lcl num 'RMST Lower Limit', /**RMST lower limit**/
                    rmst_ucl num 'RMST Upper Limit', /**RMST upper limit**/
                    rmtl char(50) 'RMTL (CI)',/**Restricted Means Time Lost**/
                    rmtl_estimate num 'RMTL Estimate', /**RMTL estimate**/
                    rmtl_lcl num 'RMTL Lower Limit', /**RMTL lower limit**/
                    rmtl_ucl num 'RMTL Upper Limit', /**RMTL upper limit**/
                    rmst_pvaln num 'RMST P-value',/**Restricted Means Survival Time P-value**/
                    rmst_pval char(50) 'RMST P-value',/**Restricted Means Survival Time P-value**/
                    totalmv num 'Total (Multivariate)',/**Number of patients**/
                    eventmv num 'Events (Multivariate)',/**Number of events**/
                    ev_nmv char(50) 'Events/N (Multivariate)',/**Combination of events and patients.  Format events/patients**/
                    n_evmv char(50) 'N (Events) (Multivariate)',/**Combination of events and patients.  Format patients (events)**/
                    medianmv char(50) 'Median (Multivariate)',/**Median time-to-event with confidence interval**/
                    medianmv_estimate num 'Median Estimate (Multivariate)', /**Median time-to-event estimate**/
                    medianmv_lcl num 'Median Lower Limit (Multivariate)', /**Median time-to-event lower limit**/
                    medianmv_ucl num 'Median Upper Limit (Multivariate)', /**Median time-to-event upper limit**/
                    hrmv char(50) 'Hazard Ratio (CI) (Multivariate)',/**Hazard ratio with confidence interval**/ 
                    hrmv_estimate num 'Hazard Ratio Estimate (Multivariate)', /**Hazard Ratio estimate**/
                    hrmv_lcl num 'Hazard Ratio Lower Limit (Multivariate)', /**Hazard Ratio lower limit**/
                    hrmv_ucl num 'Hazard Ratio Upper Limit (Multivariate)', /**Hazard Ratio upper limit**/
                    pvalnmv num 'Displayed P-value (Multivariate)',/**Chosen P-value to display for model**/
                    pvalmv char(50) 'Displayed P-value (Multivariate)',/**Chosen P-value to display for model**/
                    pval_intermv char(50) 'Adjusted Interaction P-value',/**Adjusted Interaction P-value between BY and CLASS**/
                    pvaln_intermv num 'Adjusted Interaction P-value',/**Adjusted Interaction P-value between BY and CLASS**/
                    covpvalnmv num 'Covariate Level P-value (Multivariate)',/**Covariate Level P-value**/
                    covpvalmv char(50) 'Covariate Level P-value (Multivariate)',/**Covariate Level P-value**/
                    cindexmv char(50) 'C-index (CI) (Multivariate)', /**C-index for multivariate model**/
                    cindexmv_estimate num 'C-index Estimate (Multivariate)', /**C-index estimate**/
                    cindexmv_lcl num 'C-index Lower Limit (Multivariate)', /**C-index lower limit**/
                    cindexmv_ucl num 'C-index Upper Limit (Multivariate)', /**C-index upper limit**/
                    timepoint num 'KM Time-point', /*Event Time-point*/
                    timelist char(50) 'KM Estimate (CI)',/**Event time-point estimates**/
                    timelist_estimate num 'Event time-point Estimate', /**Event time-point estimate**/
                    timelist_lcl num 'Event time-point Lower Limit', /**Event time-point lower limit**/
                    timelist_ucl num 'Event time-point Upper Limit', /**Event time-point upper limit**/
                    timelistmv char(50) 'KM Estimate (CI) (Multivariate)',/**Event time-point estimates**/
                    timelistmv_estimate num 'Event time-point Estimate (Multivariate)', /**Event time-point estimate**/
                    timelistmv_lcl num 'Event time-point Lower Limit (Multivariate)', /**Event time-point lower limit**/
                    timelistmv_ucl num 'Event time-point Upper Limit (Multivariate)' /**Event time-point upper limit**/);
                %let _model=0;           
            %end;
            %else %if &z =1 %then %do;
                /**Grab maximum modelnum from previous table to increment upon**/
                select max(modelnum) into :_model from &out;
            %end;
            /**Create temporary table to insert analysis into before inserting into summary table**/
            create table _temptable like %if %sysevalf(%superq(out)=,boolean)=0 %then %do; &out %end;
                                         %else %do; _summary %end;;    
            /***Error Check if Enough Colors were Specified***/
            %if %sysevalf(%superq(color&z)=,boolean)=0 %then %do;
                %if %superq(nclass_&z) > %sysfunc(countw(%superq(color&z),%str( ))) and
                    %sysfunc(countw(%superq(color&z),%str( ))) > 1 %then %do;
                    %put ERROR: (Model &z: COLOR) Not enough line colors specified for number of class levels (%qupcase(%superq(color&z)) vs. %sysfunc(strip(%superq(nclass_&z))) Class Levels);
                    %let nerror_run=%eval(&nerror_run+1);
                %end;
            %end;
            %else %do;
                %put ERROR: (Model &z: COLOR) No line colors specified;
                %let nerror_run=%eval(&nerror_run+1);
            %end;
            /***Error Check if Enough Plot Line Patterns were Specified***/
            %if %sysevalf(%superq(pattern&z)=,boolean)=0 %then %do;
                %if %superq(nclass_&z) > %sysfunc(countw(%superq(pattern&z),%str( ))) and
                    %sysfunc(countw(%superq(pattern&z),%str( ))) > 1 %then %do;
                    %put ERROR: (Model &z: PATTERN) Not enough patterns specified for number of class levels (%qupcase(%superq(pattern&z)) vs. %sysfunc(strip(%superq(nclass_&z))) Class Levels);
                    %let nerror_run=%eval(&nerror_run+1);
                %end;
            %end;
            %else %do;
                %put ERROR: (Model &z: PATTERN) No line patterns specified;
                %let nerror_run=%eval(&nerror_run+1);
            %end;
            /**Adjusted**/
            %if %sysevalf(%superq(pattern_adjust&z)^=,boolean) and (%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT) %then %do;
                %if %superq(nclass_&z) > %sysfunc(countw(%superq(pattern_adjust&z),%str( ))) and
                    %sysfunc(countw(%superq(pattern_adjust&z),%str( ))) > 1 %then %do;
                    %put ERROR: (Model &z: PATTERN_ADJUST) Not enough patterns specified for number of class levels (%qupcase(%superq(pattern_adjust&z)) vs. %sysfunc(strip(%superq(nclass_&z))) Class Levels);
                    %let nerror_run=%eval(&nerror_run+1);
                %end;
                %else %do j=1 %to %superq(nclass_&z);
                    %local pattern_adjust&j._&z;
                    %if %sysfunc(countw(%superq(pattern_adjust&z),%str( ))) > 1 %then %let pattern_adjust&j._&z=%scan(%superq(pattern_adjust),&j,%str( ));
                    %else %Let pattern_adjust&j._&z=%superq(pattern_adjust);
                %end;
            %end;
            %else %if (%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT) %then %do;
                %put ERROR: (Model &z: PATTERN_ADJUST) No line patterns specified;
                %let nerror_run=%eval(&nerror_run+1);
            %end;
            %if %sysevalf(%superq(plotci&z)=1,boolean) %then %do;
                /**If Plot CI background fill is enabled**/
                %if %sysevalf(%superq(plotcifill&z)=1,boolean) %then %do;
                    /***Error Check if Enough Colors were Specified for Confidence Bounds Fill***/
                    %if %sysevalf(%superq(plotcifillcolor&z)=,boolean)=0 %then %do;
                        %if %superq(nclass_&z) > %sysfunc(countw(%superq(plotcifillcolor&z),%str( ))) and
                            %sysfunc(countw(%superq(plotcifillcolor&z),%str( ))) > 1 %then %do;
                            %put ERROR: (Model &z: PLOTCIFILLCOLOR) Not enough line colors specified for number of class levels (%qupcase(%superq(plotcifillcolor&z)) vs. %sysfunc(strip(%superq(nclass_&z))) Class Levels);
                            %let nerror_run=%eval(&nerror_run+1);
                        %end;
                    %end;
                    %else %do;
                        %put ERROR: (Model &z: PLOTCIFILLCOLOR) No line colors specified;
                        %let nerror_run=%eval(&nerror_run+1);
                    %end;
                %end;
                /***Error Check if Enough Colors were Specified for Confidence Bounds Lines***/
                %if %sysevalf(%superq(plotcilinecolor&z)=,boolean)=0 %then %do;
                    %if %superq(nclass_&z) > %sysfunc(countw(%superq(plotcilinecolor&z),%str( ))) and
                        %sysfunc(countw(%superq(plotcilinecolor&z),%str( ))) > 1 %then %do;
                        %put ERROR: (Model &z: PLOTCILINECOLOR) Not enough line colors specified for number of class levels (%qupcase(%superq(plotcilinecolor&z)) vs. %sysfunc(strip(%superq(nclass_&z))) Class Levels);
                        %let nerror_run=%eval(&nerror_run+1);
                    %end;
                %end;
                %else %do;
                    %put ERROR: (Model &z: PLOTCILINECOLOR) No line colors specified;
                    %let nerror_run=%eval(&nerror_run+1);
                %end;
                /***Error Check if Enough Plot Line Patterns were Specified***/
                %if %sysevalf(%superq(color&z)=,boolean)=0 %then %do;
                    %if %superq(nclass_&z) > %sysfunc(countw(%superq(plotcilinepattern&z),%str( ))) and
                        %sysfunc(countw(%superq(plotcilinepattern&z),%str( ))) > 1 %then %do;
                        %put ERROR: (Model &z: PLOTCILINEPATTERN) Not enough patterns specified for number of class levels (%qupcase(%superq(plotcilinepattern&z)) vs. %sysfunc(strip(%superq(nclass_&z))) Class Levels);
                        %let nerror_run=%eval(&nerror_run+1);
                    %end;
                %end;
                %else %do;
                    %put ERROR: (Model &z: PLOTCILINEPATTERN) No line patterns specified;
                    %let nerror_run=%eval(&nerror_run+1);
                %end;
            %end;              
            %if &nerror_run > 0 %then %goto errhandl2;
            
            
            
            %if %superq(sreverse&z) = 1 and %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                update _timelist
                    set survival=1-survival,
                    sdf_lcl=1-sdf_ucl,
                    sdf_ucl=1-sdf_lcl;
                %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do;
                    update _timelist_adj
                        set survival=1-survival,
                        sdf_lcl=1-sdf_ucl,
                        sdf_ucl=1-sdf_lcl;
                %end;
            %end;
            %local _ntimelist&z;
            %if %sysevalf(%superq(timelist&z)=,boolean)=0 %then %do;
                %local _tl_times&z;
                select distinct timelist into :_tl_times&z separated by '|'
                    from _timelist;
                %let _ntimelist&z=%sysfunc(countw(%superq(_tl_times&z),|));
            %end;
            %else %let _ntimelist&z=1;
            /**Insert analysis into final summary dataset**/
            /**Insert BY variable Header**/
            %if %sysevalf(%superq(by&z)^=,boolean) %then %do;
                insert into _temptable 
                    (modelnum,modeltype,by_order,by_level, classlevel,title,footnote,subind,subtitle,
                     pvaln_inter,pval_inter                 
                     %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do; 
                        ,pvaln_intermv,pval_intermv
                     %end;)
                    select sum(&z,&_model) as modelnum,
                    case (upcase("&&method&z"))
                        when 'KM' then &&sreverse&z
                        when 'CIF' then 2
                        when 'INVWTS' then 3
                        when 'DIRECT' then 4
                    else . end as modeltype,
                    0 as by_order,'' as by_level,
                    0 as classlevel,
                    strip("%superq(title&z)") as title,strip(tranwrd("%superq(footnote&z)",'`','^n')) as footnote,
                    (%sysevalf(%superq(title&z)^=,boolean) and %sysfunc(find(%superq(tabledisplay),title))>0) as subind,
                    strip(put(coalescec("%superq(classdesc&z)" %if %sysevalf(%superq(class&z)^=,boolean) %then %do; ,"%superq(label&z)" %end;,""),$100.))||" by %superq(bylabel&z)" as subtitle,
                    b.pval as pvaln_inter,
                    put(ifc(^missing(b.pval),
                        put(strip(tranwrd(put(b.pval,pvalue%sysevalf(6.&&pvaldigits&z)),'<.','<0.')),$30.)||
                        case (upcase("%superq(pval_inter&z)"))
                            when "LR" then '^{super #}' 
                            when "SCORE" then '^{super $}'
                            when "WALD" then '^{super +}' 
                        else '' end,""),$50.) as pval_inter
                    %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;                        
                        ,b2.pval as pvaln_intermv,
                        put(ifc(^missing(b2.pval),
                        put(strip(tranwrd(put(b2.pval,pvalue%sysevalf(6.&&pvaldigits&z)),'<.','<0.')),$30.)||
                            case (upcase("%superq(pval_intermv&z)"))
                                when "LR" then '^{super #}' 
                                when "SCORE" then '^{super $}'
                                when "WALD" then '^{super +}' 
                            else '' end,""),$50.) as pval_intermv
                    %end;
                from 
                    %if %qupcase(%superq(pval_inter&z))=LR or %qupcase(%superq(pval_inter&z))=SCORE %then %do;
                        (select *,prob&&pval_inter&z..chisq as pval from _t3i where _multi_=0) as b 
                    %end; 
                    %else %if %qupcase(%superq(pval_inter&z))=WALD %then %do;
                        (select *,probchisq as pval from _t3i where _multi_=0) as b
                    %end;
                    %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;
                        %if %qupcase(%superq(pval_inter&z))=LR or %qupcase(%superq(pval_inter&z))=SCORE %then %do;
                            left join (select *,prob&&pval_inter&z..chisq as pval from _t3i where _multi_=1) as b2 on b._multi_=(b2._multi_-1)
                        %end; 
                        %else %if %qupcase(%superq(pval_inter&z))=WALD %then %do;
                            left join (select *,probchisq as pval from _t3i where _multi_=1) as b2 on b._multi_=(b2._multi_-1)
                        %end;
                    %end;;            
            %end;
            insert into _temptable 
                (modelnum,modeltype,by_order,%if %sysevalf(%superq(by&z)^=,boolean) %then %do; by_level, %end; classlevel,title,footnote,subind,subtitle,
                 %if %sysevalf(%superq(class&z)=,boolean) %then %do; 
                    total,event,ev_n,n_ev,median_estimate,median_lcl,median_ucl,median
                    /*%if %qupcase(&&method&z)=KM %then %do;
                        ,tau,rmst_estimate,rmst_lcl,rmst_ucl,rmst,rmtl_estimate,rmtl_lcl,rmtl_ucl,rmtl 
                    %end;*/
                    %if %sysevalf(%superq(timelist&z)=,boolean)=0 %then %do;
                        ,timepoint,timelist_estimate,timelist_lcl,timelist_ucl,timelist
                    %end;
                 %end;
                 %else %do;
                    pvaln,pval, 
                    %if %qupcase(&&method&z)=KM  and %sysfunc(exist(_rmstp))^=0 %then %do;rmst_pvaln,rmst_pval, %end; cindex_estimate,cindex_lcl,cindex_ucl,cindex                    
                    %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;   
                        ,pvalnmv,pvalmv,
                        cindexmv_estimate,cindexmv_lcl,cindexmv_ucl,cindexmv
                    %end;
                 %end;
                    )
                select sum(&z,&_model) as modelnum,
                    case (upcase("&&method&z"))
                        when 'KM' then &&sreverse&z
                        when 'CIF' then 2
                        when 'INVWTS' then 3
                        when 'DIRECT' then 4
                    else . end as modeltype,
                    a._by_order as by_order,
                    %if %sysevalf(%superq(by&z)^=,boolean) and %sysevalf(%superq(class&z)^=,boolean) %then %do;  
                        case (a._by_order)
                            %do i = 1 %to &&nby_&z;
                                when &i then "%superq(by_&z._&i)"
                            %end;
                        else '' end as by_level,
                    %end;
                    0 as classlevel,
                    strip("%superq(title&z)") as title,strip(tranwrd("%superq(footnote&z)",'`','^n')) as footnote,
                    %sysevalf(%superq(by&z)^=,boolean)+(%sysevalf(%superq(title&z)^=,boolean) and %sysfunc(find(%superq(tabledisplay),title))>0) as subind,
                    %if %sysevalf(%superq(by&z)^=,boolean) and %sysevalf(%superq(class&z)^=,boolean) %then %do; 
                        calculated by_level 
                    %end;
                    %else %do;
                        put(coalescec("%superq(classdesc&z)" %if %sysevalf(%superq(class&z)^=,boolean) %then %do; ,"%superq(label&z)" %end;,""),$100.)                     
                    %end; as subtitle,
                    %if %sysevalf(%superq(class&z)=,boolean) %then %do; 
                        a.total as total, a.failed as event, 
                        put(ifc(nmiss(a.total,a.failed)=0,strip(put(a.failed,12.))||'/'||strip(put(a.total,12.)),''),$50.) as ev_n,
                        put(ifc(nmiss(a.total,a.failed)=0,strip(put(a.total,12.))||' ('||strip(put(a.failed,12.))||')',''),$50.) as n_ev,
                        b.estimate as median_estimate,b.lowerlimit as median_lcl,b.upperlimit as median_ucl,
                        case (b.estimate)
                                when . then 'NE'
                            else strip(put(b.estimate, %sysevalf(%sysevalf(12.&&mediandigits&z)))) end || ' (' ||
                            case (b.lowerlimit)
                                when . then 'NE'
                            else strip(put(b.lowerlimit, %sysevalf(12.&&mediandigits&z))) end || '-' ||
                            case (b.upperlimit)
                                when . then 'NE'
                            else strip(put(b.upperlimit, %sysevalf(12.&&mediandigits&z))) end || ')' as median
                            
                        /*%if %qupcase(&&method&z)=KM %then %do;
                           ,c.tau as tau,
                            c.estimate as rmst_estimate,c.lowerlimit as rmst_lcl,c.upperlimit as rmst_ucl,
                            strip(put(c.estimate,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||' ('||
                                strip(put(c.lowerlimit,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||'-'||
                                strip(put(c.upperlimit,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||')' as rmst,
                            c2.estimate as rmtl_estimate,c2.lowerlimit as rmtl_lcl,c2.upperlimit as rmtl_ucl,
                            strip(put(c2.estimate,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||' ('||
                                strip(put(c2.lowerlimit,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||'-'||
                                strip(put(c2.upperlimit,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||')' as rmtl  
                        %end;*/
                        %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                            ,d.timelist as timepoint,
                            %superq(xmult_&z)*d.survival as timelist_estimate,%superq(xmult_&z)*d.sdf_lcl as timelist_lcl,%superq(xmult_&z)*d.sdf_ucl as timelist_ucl,
                            case (&&listtimepoints&z)
                              when 1 then strip(put(d.timelist,12.))||" %superq(timedx&z): "
                            else '' end||
                                case(d.survival)
                                    when . then 'NE'
                                else strip(put(%superq(xmult_&z) *d.survival,%superq(tfmt_&z))) end || ' (' ||
                                case(d.sdf_lcl)
                                    when . then 'NE'
                                else strip(put(%superq(xmult_&z) *d.sdf_lcl,%superq(tfmt_&z))) end || '-' ||
                                case (d.sdf_ucl)
                                    when . then 'NE'
                                else strip(put(%superq(xmult_&z) *d.sdf_ucl,%superq(tfmt_&z))) end ||
                                %if %sysfunc(upcase(%superq(ytype&z))) = PCT %then %do; '%' || %end; ')' as timelist
                        %end;
                    %end;                    
                    %else %do;
                        b.pval as pvaln,put(strip(tranwrd(put(b.pval,pvalue%sysevalf(6.&&pvaldigits&z)),'<.','<0.')),$30.)||
                        case (upcase("%superq(plotpval&z)"))
                            when "LR" then '^{super #}' 
                            when "SCORE" then '^{super $}'
                            when "WALD" then '^{super +}' 
                            when "LOGRANK" then '^{super *}'
                            when "LOGRANK_ONESIDED" then '^{super **}'
                            when "WILCOXON" then '^{super @}' 
                            when "GRAY" then '^{super G}' 
                        else '' end as pval,
                        /*%if %qupcase(&&method&z)=KM and %sysfunc(exist(_rmstp))^=0 %then %do;
                            c.probchisq as rmst_pvaln,put(strip(tranwrd(put(c.probchisq,pvalue%sysevalf(6.&&pvaldigits&z)),'<.','<0.')),$30.) as rmst_pval,
                        %end;*/
                        d.c as cindex_estimate,d.c-probit(1-&&alpha&z/2)*d.se as cindex_lcl,d.c+probit(1-&&alpha&z/2)*d.se as cindex_ucl,
                        strip(put(d.c,%sysevalf(12.&&cindexdigits&z)))||' ('||
                            strip(put(d.c-probit(1-&&alpha&z/2)*d.se,%sysevalf(12.&&cindexdigits&z)))||'-'||
                            strip(put(d.c+probit(1-&&alpha&z/2)*d.se,%sysevalf(12.&&cindexdigits&z)))||')' as cindex
                        %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;                        
                            ,b2.pval as pvalnmv,put(strip(tranwrd(put(b2.pval,pvalue%sysevalf(6.&&pvaldigits&z)),'<.','<0.')),$30.)||                            
                            case (upcase("%superq(plotpvalmv&z)"))
                                when "LR" then '^{super #}' 
                                when "SCORE" then '^{super $}'
                                when "WALD" then '^{super +}' 
                                when "LOGRANK" then '^{super *}'
                                when "WILCOXON" then '^{super @}' 
                                when "GRAY" then '^{super G}' 
                            else '' end as pvalmv,
                            d2.c as cindexmv_estimate,d2.c-probit(1-&&alpha&z/2)*d2.se as cindexmv_lcl,d2.c+probit(1-&&alpha&z/2)*d2.se as cindexmv_ucl,
                            strip(put(d2.c,%sysevalf(12.&&cindexdigits&z)))||' ('||
                                strip(put(d2.c-probit(1-&&alpha&z/2)*d2.se,%sysevalf(12.&&cindexdigits&z)))||'-'||
                                strip(put(d2.c+probit(1-&&alpha&z/2)*d2.se,%sysevalf(12.&&cindexdigits&z)))||')' as cindexmv
                        %end;
                    %end;
                from 
                    %if %sysevalf(%superq(class&z)=,boolean) %then %do; 
                        (select * from _sum where _multi_=0) a left join 
                            (select * from _quart where _multi_=0) b on a._multi_=b._multi_ and a._by_order=b._by_order
                            /*%if %qupcase(&&method&z)=KM %then %do;
                                left join (select * from _rmst where _multi_=0) c on a._multi_=c._multi_ and a._by_order=c._by_order
                                left join (select * from _rmtl where _multi_=0) c2 on a._multi_=c2._multi_ and a._by_order=c2._by_order
                            %end;  */                          
                             %if %sysevalf(%superq(timelist&z)=,boolean)=0 %then %do;
                                left join (select * from _timelist where _multi_=0) d on a._multi_=d._multi_ and a._by_order=d._by_order
                            %end;
                    %end;
                    %else %do;
                        (select distinct _multi_,_by_order from _sum where _multi_=0) a left join 
                            %if %qupcase(%superq(plotpval&z))=LR or %qupcase(%superq(plotpval&z))=SCORE %then %do;
                                (select *,prob&&plotpval&z..chisq as pval from _t3 where _multi_=0 and upcase(effect)=upcase(strip("_class_order"))) as b on a._multi_=b._multi_ and a._by_order=b._by_order
                            %end; 
                            %else %if %qupcase(%superq(plotpval&z))=WALD %then %do;
                                (select *,probchisq as pval from _t3 where _multi_=0 and upcase(effect)=upcase(strip("_class_order"))) as b on a._multi_=b._multi_ and a._by_order=b._by_order
                            %end;
                            %else %if %qupcase(%superq(plotpval&z))=LOGRANK_ONESIDED %then %do;
                                (select *,probchisq/2 as pval from _ltest where _multi_=0) as b on a._multi_=b._multi_ and a._by_order=b._by_order
                            %end;
                            %else %do;
                                (select *,probchisq as pval from _ltest where _multi_=0) as b on a._multi_=b._multi_ and a._by_order=b._by_order
                            %end;
                            /*%if %qupcase(&&method&z)=KM and %sysfunc(exist(_rmstp))^=0 %then %do;
                                left join (select * from _rmstp where _multi_=0) as c on a._multi_=c._multi_ and a._by_order=c._by_order
                            %end;*/
                            left join (select * from _cindex where _bylevel_=0 and strata=.) d on a._multi_=d._bylevel_ and a._by_order=d._by_order
                            %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;
                                %if %qupcase(%superq(plotpvalmv&z))=LR or %qupcase(%superq(plotpvalmv&z))=SCORE %then %do;
                                    left join (select *,prob&&plotpvalmv&z..chisq as pval from _t3 where _multi_=1 and upcase(effect)=upcase(strip("_class_order"))) as b2 on a._multi_=(b2._multi_-1) and a._by_order=b2._by_order
                                %end;
                                %else %if %qupcase(%superq(plotpvalmv&z))=WALD %then %do;
                                    left join (select *,probchisq as pval from _t3 where _multi_=1 and upcase(effect)=upcase(strip("_class_order"))) as b2 on a._multi_=(b2._multi_-1) and a._by_order=b2._by_order
                                %end;
                                %else %do;
                                    left join (select *,probchisq as pval ,0 as _multi_ from _ltestmv) as b2 on a._multi_=b2._multi_ and a._by_order=b2._by_order
                                %end;
                                left join (select * from _cindex where _bylevel_=1 and strata=.) d2 on a._multi_=(d2._bylevel_-1) and a._by_order=d2._by_order
                            %end;
                    %end;;
  
                /*Inserts 1 row for each level of class variable*/ 
                %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do;
                    insert into _temptable 
                        (modelnum,modeltype,by_order,%if %sysevalf(%superq(by&z)^=,boolean) %then %do; by_level, %end; classlevel,title,footnote,subind,subtitle,
                         total,event,ev_n,n_ev,median_estimate,median_lcl,median_ucl,median,
                         /*%if %qupcase(&&method&z)=KM %then %do;
                            tau,rmst_estimate,rmst_lcl,rmst_ucl,rmst,rmtl_estimate,rmtl_lcl,rmtl_ucl,rmtl,   
                         %end;*/
                         %if %sysevalf(%superq(timelist&z)=,boolean)=0 %then %do;
                            timepoint,timelist_estimate,timelist_lcl,timelist_ucl,timelist,
                         %end;
                         hr_estimate,hr_lcl,hr_ucl,hr,covpvaln,covpval
                         %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;
                            ,totalmv,eventmv,ev_nmv,n_evmv,
                            %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do;
                                medianmv_estimate,medianmv_lcl,medianmv_ucl,medianmv,                  
                                %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                                    timelistmv_estimate,timelistmv_lcl,timelistmv_ucl,timelistmv,
                                %end;
                            %end;
                            hrmv_estimate,hrmv_lcl,hrmv_ucl,hrmv,covpvalnmv,covpvalmv
                         %end;)
                        select sum(&z,&_model) as modelnum,
                            case (upcase("&&method&z"))
                                when 'KM' then &&sreverse&z
                                when 'CIF' then 2
                                when 'INVWTS' then 3
                                when 'DIRECT' then 4
                            else . end as modeltype,
                            a._by_order as by_order,
                            %if %sysevalf(%superq(by&z)^=,boolean) %then %do;  
                                case (a._by_order)
                                    %do i = 1 %to &&nby_&z;
                                        when &i then "%superq(by_&z._&i)"
                                    %end;
                                else '' end as by_level,
                            %end;
                            a._class_order as classlevel,
                            strip("%superq(title&z)") as title,strip(tranwrd("%superq(footnote&z)",'`','^n')) as footnote,
                            1+%sysevalf(%superq(by&z)^=,boolean)+(%sysevalf(%superq(title&z)^=,boolean) and %sysfunc(find(%superq(tabledisplay),title))>0) as subind, 
                            case (a._class_order)
                                %do i=1 %to &&nclass_&z;
                                    when &i then "%superq(class_&z._&i)"
                                %end;
                            else '' end as subtitle,
                            a.total as total, a.failed as event, 
                            put(ifc(nmiss(a.total,a.failed)=0,strip(put(a.failed,12.))||'/'||strip(put(a.total,12.)),''),$50.) as ev_n,
                            put(ifc(nmiss(a.total,a.failed)=0,strip(put(a.total,12.))||' ('||strip(put(a.failed,12.))||')',''),$50.) as n_ev,
                            b.estimate as median_estimate,b.lowerlimit as median_lcl,b.upperlimit as median_ucl,
                            case (b.estimate)
                                    when . then 'NE'
                                else strip(put(b.estimate, %sysevalf(%sysevalf(12.&&mediandigits&z)))) end || ' (' ||
                                case (b.lowerlimit)
                                    when . then 'NE'
                                else strip(put(b.lowerlimit, %sysevalf(12.&&mediandigits&z))) end || '-' ||
                                case (b.upperlimit)
                                    when . then 'NE'
                                else strip(put(b.upperlimit, %sysevalf(12.&&mediandigits&z))) end || ')' as median,
                            /*%if %qupcase(&&method&z)=KM %then %do;
                                c.tau as tau,
                                c.estimate as rmst_estimate,c.lowerlimit as rmst_lcl,c.upperlimit as rmst_ucl,
                                strip(put(c.estimate,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||' ('||
                                    strip(put(c.lowerlimit,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||'-'||
                                    strip(put(c.upperlimit,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||')' as rmst,
                                c2.estimate as rmtl_estimate,c2.lowerlimit as rmtl_lcl,c2.upperlimit as rmtl_ucl,
                                strip(put(c2.estimate,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||' ('||
                                    strip(put(c2.lowerlimit,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||'-'||
                                    strip(put(c2.upperlimit,%sysevalf(%sysevalf(12.&&rmstdigits&z))))||')' as rmtl,   
                            %end;*/
                            %if %sysevalf(%superq(timelist&z)=,boolean)=0 %then %do;
                                d.timelist as timepoint,
                                %superq(xmult_&z)*d.survival as timelist_estimate,%superq(xmult_&z)*d.sdf_lcl as timelist_lcl,%superq(xmult_&z)*d.sdf_ucl as timelist_ucl,
                                case (&&listtimepoints&z)
                                  when 1 then strip(put(d.timelist,12.))||" %superq(timedx&z): "
                                else '' end||
                                    case(d.survival)
                                        when . then 'NE'
                                    else strip(put(%superq(xmult_&z) *d.survival,%superq(tfmt_&z))) end || ' (' ||
                                    case(d.sdf_lcl)
                                        when . then 'NE'
                                    else strip(put(%superq(xmult_&z) *d.sdf_lcl,%superq(tfmt_&z))) end || '-' ||
                                    case (d.sdf_ucl)
                                        when . then 'NE'
                                    else strip(put(%superq(xmult_&z) *d.sdf_ucl,%superq(tfmt_&z))) end ||
                                    %if %sysfunc(upcase(%superq(ytype&z))) = PCT %then %do; '%' || %end; ')' as timelist,
                            %end;
                            e.hazardratio as hr_estimate,e.hrlowercl as hr_lcl,e.hruppercl as hr_ucl,
                            put(ifc(^missing(e.parameter),
                                strip(put(e.hazardratio, %sysevalf(12.&&hrdigits&z)))|| ' (' ||
                                 strip(put(e.hrlowercl, %sysevalf(12.&&hrdigits&z))) || '-' ||
                                 strip(put(e.hruppercl, %sysevalf(12.&&hrdigits&z))) || ')',"%superq(refhrtext&z)"),$50.) as hr,
                            e.probchisq as covpvaln, put(ifc(^missing(e.parameter),
                                put(strip(tranwrd(put(e.probchisq,pvalue%sysevalf(6.&&pvaldigits&z)),'<.','<0.')),$30.)||'^{super +}',"%superq(refptext&z)"),$50.) as covpval
                            %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;
                               ,a2.total as totalmv, a2.failed as eventmv, 
                                put(ifc(nmiss(a2.total,a2.failed)=0,strip(put(a2.failed,12.))||'/'||strip(put(a2.total,12.)),''),$50.) as ev_nmv,
                                put(ifc(nmiss(a2.total,a2.failed)=0,strip(put(a2.total,12.))||' ('||strip(put(a2.failed,12.))||')',''),$50.) as n_evmv,
                                
                                %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do;
                                    b2.estimate as medianmv_estimate,b2.lowerlimit as medianmv_lcl,b2.upperlimit as medianmv_ucl,
                                    case (b2.estimate)
                                            when . then 'NE'
                                        else strip(put(b2.estimate, %sysevalf(%sysevalf(12.&&mediandigits&z)))) end || ' (' ||
                                        case (b2.lowerlimit)
                                            when . then 'NE'
                                        else strip(put(b2.lowerlimit, %sysevalf(12.&&mediandigits&z))) end || '-' ||
                                        case (b2.upperlimit)
                                            when . then 'NE'
                                        else strip(put(b2.upperlimit, %sysevalf(12.&&mediandigits&z))) end || ')' as medianmv,                
                                    %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;                                        
                                        %superq(xmult_&z)*d2.survival as timelistmv_estimate,%superq(xmult_&z)*d2.sdf_lcl as timelistmv_lcl,%superq(xmult_&z)*d2.sdf_ucl as timelistmv_ucl,
                                        case (&&listtimepoints&z)
                                          when 1 then strip(put(d2.timelist,12.))||" %superq(timedx&z): "
                                        else '' end||
                                            case(d2.survival)
                                                when . then 'NE'
                                            else strip(put(%superq(xmult_&z) *d2.survival,%superq(tfmt_&z))) end || ' (' ||
                                            case(d2.sdf_lcl)
                                                when . then 'NE'
                                            else strip(put(%superq(xmult_&z) *d2.sdf_lcl,%superq(tfmt_&z))) end || '-' ||
                                            case (d2.sdf_ucl)
                                                when . then 'NE'
                                            else strip(put(%superq(xmult_&z) *d2.sdf_ucl,%superq(tfmt_&z))) end ||
                                            %if %sysfunc(upcase(%superq(ytype&z))) = PCT %then %do; '%' || %end; ')' as timelistmv,
                                    %end;
                                %end;
                                e2.hazardratio as hr_estimatemv,e2.hrlowercl as hr_lclmv,e2.hruppercl as hr_uclmv,
                                put(ifc(^missing(e2.parameter),
                                    strip(put(e2.hazardratio, %sysevalf(12.&&hrdigits&z)))|| ' (' ||
                                     strip(put(e2.hrlowercl, %sysevalf(12.&&hrdigits&z))) || '-' ||
                                     strip(put(e2.hruppercl, %sysevalf(12.&&hrdigits&z))) || ')',"%superq(refhrtext&z)"),$50.) as hrmv,
                                e2.probchisq as covpvalnmv, 
                                put(ifc(^missing(e2.parameter),
                                    put(strip(tranwrd(put(e2.probchisq,pvalue%sysevalf(6.&&pvaldigits&z)),'<.','<0.')),$30.)||'^{super +}',"%superq(refptext&z)"),$50.) as covpvalmv
                            %end;
                        from 
                            (select * from _sum where _multi_=0 and missing(control_var)) a left join 
                                (select * from _quart where _multi_=0) b on a._multi_=b._multi_ and a._class_order=b._class_order and a._by_order=b._by_order
                                /*%if %qupcase(&&method&z)=KM %then %do;
                                    left join (select * from _rmst where _multi_=0) c on a._multi_=c._multi_ and a._class_order=c._class_order and a._by_order=c._by_order
                                    left join (select * from _rmtl where _multi_=0) c2 on a._multi_=c2._multi_ and a._class_order=c2._class_order and a._by_order=c2._by_order
                                %end;  */              
                                %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                                    left join (select * from _timelist where _multi_=0) d on a._multi_=d._multi_ and a._class_order=d._class_order and a._by_order=d._by_order
                                %end;
                                left join (select * from _parm where _multi_=0 and upcase(strip(parameter))=upcase(strip("_class_order"))) e on a._multi_=e._multi_ and a._by_order=e._by_order and strip(put(a._class_order,12.))=e.classval0
                                %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;
                                    left join (select * from _sum where _multi_=1 and missing(control_var)) a2 on a._class_order=a2._class_order and a._by_order=a2._by_order
                                    %if %qupcase(&&method&z)=INVWTS %then %do;
                                        left join (select * from _quartmv) b2 on  a._class_order=b2._class_order and a._by_order=b2._by_order               
                                        %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                                            left join (select * from _timelistmv) d2 on a._class_order=d2._class_order and a._by_order=d2._by_order and d.timelist=d2.timelist
                                        %end;
                                    %end;
                                    %else %if %qupcase(&&method&z)=DIRECT %then %do;
                                        left join (select * from _quartmv) b2 on  a._class_order=b2._class_order and a._by_order=b2._by_order                    
                                        %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                                            left join (select * from _timelistmv) d2 on a._class_order=d2._class_order and d.timelist=d2.timelist and a._by_order=d2._by_order
                                        %end;
                                    %end;
                                    left join (select * from _parm where _multi_=1 and upcase(strip(parameter))=upcase(strip("_class_order"))) e2 on strip(put(a._class_order,12.))=e2.classval0 and a._by_order=e2._by_order
                                %end;
                            ;
                %end;         
        quit;
                   
        %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
            proc sort data=_temptable;
                by by_order classlevel timepoint;
            data _temptable;
                set _temptable;
                by by_order classlevel timepoint;
                if ^first.classlevel then 
                    call missing(of total--cindexmv_ucl);
            run;
            data blah3;
            set _temptable;
            run;
        %end;
        /**Create a dataset for plotting**/
        data _plot_&z;
            /**Reverses Order if DESC=1**/
            %local i2;
            length by_order&z 8. %if %sysevalf(%superq(by&z)^=,boolean) %then %do; 
                                    by_level&z %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do; by_level_adj&z  %end; $300.
                                 %end;;
            merge
                /*One set of columns for each class variable level*/
                %do j = 0 %to %sysevalf((%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT));
                    %local sfx;
                    %if &j=0 %then %let sfx=;
                    %else %let sfx=_adj;
                    _surv&sfx. (
                        rename=(
                            _by_order=by_order&sfx.&z
                            %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do; stratumnum=strat&sfx._&z %end; /**Class Variable order**/
                            t1=t&sfx._&z /**Time Variable**/
                            c1=c&sfx._&z /**Censor Variable**/
                            s1=s&sfx._&z /**Survival Estimate Variable**/
                            lcl1=lcl&sfx._&z /**Survival 95% lower bound Variable**/
                            ucl1=ucl&sfx._&z /**Survival 95% upper bound Variable**/
                            %if %sysevalf(%qupcase(%superq(method&z))=CIF,boolean) %then %do;
                                vcif=vcif_&z /**Standard Error for CIF Function**/
                            %end;))
                %end;
                ;
            %do j = 0 %to %sysevalf((%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT));
                %local sfx sfx2;
                %if &j=0 %then %do;
                    %let sfx=;
                    %let sfx2=;
                %end;
                %else %do;
                    %let sfx2=%str( (Adjusted));
                    %let sfx=_adj;
                %end;
                length cl&sfx._&z $300.;
                
                %if %sysevalf(%superq(by&z)^=,boolean) %then %do b = 1 %to &&nby_&z;
                    %if &b>1 %then %do; else %end;
                    if by_order&sfx.&z=&b then by_level&sfx.&z="%superq(by_&z._&b)";
                %end;
                %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do;   
                    cln&sfx._&z=strat&sfx._&z;
                    drop strat&sfx._&z;
                    
                    %do i=1 %to &&nclass_&z;
                        %if &i>1 %then %do; else %end;
                        if cln&sfx._&z=&i then cl&sfx._&z="%superq(class_&z._&i)"||"&sfx2.";
                    %end;
                %end;
                %else %do;
                    cl_&z=strip("%superq(classdesc&z)");
                    cln_&z=1;
                %end;
            %end;
            /**If requested to plot 1-S instead of S**/
            %do j = 0 %to %sysevalf((%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT));
                %local sfx;
                %if &j=0 %then %let sfx=;
                %else %let sfx=_adj;
                %if %superq(sreverse&z) = 1 %then %do;
                    s&sfx._&z = 1-s&sfx._&z; /**Flip Survival Estimate Variable**/
                    c&sfx._&z = 1-c&sfx._&z; /**Flip Censor Marker Survival Estimate**/
                    _temp_=ucl&sfx._&z;/**Hold UCL Value**/
                    ucl&sfx._&z = 1-lcl&sfx._&z; /**Flip Survival Confidence Interval Upper Bound**/
                    lcl&sfx._&z = 1-_temp_; /**Flip Survival Confidence Interval Lower Bound**/
                %end;
            %end;
                
            label 
                %do j = 0 %to %sysevalf((%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT));
                    %local sfx sfx2;
                    %if &j=0 %then %do;
                        %let sfx=;
                        %let sfx2=;
                    %end;
                    %else %do;
                        %let sfx=_adj;
                        %let sfx2=%str(Adjusted );
                    %end;
                    t_&z="Time: Plot &z"
                    s_&z="SDF Estimate: Plot &z"
                    lcl_&z="SDF 95% Lower Bound: Plot &z"
                    ucl_&z="SDF 95% Upper Bound: Plot &z"
                    c_&z="Censor Estimate: Plot &z"
                    %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do;
                        cl_&z="Population Description &i: Plot &z"
                        cln_&z="Population Description Level &i: Plot &z"
                    %end;
                    %if %sysevalf(%qupcase(%superq(method&z))=CIF,boolean) %then %do;
                        vcif_&z="CIF Function Standard Error: Plot &z"
                    %end;
                %end;;
            drop _multi_ ; 
        run;
    
        /**Sets up dataset for patients at risk**/ 
        %if %sysevalf(%superq(risklist&z)=,boolean)=0 %then %do;
            %local partitle_&z;
            /**Determine if patients-at-risk header is requested**/
            %if %sysevalf(%superq(parheader&z)=%str(),boolean)=0 %then %let partitle_&z = 1;
            %else %let partitle_&z = 0;
            %local i2;
            data _riskplot;
                merge
                    /**Make one set of columns per class variable**/
                    %do j = 0 %to %sysevalf((%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT));
                        %local sfx;
                        %if &j=0 %then %let sfx=;
                        %else %let sfx=_adj;
                        %do i = 1 %to %superq(nclass_&z);
                            /**Reverses Order if DESC=1**/
                            %if %superq(desc&z)=1 %then %let i2=%sysfunc(abs(%sysfunc(sum(&i,-%superq(nclass_&z),-1))));
                            %else %let i2=&i;
                            _splot (rename=(
                                time=time&i.&sfx._&z /**Time Variable**/ 
                                atrisk=risk&i.&sfx._&z /**Number of Patients-at-Risk Variable**/
                                _ncens=cens&i.&sfx._&z /**Number of Patients Censored Variable**/
                                _nevents=event&i.&sfx._&z /**Number of Patients Censored Variable**/)
                                /**Grab different class values depending on class order requested**/
                                /**Stratumnum is numbered by PROC LIFETEST in the order the class variables are displayed**/
                                where=(_multi_=&j %if %sysevalf(%superq(class&z)^=,boolean) %then %do; and stratumnum=&i2 %end;))
                        %end;
                    %end;
                    ;
                /**Sets up Variables for the header in the Patients-at-Risk table when RISKLOCATION=INSIDE**/
                %if &&partitle_&z =1 %then %do;
                    if _n_ = 1 then do;
                    partitle_&z="&&parheader&z";
                    ncenstitle_&z="&&ncensheader&z";
                    neventstitle_&z="&&neventsheader&z";
                    end;
                    label partitle_&z="PAR Subheader: Plot &z"
                        ncenstitle_&z="Censored N Subheader: Plot &z"
                        neventstitle_&z="Events N Subheader: Plot &z";
                %end;
                
                %do j = 0 %to %sysevalf((%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT));
                    %local sfx;
                    %if &j=0 %then %let sfx=;
                    %else %let sfx=_adj;
                    %do i = 1 %to %superq(nclass_&z);
                        length atrisk&i.&sfx._&z ncens&i.&sfx._&z nevent&i.&sfx._&z $12.;
                        /**To be used in a BLOCKPLOT Value option**/
                        atrisk&i.&sfx._&z=strip(put(risk&i.&sfx._&z, 12.));
                        ncens&i.&sfx._&z=strip(put(cens&i.&sfx._&z, 12.));
                        nevent&i.&sfx._&z=strip(put(event&i.&sfx._&z, 12.));
                        drop risk&i.&sfx._&z cens&i.&sfx._&z event&i.&sfx._&z;
                    %end;
                %end;
                label
                    %do j = 0 %to %sysevalf((%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT));
                        %local sfx sfx2;
                        %if &j=0 %then %do;
                            %let sfx=;
                            %let sfx2=;
                        %end;
                        %else %do;
                            %let sfx=_adj;
                            %let sfx2=%str(Adjusted );
                        %end;
                        %do i = 1 %to %superq(nclass_&z);
                            %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do; 
                                time&i.&sfx._&z="&sfx2.PAR Time Class Level &i: Plot &z"
                                atrisk&i.&sfx._&z="&sfx2.PAR N Class Level &i: Plot &z"
                                ncens&i.&sfx._&z="&sfx2.Censored N Class Level &i: Plot &z"
                                nevent&i.&sfx._&z="&sfx2.Events N Class Level &i: Plot &z"
                            %end;
                            %else %do;
                                time&i._&z="PAR Time: Plot &z"
                                atrisk&i._&z="PAR N: Plot &z"
                                ncens&i._&z="Censored N: Plot &z"
                                nevent&i._&z="Events N: Plot &z"
                            %end;
                        %end;  
                    %end;
                    ;                          
                drop %if %sysevalf(%superq(class&z)^=,boolean) %then %do; stratumnum %end; tatrisk _multi_;   
            run;
            /**Merges patients-at-risk dataset with plot datset**/
            data _plot_&z;
                merge _plot_&z _riskplot;
            run;
        %end;
        
        /**Reference lines **/
        %if %sysevalf(%superq(reflines&z)^=,boolean) %then %do;
            %if %sysevalf(%qupcase(%superq(reflines&z))=TIMEPOINTS,boolean) and 
                %sysevalf(%superq(timelist&z)^=,boolean) %then %do;    
                /**Create at table of values for the plot**/
                proc sql noprint;
                        %if %sysevalf(%qupcase(%superq(reflinemethod&z))=FULL,boolean) %then %do;
                            create table _reflines_t as
                                select distinct timelist as ref_t_&z "Reference Lines Time Coordinates &z"
                                from (%if %qupcase(&&method&z)=KM or %qupcase(&&method&z)=CIF %then %do;
                                          select * from _timelist
                                      %end;
                                      %else %do;
                                          select timelist from _timelistmv 
                                          %if &&plot_unadjust&z=1 %then %do;
                                              OUTER UNION CORR 
                                              select timelist from _timelist
                                          %end;
                                      %end;);
                            create table _reflines_y as
                                select distinct %superq(xmult_&z)*survival as ref_y_&z "Reference Lines Y-Coordinate Coordinates &z"
                                from (%if %qupcase(&&method&z)=KM or %qupcase(&&method&z)=CIF %then %do;
                                          select * from _timelist
                                      %end;
                                      %else %do;
                                          select survival from _timelistmv 
                                          %if &&plot_unadjust&z=1 %then %do;
                                              OUTER UNION CORR 
                                              select survival from _timelist
                                          %end;
                                      %end;);
                            data _reflines;
                                merge _reflines_t _reflines_y;
                            run;
                        %end;
                        %else %if %sysevalf(%qupcase(%superq(reflinemethod&z))=DROP,boolean) %then %do;
                            create table _reflines as
                                select timelist as ref_t_&z "Reference Lines Time Coordinates &z",
                                       %superq(xmult_&z)*SURVIVAL as ref_y_&z "Reference Line Y-Coordinate &z"
                                from (%if %qupcase(&&method&z)=KM or %qupcase(&&method&z)=CIF %then %do;
                                          select * from _timelist
                                      %end;
                                      %else %do;
                                          select timelist,survival from _timelistmv 
                                          %if &&plot_unadjust&z=1 %then %do;
                                              OUTER UNION CORR 
                                              select timelist,survival from _timelist
                                          %end;
                                      %end;);
                        %end;
                quit;
            %end;
            %else %if %sysevalf(%qupcase(%superq(reflines&z))=MEDIANS,boolean) %then %do;    
                /**Create at table of values for the plot**/
                proc sql noprint;
                    create table _reflines as
                        %if %sysevalf(%qupcase(%superq(reflinemethod&z))=FULL,boolean) %then %do;
                            select distinct estimate as ref_t_&z "Reference Lines Time Coordinates &z",
                                   %superq(xmult_&z)*0.50 as ref_y_&z "Reference Line Y-Coordinate &z"
                            from (%if %qupcase(&&method&z)=KM or %qupcase(&&method&z)=CIF %then %do;
                                      select * from _quart
                                  %end;
                                  %else %do;
                                      select estimate from _quartmv 
                                      %if &&plot_unadjust&z=1 %then %do;
                                          OUTER UNION CORR 
                                          select estimate from _quart
                                      %end;
                                  %end;);
                        %end;
                        %else %if %sysevalf(%qupcase(%superq(reflinemethod&z))=DROP,boolean) %then %do;
                            select distinct estimate as ref_t_&z "Reference Lines Time Coordinates &z",
                                   %superq(xmult_&z)*0.50 as ref_y_&z "Reference Line Y-Coordinate &z"
                            from (%if %qupcase(&&method&z)=KM or %qupcase(&&method&z)=CIF %then %do;
                                      select * from _quart
                                  %end;
                                  %else %do;
                                      select estimate from _quartmv 
                                      %if &&plot_unadjust&z=1 %then %do;
                                          OUTER UNION CORR 
                                          select estimate from _quart
                                      %end;
                                  %end;);
                        %end;
                quit;
            %end;
            %if (%sysevalf(%qupcase(%superq(reflines&z))=TIMEPOINTS,boolean) and %sysevalf(%superq(timelist&z)^=,boolean)) or
                %sysevalf(%qupcase(%superq(reflines&z))=MEDIANS,boolean) %then %do; 
                /**Merges reference lines dataset with plot datset**/
                data _plot_&z;
                    merge _plot_&z _reflines;
                run;
            %end;
        %end;
    
        /**Determine which statistics are displayed in the plot**/
        /*Class Level Gridded Block*/
        %local _ndisplay_class_&z classcolumns hrcolumns;
        %let _ndisplay_class_&z=0;
        %let classcolumns=TOTAL|EVENT|MEDIAN|EV_N|N_EV|TOTALMV|EVENTMV|EV_NMV|N_EVMV|RMST|RMTL;
        %let hrcolumns=HR|HRMV|COVPVAL|COVPVALMV;
        /*Model Level Gridded Block*/
        %local _ndisplay_model_&z _ndisplay_mstats_&z _pvalcol_&z;
        %let _ndisplay_model_&z=0;
        %let _ndisplay_mstats_&z=0;
        /*Interaction p-value Gridded Block*/
        %local _ndisplay_interp_&z _ndisplay_istats_&z _ipval_&z;
        %let _ndisplay_interp_&z=0;
        %let _ndisplay_istats_&z=0;
        %if %sysevalf(%superq(display&z)=,boolean)=0 %then %do i = 1 %to %sysfunc(countw(%superq(display&z),%str( )));
            %if (%qupcase(%scan(%superq(display&z),&i,%str( )))=PVAL or %qupcase(%scan(%superq(display&z),&i,%str( )))=PVALMV)
                and &&model_stats_display&z=1 %then %let _pvalcol_&z=1;
            %if %sysevalf(%superq(by&z)^=,boolean) %then %let _ipval_&z=1;
        %end;
        %if %index(%qupcase(%superq(display&z)),LEGEND) > 0 %then %do;
            %let _ndisplay_class_&z=%sysevalf(&&_ndisplay_class_&z+1);
            %let _display_class&&_ndisplay_class_&z.._&z = LEGEND;
            %let _ndisplay_class_&z=%sysevalf(&&_ndisplay_class_&z+1);
            %let _display_class&&_ndisplay_class_&z.._&z = LEGEND2;
            %if &&_pvalcol_&z=1 or %index(%qupcase(%superq(display&z)),CINDEX) > 0 or
                (%sysevalf(%superq(by&z)=,boolean) and %index(%qupcase(%superq(display&z)),TABLECOMMENTS) > 0) or 
                (%sysevalf(%superq(by&z)=,boolean) and %superq(censormarkers&z)=1) %then %do;
                %let _ndisplay_model_&z=%sysevalf(&&_ndisplay_model_&z+1);
                %let _display_model&&_ndisplay_model_&z.._&z = LEGEND;
            %end;  
        %end;
        %if &&_ipval_&z=1 %then %do;
            %if &&_pvalcol_&z=1 
                or (&&model_stats_display&z=1 and %index(%qupcase(%superq(display&z)),CINDEX) > 0) %then %do;
                %let _ndisplay_model_&z=%sysevalf(&&_ndisplay_model_&z+1);
                %let _display_model&&_ndisplay_model_&z.._&z = STATS;
            %end;
            /**Interaction Test**/
            %if &&model_stats_display&z=1 or 1 or 
                (%index(%qupcase(%superq(display&z)),TABLECOMMENTS) > 0 and %sysevalf(%superq(tablecomments&z)^=,boolean)) %then %do;
                %let _ndisplay_interp_&z=%sysevalf(&&_ndisplay_interp_&z+1);
                %let _display_interp&&_ndisplay_interp_&z.._&z = STATS;
            %end;
            %if %superq(censormarkers&z)=1 %then %do;
                %let _ndisplay_interp_&z=%sysevalf(&&_ndisplay_interp_&z+1);
                %let _display_interp&&_ndisplay_interp_&z.._&z = CENSORS;
            %end;
        %end;
        %else %do;
            %if &&_pvalcol_&z=1 or %index(%qupcase(%superq(display&z)),TABLECOMMENTS) > 0 
                or (&&model_stats_display&z=1 and %index(%qupcase(%superq(display&z)),CINDEX) > 0) %then %do;
                %let _ndisplay_model_&z=%sysevalf(&&_ndisplay_model_&z+1);
                %let _display_model&&_ndisplay_model_&z.._&z = STATS;
            %end;
            %if %superq(censormarkers&z)=1 %then %do;
                %let _ndisplay_model_&z=%sysevalf(&&_ndisplay_model_&z+1);
                %let _display_model&&_ndisplay_model_&z.._&z = CENSORS;
            %end;
        %end;
        
        
        %if %sysevalf(%superq(display&z)=,boolean)=0 %then %do i = 1 %to %sysfunc(countw(%superq(display&z),%str( )));
            %local _display_current;
            %let _display_current=%qupcase(%scan(%superq(display&z),&i,%str( )));
            %if %sysevalf(%superq(class&z)^=,boolean) %then %do j = 1 %to %sysfunc(countw(&hrcolumns,|));
                %if &_display_current=%scan(&hrcolumns,&j,|) %then %do;
                    %let _ndisplay_class_&z=%sysevalf(&&_ndisplay_class_&z+1);
                    %let _display_class&&_ndisplay_class_&z.._&z = %scan(&hrcolumns,&j,|);
                %end;
            %end;
            %do j = 1 %to %sysfunc(countw(&classcolumns,|));
                %if &_display_current=%scan(&classcolumns,&j,|) %then %do;
                    %let _ndisplay_class_&z=%sysevalf(&&_ndisplay_class_&z+1);
                    %let _display_class&&_ndisplay_class_&z.._&z = %scan(&classcolumns,&j,|);
                %end;
            %end;
            %if &_display_current=TIMELIST %then %do;
                %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                    %let _ndisplay_class_&z=%sysevalf(&&_ndisplay_class_&z+1);
                    %let _display_class&&_ndisplay_class_&z.._&z = TIMELIST;
                    %let _ndisplay_class_&z=%sysevalf(&&_ndisplay_class_&z+%superq(listtimepoints&z));
                    %if %superq(listtimepoints&z) %then %do;
                        %let _display_class&&_ndisplay_class_&z.._&z = TIMEPOINTS;
                    %end;
                %end;                            
                %if  %qupcase(%superq(risklocation&z))= TIMELIST and %sysevalf(%superq(risklist&z)^=,boolean) %then %do;
                    %let _ndisplay_class_&z=%sysevalf(&&_ndisplay_class_&z+1);
                    %let _display_class&&_ndisplay_class_&z.._&z = RISKTABLE;
                %end;
            %end; 
            %do j = 1 %to %sysfunc(countw(PVAL|PVALMV,|));
                %if &_display_current=%scan(PVAL|PVALMV,&j,|) and %sysevalf(%superq(class&z)^=,boolean) %then %do;
                    %if &&model_stats_display&z^=1 %then %do;
                        %let _ndisplay_class_&z=%sysevalf(&&_ndisplay_class_&z+1);
                        %let _display_class&&_ndisplay_class_&z.._&z = %scan(PVAL|PVALMV,&j,|);
                    %end;
                    %else %do;
                        %let _ndisplay_mstats_&z=%sysevalf(&&_ndisplay_mstats_&z+1);
                        %let _display_mstats&&_ndisplay_mstats_&z.._&z = %scan(PVAL|PVALMV,&j,|);
                    %end;
                %end;
            %end;
            %do j = 1 %to %sysfunc(countw(PVAL_INTER|PVAL_INTERMV,|));
                %if &_display_current=%scan(PVAL_INTER|PVAL_INTERMV,&j,|) and %sysevalf(%superq(class&z)^=,boolean) %then %do;
                    %if &&model_stats_display&z=2 %then %do;
                        %let _ndisplay_class_&z=%sysevalf(&&_ndisplay_class_&z+1);
                        %let _display_class&&_ndisplay_class_&z.._&z = %scan(PVAL_INTER|PVAL_INTERMV,&j,|);
                    %end;
                    %else %do;
                        %let _ndisplay_istats_&z=%sysevalf(&&_ndisplay_istats_&z+1);
                        %let _display_istats&&_ndisplay_istats_&z.._&z = %scan(PVAL_INTER|PVAL_INTERMV,&j,|);
                    %end;
                %end;
            %end;
            %do j = 1 %to %sysfunc(countw(CINDEX|CINDEXMV,|));
                %if &_display_current=%scan(CINDEX|CINDEXMV,&j,|) and %sysevalf(%superq(class&z)^=,boolean) %then %do;
                    %if &&model_stats_display&z^=1 %then %do;
                        %let _ndisplay_class_&z=%sysevalf(&&_ndisplay_class_&z+1);
                        %let _display_class&&_ndisplay_class_&z.._&z = %scan(CINDEX|CINDEXMV,&j,|);
                    %end;
                    %else %do;
                        %let _ndisplay_mstats_&z=%sysevalf(&&_ndisplay_mstats_&z+1);
                        %let _display_mstats&&_ndisplay_mstats_&z.._&z = %scan(CINDEX|CINDEXMV,&j,|);
                    %end;
                %end;
            %end;
            %if &_display_current=RMST_PVAL and %sysevalf(%superq(class&z)^=,boolean) %then %do;
                %if &&model_stats_display&z^=1 %then %do;
                    %let _ndisplay_class_&z=%sysevalf(&&_ndisplay_class_&z+1);
                    %let _display_class&&_ndisplay_class_&z.._&z = RMST_PVAL;
                %end;
                %else %do;
                    %let _ndisplay_mstats_&z=%sysevalf(&&_ndisplay_mstats_&z+1);
                    %let _display_mstats&&_ndisplay_mstats_&z.._&z = RMST_PVAL; 
                %end;
            %end;
            %if &_display_current=TABLECOMMENTS and %sysevalf(%superq(tablecomments&z)^=,boolean) %then %do;
                %if &&_ipval_&z=1 %then %do;
                    %let _ndisplay_istats_&z=%sysevalf(&&_ndisplay_istats_&z+1);
                    %let _display_istats&&_ndisplay_istats_&z.._&z = TABLECOMMENTS; 
                    %let _ndisplay_istats_&z=%sysevalf(&&_ndisplay_istats_&z-1+%sysfunc(countw(%superq(tablecomments&z),`,m)));        
                %end;
                %else %do;
                    %let _ndisplay_mstats_&z=%sysevalf(&&_ndisplay_mstats_&z+1);
                    %let _display_mstats&&_ndisplay_mstats_&z.._&z = TABLECOMMENTS; 
                    %let _ndisplay_mstats_&z=%sysevalf(&&_ndisplay_mstats_&z-1+%sysfunc(countw(%superq(tablecomments&z),`,m)));        
                %end;
            %end;
        %end;
        proc sql noprint;
            /**saves values of each metric to macro variables**/
            create table _no_tl as
                select distinct * from _temptable (drop=timelist: timepoint) 
                order by classlevel;
            /**Hazard Ratios and Confidence Bounds**/
            %if %sysevalf(%superq(class&z)^=,boolean) %then %do;            
                %local by_level&z;
                select distinct by_order,by_level into :null,:by_level&z separated by '|'
                    from _no_tl
                    where ^missing(by_level);            
                %local hr&z;
                select hr into :hr&z separated by '|'
                    from _no_tl
                    where ^missing(hr);
                %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;
                    %local hrmv&z;
                    select hrmv into :hrmv&z separated by '|'
                        from _no_tl
                        where ^missing(hrmv);
                    /**Total number of patients**/
                    %local totalmv&z;
                    select totalmv into :totalmv&z separated by '|'
                        from _no_tl
                        where ^missing(totalmv);
                    /**Total number of events**/
                    %local eventmv&z;
                    select eventmv into :eventmv&z separated by '|'
                        from _no_tl
                        where ^missing(eventmv);
                    /**Formatted events/total**/
                    %local ev_nmv&z;
                    select ev_nmv into :ev_nmv&z separated by '|'
                        from _no_tl
                        where ^missing(ev_nmv);
                    /**Formatted total (events)**/
                    %local ev_nmv&z;
                    select strip(put(totalmv,12.))||' ('||strip(put(eventmv,12.))||')' into :n_evmv&z separated by '|'
                        from _no_tl
                        where ^missing(totalmv) and ^missing(eventmv);
                    %local pvalmv&z;                
                    select scan(pvalmv,1,'^') into :pvalmv&z separated by '|'
                        from _no_tl
                        where ^missing(pvalmv);
                    %local ipvalmv&z;                
                    select scan(pval_intermv,1,'^') into :ipvalmv&z separated by '|'
                        from _no_tl
                        where ^missing(pval_intermv);
                    %local cindexmv&z;                
                    select cindexmv into :cindexmv&z separated by '|'
                        from _no_tl
                        where ^missing(cindexmv);
                    %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do;
                        %local medianmv&z;
                        select medianmv into :medianmv&z separated by '|'
                            from _no_tl
                            where ^missing(medianmv);
                        /**Adjusted Time-point estimates**/
                        %local ntlmv_&z;
                        %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                            select count(distinct timelist)  into :ntlmv_&z
                                from _timelistmv;
                            %local timelistmvn&z timelistmvv&z;
                            select timelist,strip(put(timelist, best12.3)) || " %superq(timedx&z)"
                                into :timelistmvn&z separated by '|',:timelistmvv&z separated by '|'
                                from (select distinct timelist from _timelistmv);
                            %do k = 1 %to %superq(ntlmv_&z);
                                %local timelistmv_&k._&z;
                                select 
                                    case(survival)
                                        when . then 'NE'
                                    else strip(put(%superq(xmult_&z) *survival,%superq(tfmt_&z)))end || ' (' ||
                                    case(sdf_lcl)
                                        when . then 'NE'
                                    else strip(put(%superq(xmult_&z) *sdf_lcl,%superq(tfmt_&z))) end || '-' ||
                                    case (sdf_ucl)
                                        when . then 'NE'
                                    else strip(put(%superq(xmult_&z) *sdf_ucl,%superq(tfmt_&z))) end ||
                                %if %sysfunc(upcase(%superq(ytype&z))) = PCT %then %do; '%' || %end; ')'
                                into :timelistmv_&k._&z separated by '|'
                                from (
                                    %if %sysevalf(%superq(class&z)=,boolean) %then %do;
                                        select * from _timelistmv
                                    %end;
                                    %else %do i = 1 %to %superq(nclass_&z);
                                        select * from _timelistmv
                                        %if %sysevalf(%superq(class&z)=,boolean)=0 %then %do;
                                            where _class_order=&i
                                        %end;
                                        %if &i < %superq(nclass_&z) %then %do; OUTER UNION CORR %end;
                                    %end;)
                                where timelist=%scan(%superq(timelistmvn&z),&k,|);         
                            %end;
                        %end;
                    %end;
                %end;          
            %end;
            /**Median time-to-events and confidence bounds**/
            %local median&z;
            select median into :median&z separated by '|'
                from _no_tl
                where ^missing(median);
            /**Restricted Means Survival Time**/
            %local rmst&z;
            select rmst into :rmst&z separated by '|'
                from _no_tl
                where ^missing(rmst);
            /**Restricted Means Time Lost**/
            %local rmtl&z;
            select rmtl into :rmtl&z separated by '|'
                from _no_tl
                where ^missing(rmtl);
            /**Total number of patients**/
            %local total&z;
            select total into :total&z separated by '|'
                from _no_tl
                where ^missing(total);
            /**Total number of events**/
            %local event&z;
            select event into :event&z separated by '|'
                from _no_tl
                where ^missing(event);
            /**Formatted events/total**/
            %local ev_n&z;
            select ev_n into :ev_n&z separated by '|'
                from _no_tl
                where ^missing(ev_n);
            /**Formatted total (events)**/
            %local ev_n&z;
            select strip(put(total,12.))||' ('||strip(put(event,12.))||')' into :n_ev&z separated by '|'
                from _no_tl
                where ^missing(total) and ^missing(event);
            /**P-values**/
            %if %sysevalf(%superq(class&z)^=,boolean) %then %do;
                /**Model Level P-values**/
                %local pval&z;                
                select scan(pval,1,'^') into :pval&z separated by '|'
                    from _no_tl
                    where ^missing(pval);
                %local ipval&z;                
                select scan(pval_inter,1,'^') into :ipval&z separated by '|'
                    from _no_tl
                    where ^missing(pval_inter);
                /**RMST P-values**/
                %local rmst_pval&z;                
                select scan(rmst_pval,1,'^') into :rmst_pval&z separated by '|'
                    from _no_tl
                    where ^missing(rmst_pval);
                /**Univariate C-index**/  
                %local cindex&z;                
                select cindex into :cindex&z separated by '|'
                    from _no_tl
                    where ^missing(cindex);
                /**Covariate Level P-values**/                
                %local covpval&z;
                select scan(covpval,1,'^') into :covpval&z separated by '|'
                    from _no_tl
                    where ^missing(covpval);
                %if %sysevalf(%superq(classcov&z)^=,boolean) or %sysevalf(%superq(contcov&z)^=,boolean) %then %do;
                    %local covpvalmv&z;
                    select scan(covpvalmv,1,'^') into :covpvalmv&z separated by '|'
                        from _no_tl
                        where ^missing(covpvalmv);
                %end;
            %end;
            /**Time-point estimates**/
            %local ntl_&z;
            %if %sysevalf(%superq(timelist&z)^=,boolean) %then %do;
                select count(distinct timelist)  into :ntl_&z
                    from _timelist;
                %local timelistn&z timelistv&z;
                select timelist,strip(put(timelist, best12.3)) || " %superq(timedx&z)"
                    into :timelistn&z separated by '|',:timelistv&z separated by '|'
                    from (select distinct timelist from _timelist where _multi_=0);
                %do k = 1 %to %superq(ntl_&z);
                    %local timelist_&k._&z;
                    select 
                        case(survival)
                            when . then 'NE'
                        else strip(put(%superq(xmult_&z) *survival,%superq(tfmt_&z)))end || ' (' ||
                        case(sdf_lcl)
                            when . then 'NE'
                        else strip(put(%superq(xmult_&z) *sdf_lcl,%superq(tfmt_&z))) end || '-' ||
                        case (sdf_ucl)
                            when . then 'NE'
                        else strip(put(%superq(xmult_&z) *sdf_ucl,%superq(tfmt_&z))) end ||
                    %if %sysfunc(upcase(%superq(ytype&z))) = PCT %then %do; '%' || %end; ')'
                    into :timelist_&k._&z separated by '|'
                    from (
                        %if %sysevalf(%superq(class&z)=,boolean) %then %do;
                            select * from _timelist where _multi_=0
                        %end;
                        %else %do i = 1 %to %superq(nclass_&z);
                            select distinct timelist,survival,sdf_lcl,sdf_ucl from _timelist
                                where _multi_=0 and _class_order=&i
                            %if &i < %superq(nclass_&z) %then %do; OUTER UNION CORR %end;
                        %end;)
                    where timelist=%scan(%superq(timelistn&z),&k,|);         
                %end;
            %end;
            /**Patients-at-Risk outside summary table**/
            %if %sysevalf(%superq(risklist&z)^=,boolean) %then %do i = 1 %to %superq(nclass_&z);
                %local risklist_t_&i._&z risklist_v_&i._&z risklist_c_&i._&z risklist_e_&i._&z;
                select time&i._&z,atrisk&i._&z,ncens&i._&z,nevent&i._&z
                    into :risklist_t_&i._&z separated by '|',:risklist_v_&i._&z separated by '|',
                         :risklist_c_&i._&z separated by '|',:risklist_e_&i._&z separated by '|'
                    from _riskplot
                    where ^missing(atrisk&i._&z);
                %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do;
                    %local risklist_t_&i._adj_&z risklist_v_&i._adj_&z risklist_c_&i._adj_&z risklist_e_&i._adj_&z;
                    select time&i._adj_&z,atrisk&i._adj_&z,ncens&i._adj_&z,nevent&i._adj_&z
                        into :risklist_t_&i._adj_&z separated by '|',:risklist_v_&i._adj_&z separated by '|',
                             :risklist_c_&i._adj_&z separated by '|',:risklist_e_&i._adj_&z separated by '|'
                        from _riskplot
                        where ^missing(atrisk&i._adj_&z);
                %end;
            %end;
            /**Inserts values from temporary analysis summary into output dataset**/
            %if %sysevalf(%superq(out)=,boolean)=0 %then %do;
                insert into &out
                    select * from _temptable;
            %end;
            %else %do;
                insert into _summary
                select * from _temptable;
            %end;
        quit;
                
        /**Run-time Errors are sent there to delete temporary datasets before being sent to 
        errhandl, which stops the macro**/
        %errhandl2:
        proc datasets nodetails nolist;
           %if &debug=0 %then %do;
               delete _temptable _ltest _parm _parmmv _quart _sum _summv _surv _t3 _t3mv _t3i _reflines _reflines_t _reflines_y 
                    _riskplot _splot _timelist _tempdsn&z _tempcif _tempcif2 _score _variance _stat _cindex
                    _adjwts_prep _adjwts _surv_adj _ltestmv _quartmv _timelistmv _surv_adj2 _covs _splot_adj _dclass
                    _rmst _rmtl _rmstp _no_tl;
           %end;
        quit; 
        /**If errors occurred then throw message and end macro**/
        %if &nerror_run > 0 %then %do;
            %put ERROR: &nerror_run run-time errors listed;
            %put ERROR: Macro NEWSURV will cease;           
            %goto errhandl;
        %end;
    %end;/**Ends Analysis Loop**/
    
    /**Put all model plot datasets together for final plot dataset**/
    data _plot;
        merge
            %do z = 1 %to &nmodels;
                _plot_&z
            %end; ;
    run;
    
    /***Delete tails after xmax or below ymin for plotting***/
    data _plot;
        set _plot;
        
        %do z = 1 %to &nmodels;
            %if %sysevalf(%superq(class&z)^=,boolean) %then %do;
                length cl_&z._lag $300.;
                cl_&z._lag=lag1(cl_&z.);
                if cl_&z^=cl_&z._lag then flag_&z=0;
            %end;
            s_&z._lag=lag1(s_&z);
            lcl_&z._lag=lag1(lcl_&z);
            ucl_&z._lag=lag1(ucl_&z);
            retain flag_&z;
            if flag_&z^=1 then do;
                if t_&z > %superq(xmax&z) then do;
                    t_&z = %superq(xmax&z);
                    s_&z=s_&z._lag;
                    lcl_&z=lcl_&z._lag;
                    ucl_&z=ucl_&z._lag;
                    flag_&z=1;
                end;
                if s_&z*%superq(xmult_&z) lt %superq(ymin&z) then do;
                    s_&z = %superq(ymin&z)/%superq(xmult_&z);
                    flag_&z=1;
                end;
            end;
            else do;
                t_&z = .;
                s_&z = .;
                lcl_&z = .;
                ucl_&z = .;
                end;
            drop flag_&z s_&z._lag lcl_&z._lag ucl_&z._lag
                %if %sysevalf(%superq(class&z)^=,boolean) %then %do; cl_&z._lag %end;;
        %end;
    run;
               
    /**Creates template for Kaplan-Meier curve**/
    /*Calculate how many rows for each model*/
    %do z = 1 %to &nmodels;
        %local _rows&z _rowweights&z;
        %let _rows&z=1;
        %if %sysevalf(%superq(risklist&z)^=,boolean) and %sysevalf(%qupcase(%superq(risklocation&z))=BOTTOM,boolean) %then %do;
            /**Add in patients at risk rows**/
            %if %qupcase(&&method&z)=KM or %qupcase(&&method&z)=CIF or &&plot_unadjust&z=0 %then %let _rows&z=%sysevalf(&&_rows&z+%sysfunc(countw(%superq(pardisplay&z),%str( )))*%superq(nclass_&z));
            %else %let _rows&z=%sysevalf(&&_rows&z+%sysfunc(countw(%superq(pardisplay&z),%str( )))*2*%superq(nclass_&z));
            /**Add in rows for patients at risk labels when set to above**/
            %if %sysevalf(%qupcase(%superq(risklabellocation&z))=ABOVE,boolean) %then %do;
                %if %qupcase(&&method&z)=KM or %qupcase(&&method&z)=CIF or &&plot_unadjust&z=0 %then %let _rows&z=%sysevalf(&&_rows&z+%sysfunc(countw(%superq(pardisplay&z),%str( )))*%superq(nclass_&z));
                %else %let _rows&z=%sysevalf(&&_rows&z+%sysfunc(countw(%superq(pardisplay&z),%str( )))*2*%superq(nclass_&z));
            %end;
            /**Add in a row for each patients at risk header displayed**/
            %do j = 1 %to %sysfunc(countw(%superq(pardisplay&z),%str( )));
                %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&j,%str( )))=PAR,boolean) and
                    %sysevalf(%superq(parheader&z)^=,boolean) %then %let _rows&z=%sysevalf(&&_rows&z+1);
                %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&j,%str( )))=NCENS,boolean) and
                    %sysevalf(%superq(ncensheader&z)^=,boolean) %then %let _rows&z=%sysevalf(&&_rows&z+1);
                %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&j,%str( )))=NEVENTS,boolean) and
                    %sysevalf(%superq(neventsheader&z)^=,boolean) %then %let _rows&z=%sysevalf(&&_rows&z+1);
                %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&j,%str( )))=PAR_NEVENTS,boolean) and
                    %sysevalf(%superq(neventsheader&z)^=,boolean) and
                    %sysevalf(%superq(parheader&z)^=,boolean) %then %let _rows&z=%sysevalf(&&_rows&z+1);
                %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&j,%str( )))=PAR_NCENS,boolean) and
                    %sysevalf(%superq(ncensheader&z)^=,boolean) and
                    %sysevalf(%superq(parheader&z)^=,boolean) %then %let _rows&z=%sysevalf(&&_rows&z+1);
            %end;
            /**Add in a rows for each potential BY level**/
            %if %sysevalf(%superq(by&z)^=,boolean) %then %let _rows&z=%sysevalf(&&_rows&z+&&nby_&z*%sysfunc(countw(%superq(pardisplay&z),%str( ))));
        %end;
    %end;
    /*Calculate row weights for uniform height lattice*/
    %if &uniformheight=1 and &nmodels > 1 and %sysfunc(find(%superq(risklocation),bottom,i))>0 and 
        %sysevalf(%superq(risklist)^=,boolean) %then %do;
        %local _maxrows _maxrowweights;
        %let _maxrows=%superq(_rows1);
        %let _maxrowweights=%superq(riskrowweights1);
        %do i=2 %to &nmodels;
            %let _maxrows=%sysfunc(max(&_maxrows,%superq(_rows&i)));
            %let _maxrowweights=%sysfunc(max(&_maxrowweights,%superq(riskrowweights&i)));
        %end;           
    %end;
    proc template;
        /*Template for Summary Table: HTML, EXCEL, PDF*/
        define style _newsurvtable;
            parent=styles.rtf;
            style Table /
               color=black
               cellpadding = 0
               borderspacing = 0
               cellspacing=0
               frame = void
               rules = groups
               bordercollapse = separate
               borderleftstyle = none
               borderrightstyle = none
               bordertopstyle = none
               borderbottomstyle = none;
            style Header /
               color=black
               vjust=bottom
               backgroundcolor = white
               bordercollapse = separate
               borderleftstyle = none
               borderrightstyle = none
               bordertopstyle = solid
               borderbottomstyle = solid
               borderbottomcolor=black
               bordertopcolor=black
               bordertopwidth=0.1
               borderbottomwidth=0.1
               fontfamily="&TABLEHEADERFAMILY" 
               fontsize=&TABLEHEADERSIZE
               fontweight=&tableheaderweight;
            style Data /
               color=black
               backgroundcolor = white
               bordercolor = white
               borderstyle = none
               fontfamily="&TABLEDATAFAMILY" 
               fontsize=&TABLEDATASIZE
               fontweight=&tabledataweight
               vjust=top;
            class linecontent / 
                background=white 
                fontsize=&tablefootnotesize 
                color=black 
                fontfamily="&tablefootnotefamily"
                fontweight=&tablefootnoteweight;
        End;
        /*Template for Summary Table: POWERPOINT*/
        define style _newsurvtableppt;
            parent=styles.powerpointlight;
            class Header / 
                background=white 
               fontfamily="&TABLEHEADERFAMILY" 
               fontsize=&TABLEHEADERSIZE
               fontweight=&tableheaderweight
                color=black 
                vjust=bottom 
                borderstyle=solid 
                bordercolor=black 
                borderwidth=0.1 ;
            class Data / 
                background=white 
                color=black 
               fontfamily="&TABLEDATAFAMILY" 
               fontsize=&TABLEDATASIZE
               fontweight=&tabledataweight
                vjust=top
                borderstyle=hidden;
            class linecontent / 
                background=white 
                fontsize=&tablefootnotesize 
                color=black 
                fontfamily="&tablefootnotefamily"
                fontweight=&tablefootnoteweight
                borderstyle=solid 
                bordercolor=black 
                borderwidth=0.1;
            class Table / 
                color=black 
                cellpadding=0 
                borderspacing=0 
                cellspacing=0 
                frame=void 
                rules=rows 
                borderstyle=solid 
                bordercolor=black 
                borderwidth=0.1pt;
        End;
        /*Template for Graph*/
        define statgraph _km;
            begingraph / designheight=&height designwidth=&width    
                backgroundcolor=&background        
                %if %superq(transparent)=1 %then %do;
                    opaque=false 
                %end;
                /**Turns the border around the plot off if border=0**/
                %if %superq(border)=0 %then %do;
                    border=false 
                    %if %superq(transparent)=1 %then %do;
                        pad=0px    
                    %end;
                %end;;
                
                /**Create overall plot title**/
                %if %sysevalf(%superq(ovtitle)=,boolean)=0 %then %do i = 1 %to %sysfunc(countw(%superq(ovtitle),`,m));
                    entrytitle halign=&ovtitlealign "%scan(%superq(ovtitle),&i,`,m)" / 
                        textattrs=(color=&fontcolor weight=&ovtweight size=&ovtsize family="&ovtfamily" style=normal);
                %end;
                /**Create overall plot footnote**/
                %if %sysevalf(%superq(ovfootnote)=,boolean)=0 %then %do;
                    entryfootnote halign=&ovfootnotealign "&ovfootnote" / 
                        textattrs=(color=&fontcolor weight=&ovfnweight size=&ovfnsize family="&ovfnfamily" style=normal);
                %end;
                
                /*Creates Discrete Attribute Map*/
                %do z = 1 %to &nmodels;
                    discreteattrmap name="UA&z";
                        %do i = 1 %to %superq(nclass_&z);
                            value "&i" / 
                                lineattrs=(thickness=&&linesize&z
                                    %if %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                        color=%superq(color&z)
                                        %if %qupcase(%superq(pattern&z)) = AUTO %then %do;
                                            /**When all colors are the same, AUTO makes each pattern different**/
                                            pattern=&i
                                        %end;
                                        %else %do;
                                            %if %sysfunc(countw(%superq(pattern&z))) > 1 %then %do; pattern=%scan(%superq(pattern&z),&i) %end;
                                            %else %do; pattern=%superq(pattern&z) %end;
                                        %end;
                                    %end;
                                    %else  %do;
                                        %if %upcase(%superq(pattern&z)) = AUTO %then %do;
                                            /**When all colors are the different, AUTO makes each pattern solid**/
                                            pattern=solid
                                        %end;
                                        %else %do;
                                            %if %sysfunc(countw(%superq(pattern&z))) > 1 %then %do; pattern=%scan(%superq(pattern&z),&i) %end;
                                            %else %do; pattern=%superq(pattern&z) %end;
                                        %end;                            
                                        color=%scan(%superq(color&z), &i)                                        
                                    %end;)
                                markerattrs=(
                                            %if %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                color=%superq(color&z)
                                            %end;
                                            %else %do;
                                                color=%scan(%superq(color&z), &i)
                                            %end;
                                            symbol=plus
                                            size=%superq(symbolsize&z) weight=%superq(symbolweight&z));
                        %end;                    
                    enddiscreteattrmap;
                    discreteattrvar attrvar=_ua_&z var=cln_&z attrmap="UA&z";
                    %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do;
                        discreteattrmap name="A&z";
                            %do i = 1 %to %superq(nclass_&z);
                                value "&i" / 
                                    lineattrs=(thickness=&&linesize&z
                                        %if %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                            color=%superq(color&z)
                                            %if %qupcase(%superq(pattern_adjust&z)) = AUTO %then %do;
                                                /**When all colors are the same, AUTO makes each pattern different**/
                                                %if &&plot_unadjust&z=1 %then %do;
                                                    pattern=%sysevalf(&&nclass_&z + &i)
                                                %end;
                                                %else %do;
                                                    pattern=&i
                                                %end;
                                            %end;
                                            %else %do;
                                                %if %sysfunc(countw(%superq(pattern_adjust&z))) > 1 %then %do; pattern=%scan(%superq(pattern_adjust&z),&i) %end;
                                                %else %do; pattern=%superq(pattern_adjust&z) %end;
                                            %end;
                                        %end;
                                        %else  %do;
                                            %if %upcase(%superq(pattern_adjust&z)) = AUTO %then %do;
                                                /**When all colors are the different, AUTO makes each pattern solid**/
                                                %if &&plot_unadjust&z=1 %then %do; pattern=2 %end;
                                                %else %do; pattern=solid %end;
                                            %end;
                                            %else %do;
                                                %if %sysfunc(countw(%superq(pattern_adjust&z))) > 1 %then %do; pattern=%scan(%superq(pattern_adjust&z),&i) %end;
                                                %else %do; pattern=%superq(pattern_adjust&z) %end;
                                            %end;                            
                                            color=%scan(%superq(color&z), &i)                                        
                                        %end;)
                                    markerattrs=(
                                                %if %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                    color=%superq(color&z)
                                                %end;
                                                %else %do;
                                                    color=%scan(%superq(color&z), &i)
                                                %end;
                                                symbol=plus
                                                size=%superq(symbolsize&z) weight=%superq(symbolweight&z));
                            %end;                    
                        enddiscreteattrmap;
                        discreteattrvar attrvar=_a_&z var=cln_adj_&z attrmap="A&z";                    
                    %end;
                    /**Confidence Bounds if flagged**/
                    %if %superq(plotci&z) = 1 or 
                        (%superq(plotci&z)=2 and %sysevalf(%superq(class&z)=,boolean)) %then %do;
                        discreteattrmap name="UA_CI_&z";
                            %do i = 1 %to %superq(nclass_&z);
                                value "&i" / 
                                    fillattrs=(transparency=%superq(plotcifilltransparency&z) 
                                               %if %sysfunc(countw(%superq(plotcifillcolor&z))) = 1 %then %do;
                                                  color=%superq(plotcifillcolor&z)
                                               %end;
                                               %else %do;
                                                   color=%scan(%superq(plotcifillcolor&z), &i)
                                               %end;)
                                    lineattrs=(thickness=%superq(plotcilinesize&z)
                                                  %if %sysfunc(countw(%superq(plotcilinecolor&z))) = 1 %then %do;
                                                      color=%superq(plotcilinecolor&z)
                                                  %end;
                                                  %else %do;
                                                      color=%scan(%superq(plotcilinecolor&z), &i)
                                                  %end;
                                                  %if %sysfunc(countw(%superq(plotcilinepattern&z))) = 1 %then %do;
                                                      pattern=%superq(plotcilinepattern&z)
                                                  %end;
                                                  %else %do;
                                                      pattern=%scan(%superq(plotcilinepattern&z), &i)
                                                  %end;);
                              %end;
                        enddiscreteattrmap;
                        discreteattrvar attrvar=_ua_ci_&z var=cln_&z attrmap="UA_CI_&z";  
                        %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do;
                            discreteattrmap name="A_CI_&z";
                                %do i = 1 %to %superq(nclass_&z);
                                    value "&i" / 
                                        fillattrs=(transparency=%superq(plotcifilltransparency&z) 
                                                   %if %sysfunc(countw(%superq(plotcifillcolor&z))) = 1 %then %do;
                                                      color=%superq(plotcifillcolor&z)
                                                   %end;
                                                   %else %do;
                                                       color=%scan(%superq(plotcifillcolor&z), &i)
                                                   %end;)
                                        lineattrs=(thickness=%superq(plotcilinesize&z)
                                                      %if %sysfunc(countw(%superq(plotcilinecolor&z))) = 1 %then %do;
                                                          color=%superq(plotcilinecolor&z)
                                                      %end;
                                                      %else %do;
                                                          color=%scan(%superq(plotcilinecolor&z), &i)
                                                      %end;
                                                      %if %sysfunc(countw(%superq(plotcilinepattern&z))) = 1 %then %do;
                                                          pattern=%superq(plotcilinepattern&z)
                                                      %end;
                                                      %else %do;
                                                          pattern=%scan(%superq(plotcilinepattern&z), &i)
                                                      %end;);
                                  %end;
                            enddiscreteattrmap;
                            discreteattrvar attrvar=_a_ci_&z var=cln_adj_&z attrmap="A_CI_&z";
                          %end;  
                    %end;
                %end;
                %do z = 1 %to &nmodels;
                    /**Legend Spacer**/
                    legendItem type=line name="spacer" / label=' ' lineattrs=(thickness=0pt color=&background);
                    /*Creates Censor Legend*/
                    legendItem type=marker name="cens&z" / markerattrs=(color=&fontcolor symbol=plus) label='Censor'
                        labelattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                    %do i = 1 %to %superq(nclass_&z);
                        legendItem type=line name="plot&i._&z" /  
                            lineattrs=(thickness=&&linesize&z
                                %if %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                    color=%superq(color&z)
                                    %if %qupcase(%superq(pattern&z)) = AUTO %then %do;
                                        /**When all colors are the same, AUTO makes each pattern different**/
                                        pattern=&i
                                    %end;
                                    %else %do;
                                        %if %sysfunc(countw(%superq(pattern&z))) > 1 %then %do; pattern=%scan(%superq(pattern&z),&i) %end;
                                        %else %do; pattern=%superq(pattern&z) %end;
                                    %end;
                                %end;
                                %else  %do;
                                    %if %upcase(%superq(pattern&z)) = AUTO %then %do;
                                        /**When all colors are the different, AUTO makes each pattern solid**/
                                        pattern=solid
                                    %end;
                                    %else %do;
                                        %if %sysfunc(countw(%superq(pattern&z))) > 1 %then %do; pattern=%scan(%superq(pattern&z),&i) %end;
                                        %else %do; pattern=%superq(pattern&z) %end;
                                    %end;                            
                                    color=%scan(%superq(color&z), &i)                                        
                                %end;) 
                            label=' ' labelattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                        
                        %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do;
                            legendItem type=line name="plot&i.mv_&z" /  
                                lineattrs=(thickness=&&linesize&z
                                        %if %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                            color=%superq(color&z)
                                            %if %qupcase(%superq(pattern_adjust&z)) = AUTO %then %do;
                                                /**When all colors are the same, AUTO makes each pattern different**/
                                                %if &&plot_unadjust&z=1 %then %do;
                                                    pattern=%sysevalf(&&nclass_&z + &i)
                                                %end;
                                                %else %do;
                                                    pattern=&i
                                                %end;
                                            %end;
                                            %else %do;
                                                %if %sysfunc(countw(%superq(pattern_adjust&z))) > 1 %then %do; pattern=%scan(%superq(pattern_adjust&z),&i) %end;
                                                %else %do; pattern=%superq(pattern_adjust&z) %end;
                                            %end;
                                        %end;
                                        %else  %do;
                                            %if %upcase(%superq(pattern_adjust&z)) = AUTO %then %do;
                                                /**When all colors are the different, AUTO makes each pattern solid**/
                                                %if &&plot_unadjust&z=1 %then %do; pattern=2 %end;
                                                %else %do; pattern=solid %end;
                                            %end;
                                            %else %do;
                                                %if %sysfunc(countw(%superq(pattern_adjust&z))) > 1 %then %do; pattern=%scan(%superq(pattern_adjust&z),&i) %end;
                                                %else %do; pattern=%superq(pattern_adjust&z) %end;
                                            %end;                            
                                            color=%scan(%superq(color&z), &i)                                        
                                        %end;)
                                label=' ' labelattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                        %end;
                    %end;
                %end;
                /**Creates outer lattice block to contain all model plots**/
                layout lattice / columns=&columns rows=&rows opaque=false columngutter=&columngutter rowgutter=&rowgutter
                    order=&order columndatarange=DATA rowdatarange=DATA opaque=false;
                    /**Begins to fill in each cell of outer lattice block**/
                    %do z = 1 %to &nmodels;
                        layout lattice / columns=1 rows=1 opaque=false;
                            /**Creates footnotes at bottom of inner lattice block**/
                            %if %sysevalf(%superq(footnote&z)=,boolean)=0 %then %do;
                                /**SIDEBAR block extends entire bottom length of inner lattice block**/
                                sidebar / align=bottom;
                                    /**Layout Gridded allows multiple rows of ENTRY statements within one SIDEBAR**/
                                    layout gridded / rows=%sysfunc(countw(%superq(footnote&z),`,m)) border=false;
                                        %do i = 1 %to %sysfunc(countw(%superq(footnote&z),`,m));
                                            entry halign=%superq(footnotealign&z) "%scan(%superq(footnote&z),&i,`,m)" /
                                                textattrs=(color=&fontcolor weight=%superq(fnweight&z) size=%superq(fnsize&z) family="%superq(fnfamily&z)");
                                        %end;
                                    endlayout;
                                endsidebar;
                            %end;
                            
                            /**Creates a SIDEBAR block in the plot layout block to allow for individual model titles**/
                            sidebar /align=top;
                                %if %sysevalf(%superq(title&z)=,boolean)=0 %then %do;
                                    layout gridded / rows=%sysfunc(countw(%superq(title&z),`,m)) border=false;
                                        %do i = 1 %to %sysfunc(countw(%superq(title&z),`,m));
                                            entry halign=%superq(titlealign&z) "%scan(%superq(title&z),&i,`,m)" / 
                                                textattrs=(color=&fontcolor weight=%superq(tweight&z) size=%superq(tsize&z) family="%superq(tfamily&z)");
                                        %end;
                                    endlayout;
                                %end;
                            endsidebar;
                            /**Creates inner lattice block.  Adds a row if RISKLOCATION=BOTTOM.  Row Weights determined by RISKROWWEIGHTS**/
                            layout lattice / columns=1 columndatarange=union opaque=false 
                                %if &uniformheight=1 and &nmodels > 1 %then %do;
                                    rows=&_maxrows
                                    rowweights=(%sysevalf(1-%sysevalf((&_maxrows-1)*&_maxrowweights)) 
                                    %do i=2 %to %sysevalf(&_maxrows);
                                        &_maxrowweights
                                    %end;)
                                %end;
                                %else %do;
                                    rows=%superq(_rows&z)
                                    rowweights=(%sysevalf(1-(%superq(_rows&z)-1)*%superq(riskrowweights&z))
                                                %do i = 2 %to %superq(_rows&z);
                                                    %superq(riskrowweights&z)
                                                %end;)
                                %end;
                                rowgutter=0;
                                
                                rowheaders;
                                   layout gridded / columns=%sysfunc(countw(%superq(ylabel&z),`)) rows=1;
                                        %do k = 1 %to %sysfunc(countw(%superq(ylabel&z),`));
                                            entry halign=right "%scan(%superq(ylabel&z),&k,`)" / rotate=90 valign=center
                                                textattrs=(color=&fontcolor size=%superq(lsize&z) weight=%superq(lweight&z) family="%superq(lfamily&z)");
                                        %end;
                                    endlayout;
                                    %if %sysevalf(%superq(risklist&z)^=,boolean) and %sysevalf(%qupcase(%superq(risklocation&z))=BOTTOM,boolean) %then 
                                        %do k = 1 %to %sysfunc(countw(%superq(pardisplay&z),%str( )));
                                        %if (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) and %sysevalf(%superq(parheader&z)^=,boolean)) or
                                             (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) and %sysevalf(%superq(ncensheader&z)^=,boolean)) or
                                             (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) and %sysevalf(%superq(neventsheader&z)^=,boolean)) or
                                             (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) and 
                                                (%sysevalf(%superq(parheader&z)^=,boolean) and %sysevalf(%superq(ncensheader&z)^=,boolean))) or
                                             (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) and 
                                                (%sysevalf(%superq(parheader&z)^=,boolean) and %sysevalf(%superq(neventsheader&z)^=,boolean))) %then %do; 
                                            %if %sysevalf(%qupcase(%superq(paralign&z))=LABELS,boolean) %then %do;
                                                layout gridded;
                                                    entry ' ' / textattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                    drawtext  textattrs=(color=&fontcolor weight=%superq(parweight&z) size=%superq(parsize&z) family="%superq(parfamily&z)")
                                                        %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) %then %do;
                                                            "%superq(parheader&z)"
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) %then %do;
                                                            "%superq(ncensheader&z)"
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) %then %do;
                                                            "%superq(neventsheader&z)"
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) %then %do;
                                                            "%superq(parheader&z) (%superq(ncensheader&z))"
                                                        %end; 
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) %then %do;
                                                            "%superq(parheader&z) (%superq(neventsheader&z))"
                                                        %end; / 
                                                        drawspace=layoutpercent y=50 x=0 anchor=left justify=left width=10000;
                                                endlayout;
                                            %end;
                                            %else %do;
                                                entry ' ' / textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                            %end;
                                        %end;
                                        %local adj;
                                        %do adj=(&&plot_unadjust&z=0) %to %sysevalf(%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT);
                                            %if &adj=1 and &&plot_unadjust&z=1 %then %let sfx=%str( (Adjusted));
                                            %else %let sfx=;
                                            %let b=1;%let c_tally=0;
                                            %do i = 1 %to %superq(nclass_&z);
                                                %let c_tally=%sysevalf(&c_tally + 1);
                                                %if %sysevalf(%superq(by&z)^=,boolean) and ((&b=1 and &i=1) or (&c_tally>%superq(nclass_&b._&z.))) %then %do;
                                                    %if &c_tally>%superq(nclass_&b._&z.) %then %do;
                                                        %let b=%sysevalf(&b+1);
                                                        %let c_tally=0;
                                                    %end;
                                                    layout gridded / rows=1 columns=1 border=false opaque=false ;
                                                        /**Draws line above risk list headers when set to bottom**/
                                                        %if &&draw_underlines&z=1 %then %do;
                                                            drawline x1=0 x2=100 y1=100 y2=100 / drawspace=layoutpercent lineattrs=(color=&&ul_color&z pattern=&&ul_pattern&z thickness=&&ul_size&z);
                                                        %end;
                                                        entry " "/ valign=top 
                                                            textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=white);
                                                        drawtext textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=black weight=bold) "%superq(by_&z._&b)"/ 
                                                            drawspace=layoutpercent justify=left y=50 x=0 width=10000 widthunit=percent anchor=left;
                                                    endlayout;
                                                %end;
                                                %if %sysevalf(%qupcase(%superq(risklabellocation&z))^=LEFT,boolean) %then %do;
                                                    %if %sysevalf(%qupcase(%superq(risklabellocation&z))=ABOVE,boolean) %then %do;
                                                        entry ' ' / textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                    %end;
                                                    entry ' ' / textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                %end;
                                                %else %do;
                                                    layout gridded / rows=1 columns=1 border=false opaque=false ;
                                                        /**Draws line above risk list headers when set to bottom**/
                                                        %if &&draw_underlines&z=1 and %sysevalf(%superq(by&z)=,boolean) and &i=1 %then %do;
                                                            drawline x1=0 x2=100 y1=105 y2=105 / drawspace=layoutpercent lineattrs=(color=&&ul_color&z pattern=&&ul_pattern&z thickness=&&ul_size&z);
                                                        %end;
                                                        entry halign=right %unquote("%superq(class_&z._&i)&sfx.%superq(risklabeldlm&z)")
                                                        / textattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                            %if %superq(riskcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                color=%superq(color&z)
                                                            %end;
                                                            %else %if %superq(riskcolor&z)=1 %then %do;
                                                                color=%scan(%superq(color&z), &i)
                                                            %end;) valign=center;
                                                    endlayout;
                                                %end;
                                            %end;
                                        %end;
                                    %end;
                                endrowheaders;
                                /**Sets up axes for the Kaplan-Meier Curves**/
                                layout overlay /  
                                    %if &showwalls=1 %then %do;
                                        walldisplay=(outline)
                                    %end;
                                    %else %do;
                                        walldisplay=none
                                    %end;
                                    /**Y-Axis**/
                                    yaxisopts=(                
                                        display=(line ticks tickvalues )     
                                        %if &&gridlines&z=1 and (%sysevalf(%qupcase(&&gridline_axis&z)=Y,boolean) or %sysevalf(%qupcase(&&gridline_axis&z)=BOTH,boolean)) %then %do;
                                            griddisplay=on gridattrs=(thickness=&&gridline_size&z color=&&gridline_color&z pattern=&&gridline_pattern&z)
                                        %end;
                                        LABELSPLITCHAR='`' LABELFITPOLICY=SPLITALWAYS
                                        label="&&ylabel&z" labelattrs=(color=&fontcolor size=%superq(lsize&z) weight=%superq(lweight&z) family="%superq(lfamily&z)")
                                        type=linear tickvalueattrs=(color=&fontcolor size=%superq(ytickvalsize&z) weight=%superq(ytickvalweight&z) family="%superq(ytickvalfamily&z)")
                                        /**Offset creates space at the top or bottom of the window that the plot cannot use, proportion from 0 to 1**/
                                        %if %sysevalf(%superq(ymaxoffset&z)^=,boolean) %then %do; offsetmax=%superq(ymaxoffset&z)%end;
                                        %if %sysevalf(%superq(yminoffset&z)^=,boolean) %then %do; offsetmin=%superq(yminoffset&z) %end;   
                                        linearopts=(tickvaluesequence=(start=%superq(ymin&z) end=%superq(ymax&z) increment=%superq(yincrement&z))
                                        /**VIEWMAX and VIEWMIN are also required to show the desired range**/
                                        viewmin=%superq(ymin&z) viewmax=%superq(ymax&z)))
                                    /**X-Axis**/
                                    xaxisopts=(display=(line ticks tickvalues label)  label="&&xlabel&z"
                                        LABELSPLITCHAR='`' LABELFITPOLICY=SPLITALWAYS
                                        %if &&gridlines&z=1 and (%sysevalf(%qupcase(&&gridline_axis&z)=X,boolean) or %sysevalf(%qupcase(&&gridline_axis&z)=BOTH,boolean)) %then %do;
                                            griddisplay=on gridattrs=(thickness=&&gridline_size&z color=&&gridline_color&z pattern=&&gridline_pattern&z)
                                        %end;
                                        type=linear labelattrs=(color=&fontcolor size=%superq(lsize&z) weight=%superq(lweight&z) family="%superq(lfamily&z)")
                                        tickvalueattrs=(color=&fontcolor size=%superq(xtickvalsize&z) weight=%superq(xtickvalweight&z) family="%superq(xtickvalfamily&z)")               
                                        /**Offset creates space at the top or bottom of the window that the plot cannot use, proportion from 0 to 1**/
                                        %if %sysevalf(%superq(xmaxoffset&z)^=,boolean) %then %do; offsetmax=%superq(xmaxoffset&z)%end;
                                        %if %sysevalf(%superq(xminoffset&z)^=,boolean) %then %do; offsetmin=%superq(xminoffset&z)%end;
                                        %else %if %sysevalf(%superq(risklist&z)=,boolean) or %sysevalf(%qupcase(&&risklocation&z)=BOTTOM,boolean) %then %do; offsetmin=0 %end;
                                        /**TICKVALUESEQUENCE automatically calculates tick marks**/
                                        linearopts=(tickvaluesequence=(start=%superq(xmin&z) end=%superq(xmax&z) increment=%superq(xincrement&z))
                                        /**VIEWMAX and VIEWMIN are also required to show the desired range**/
                                        viewmin=%superq(xmin&z) viewmax=%superq(xmax&z)));
                                
                                    /**Confidence Bounds if flagged**/
                                    %if %superq(plotci&z) = 1 or 
                                        (%superq(plotci&z)=2 and %sysevalf(%superq(class&z)=,boolean)) %then %do;
                                        
                                        %if &&plot_unadjust&z=1 %then %do;
                                            bandplot x=t_&z limitlower=eval(%superq(xmult_&z) *lcl_&z) 
                                                limitupper=eval(%superq(xmult_&z) *ucl_&z) / type=step group=_ua_ci_&z
                                                display=(%if %sysevalf(%superq(plotcifill&z)=1,boolean) %then %do;
                                                            fill
                                                        %end;
                                                        %else %do;
                                                            outline 
                                                        %end;);
                                        %end;
                                        %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do;
                                            bandplot x=t_adj_&z limitlower=eval(%superq(xmult_&z) *lcl_adj_&z) 
                                                limitupper=eval(%superq(xmult_&z) *ucl_adj_&z) / type=step group=_a_ci_&z
                                                display=(%if %sysevalf(%superq(plotcifill&z)=1,boolean) %then %do;
                                                            fill
                                                        %end;
                                                        %else %do;
                                                            outline 
                                                        %end;);
                                        %end;
                                    %end;
                                    
                                    %if %superq(censormarkers&z) ^= 0 %then %do;           
                                        /**Draws censor symbols*/
                                        /**xmult_&z is a factor based on YTYPE (PPT vs PCT)**/
                                        %if &&plot_unadjust&z=1 %then %do;
                                            if (max(c_&z)>.)
                                                scatterplot x=t_&z y=eval(%superq(xmult_&z) * c_&z) / 
                                                    group=_ua_&z;
                                            endif;
                                        %end;
                                        %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do;
                                            if (max(c_&z)>.)
                                                scatterplot x=t_adj_&z y=eval(%superq(xmult_&z) * c_adj_&z) / 
                                                    group=_a_&z;
                                            endif;                                    
                                        %end;
                                    %end;
                                    
                                    
                                    /**Generates the Kaplan-Meier Curves**/
                                    /**Legend labels are added manually with ENTRY statements later**/ 
                                    %if &&plot_unadjust&z=1 %then %do;
                                        stepplot x=t_&z y=eval(%superq(xmult_&z) * s_&z) / legendlabel=' ' group=_ua_&z;
                                    %end;
                                    %if %qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT %then %do;
                                        stepplot x=t_adj_&z y=eval(%superq(xmult_&z) * s_adj_&z) / legendlabel=' ' group=_a_&z;
                                    %end;;
    
                                    /**Draw X axis Reference Lines**/
                                    %if (%sysevalf(%qupcase(%superq(reflines&z))=TIMEPOINTS,boolean) and %sysevalf(%superq(timelist&z)^=,boolean)) or 
                                        %sysevalf(%qupcase(%superq(reflines&z))=MEDIANS,boolean) %then %do;  
                                        /**Full lines**/
                                        %if %sysevalf(%qupcase(%superq(reflinemethod&z))=FULL,boolean) %then %do;
                                            %if %sysevalf(%qupcase(%superq(reflineaxis&z))=X,boolean) or 
                                                %sysevalf(%qupcase(%superq(reflineaxis&z))=BOTH,boolean) %then %do;
                                                referenceline x=ref_t_&z / 
                                                    lineattrs=(thickness=%superq(reflinesize&z) color=%superq(reflinecolor&z) pattern=%superq(reflinepattern&z));
                                                %end;                       
                                            %if %sysevalf(%qupcase(%superq(reflineaxis&z))=Y,boolean) or 
                                                %sysevalf(%qupcase(%superq(reflineaxis&z))=BOTH,boolean) %then %do;
                                                referenceline y=ref_y_&z / 
                                                    lineattrs=(thickness=%superq(reflinesize&z) color=%superq(reflinecolor&z) pattern=%superq(reflinepattern&z));
                                            %end;
                                        %end; 
                                        /**Drop lines**/ 
                                        %else %if %sysevalf(%qupcase(%superq(reflinemethod&z))=DROP,boolean) %then %do; 
                                            %if %sysevalf(%qupcase(%superq(reflineaxis&z))=X,boolean) or 
                                                %sysevalf(%qupcase(%superq(reflineaxis&z))=BOTH,boolean) %then %do;
                                                dropline x=ref_t_&z y=ref_y_&z / 
                                                    lineattrs=(thickness=%superq(reflinesize&z) color=%superq(reflinecolor&z) pattern=%superq(reflinepattern&z))
                                                    dropto=X;    
                                            %end; 
                                            %if %sysevalf(%qupcase(%superq(reflineaxis&z))=Y,boolean) or 
                                                %sysevalf(%qupcase(%superq(reflineaxis&z))=BOTH,boolean) %then %do;
                                                dropline x=ref_t_&z y=ref_y_&z / 
                                                    lineattrs=(thickness=%superq(reflinesize&z) color=%superq(reflinecolor&z) pattern=%superq(reflinepattern&z))
                                                    dropto=Y;    
                                            %end;
                                        %end;
                                    %end;/**End Reference Lines**/
                                    /**Patients-at-Risk table when RISKLOCATION=INSIDE**/
                                    %if %sysevalf(%superq(risklist&z)=,boolean)=0 and %qupcase(%superq(risklocation&z))=INSIDE %then %do k=1 %to %sysfunc(countw(%superq(pardisplay&z),%str( ))) %by 1;
                                        innermargin / align=bottom;
                                            /**Makes a header for the Patients-at-Risk table**/
                                             %if (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) and %sysevalf(%superq(parheader&z)^=,boolean)) or
                                                 (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) and %sysevalf(%superq(ncensheader&z)^=,boolean)) or
                                                 (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) and %sysevalf(%superq(neventsheader&z)^=,boolean)) or
                                                 (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) and 
                                                    (%sysevalf(%superq(parheader&z)^=,boolean) and %sysevalf(%superq(ncensheader&z)^=,boolean))) or
                                                 (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) and 
                                                    (%sysevalf(%superq(parheader&z)^=,boolean) and %sysevalf(%superq(neventsheader&z)^=,boolean))) %then %do;
                                                blockplot x=eval(ifn(^missing(partitle_&z),%superq(xmin&z),.)) 
                                                    %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) %then %do;
                                                        block=partitle_&z 
                                                    %end;
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) %then %do;
                                                        block=ncenstitle_&z 
                                                    %end;
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) %then %do;
                                                        block=neventstitle_&z 
                                                    %end;
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) %then %do;
                                                        block=eval(strip(partitle_&z)||' ('||strip(ncenstitle_&z)||')')
                                                    %end; 
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) %then %do;
                                                        block=eval(strip(partitle_&z)||' ('||strip(neventstitle_&z)||')')
                                                    %end;
                                                    /
                                                    %if %sysevalf(%qupcase(%superq(paralign&z))=LABELS,boolean) %then %do;
                                                        display=(label) 
                                                        %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) %then %do;
                                                            label="%superq(parheader&z)"
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) %then %do;
                                                            label="%superq(ncensheader&z)"
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) %then %do;
                                                            label="%superq(neventsheader&z)"
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) %then %do;
                                                            label="%superq(parheader&z) (%superq(ncensheader&z))"
                                                        %end; 
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) %then %do;
                                                            label="%superq(parheader&z) (%superq(neventsheader&z))"
                                                        %end;
                                                        labelattrs=(color=&fontcolor size=%superq(parsize&z) weight=%superq(parweight&z) family="%superq(parfamily&z)")
                                                    %end;
                                                    %else %do;
                                                        display=(values)
                                                        valuehalign=%superq(paralign&z) 
                                                    %end;
                                                    valueattrs=(color=&fontcolor size=%superq(parsize&z) weight=%superq(parweight&z) family="%superq(parfamily&z)");
                                            %end;
                                            /**Makes one block plot per class level**/
                                            %local adj sfx sfx2;
                                            %do adj=(&&plot_unadjust&z=0) %to %sysevalf(%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT);
                                                %if &adj=1 %then %do;
                                                    %let sfx=_adj;
                                                    %if &&plot_unadjust&z=1 %then %let sfx2=%str((Adjusted));
                                                    %else %let sfx2=;
                                                %end;
                                                %else %do;
                                                    %let sfx=;
                                                    %let sfx2=;
                                                %end;
                                                %let b=1;%let c_tally=0;
                                                %do i = 1 %to %superq(nclass_&z) %by 1;
                                                    %let c_tally=%sysevalf(&c_tally + 1);
                                                    %if %sysevalf(%superq(by&z)^=,boolean) and ((&b=1 and &i=1) or (&c_tally>%superq(nclass_&b._&z.))) %then %do;
                                                        %if &c_tally>%superq(nclass_&b._&z.) %then %do;
                                                            %let b=%sysevalf(&b+1);
                                                            %let c_tally=1;
                                                        %end;
                                                        blockplot x=eval(ifn(by_order&z=&b and cl_&z="%superq(class_&z._&i)",%superq(xmin&z),.)) 
                                                            block=by_level&z / 
                                                            %if %qupcase(%superq(risklabellocation&z))^=LEFT %then %do;
                                                                display=(values) 
                                                            %end;
                                                            %else %do;
                                                                display=(label) 
                                                            %end;
                                                            label="%superq(by_&z._&b)"
                                                            valuehalign=%superq(risklabelalign&z)
                                                            valueattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" weight=bold)
                                                            labelattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" weight=bold);
                                                    %end;
                                                    %if %qupcase(%superq(risklabellocation&z))=ABOVE and
                                                        (%sysevalf(%superq(class&z)^=,boolean) or %sysevalf(%superq(classdesc&z)^=,boolean)) %then %do;
                                                        blockplot x=eval(ifn(t&sfx._&z>. and cl_&z="%superq(class_&z._&i)",%superq(xmin&z),.)) 
                                                            block=eval(catx(' ',cl_&z,"&sfx2.")) / display=(values)
                                                            valuehalign=%superq(risklabelalign&z)
                                                            valueattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" weight=%superq(risklabelweight));
                                                    %end; 
                                                    blockplot x=eval(ifn(time&i.&sfx._&z ge %superq(xmin&z),time&i.&sfx._&z,.)) 
                                                        %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) %then %do;
                                                            block=atrisk&i.&sfx._&z 
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) %then %do;
                                                            block=ncens&i.&sfx._&z 
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) %then %do;
                                                            block=nevent&i.&sfx._&z 
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) %then %do;
                                                            block=eval(strip(atrisk&i.&sfx._&z)||' ('||strip(ncens&i.&sfx._&z)||')')
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) %then %do;
                                                            block=eval(strip(atrisk&i.&sfx._&z)||' ('||strip(ncens&i.&sfx._&z)||')')
                                                        %end;
                                                        /repeatedvalues=TRUE
                                                        %if %qupcase(%superq(risklabellocation&z))=LEFT %then %do;
                                                            display=(values label)                
                                                            label=
                                                                %if %sysevalf(%superq(class_&z._&i)^=,boolean) %then %do;
                                                                    "%superq(class_&z._&i) &sfx2."
                                                                %end;
                                                                %else %do;
                                                                    " "
                                                                %end;
                                                            labelattrs=(color=&fontcolor size=%superq(ptabsize&z) weight=%superq(risklabelweight&z) family="%superq(ptabfamily&z)")
                                                        %end;
                                                        %else %do;
                                                            display=(values)
                                                        %end;
                                                        valuehalign=start
                                                        valueattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                            %if %superq(riskcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                color=%superq(color&z)
                                                            %end;
                                                            %else %if %superq(riskcolor&z)=1 %then %do;
                                                                color=%scan(%superq(color&z), &i)
                                                            %end;
                                                            %else %do;
                                                                color=&fontcolor 
                                                            %end;);  
                                                %end;
                                            %end;
                                        endinnermargin;                         
                                        /**Places a reference line at the minimum Y-axis value to separate the patients-at-risk table from the plot**/               
                                        %if %superq(riskdivider&z)=1 %then %do;
                                            referenceline y=%superq(ymin&z) / lineattrs=(color=%superq(riskdivcolor) pattern=%superq(riskdivstyle&z));
                                        %end;
                                    %end;/**Ends Inner 
                                    
                                    /**Design the Statistics Summary Table**/
                                    %if %superq(_ndisplay_class_&z) gt 0 or %superq(_ndisplay_model_&z)>0 or %superq(_ndisplay_interp_&z)>0 %then %do; 
                                        /**Creates the outer gridded block with up to 3 rows**/ 
                                        layout gridded / 
                                            rows=%sysevalf(%sysevalf(%superq(_ndisplay_class_&z)>0,boolean) + %sysevalf(%superq(_ndisplay_mstats_&z)>0,boolean)
                                                             + %sysevalf(%superq(_ndisplay_istats_&z)>0,boolean)) columns=1 border=false
                                            location=%superq(location&z) autoalign=(%superq(autoalign&z)) %if &&gridlines&z=1 and &transparent=0 %then %do; opaque=true %end; ;
                                            /**Creates the class-level gridded block**/
                                            %if %sysevalf(%superq(_ndisplay_class_&z)>0,boolean) %then %do;
                                                layout gridded / 
                                                    %if %qupcase(&&method&z)=KM or %qupcase(&&method&z)=CIF %then %do;
                                                        rows=%sysevalf(%superq(nclass_&z)+1+%sysevalf(%superq(by&z)^=,boolean)*&&nby_&z +
                                                                        %sysevalf(%superq(_ndisplay_mstats_&z)>0,boolean)*(&&nby_&z-1)) /**One row per class level**/                                                
                                                    %end;
                                                    %else %do;
                                                        rows=%sysevalf((1+&&plot_unadjust&z)*%superq(nclass_&z)+1+%sysevalf(%superq(by&z)^=,boolean)*&&nby_&z +
                                                                        %sysevalf(%superq(_ndisplay_mstats_&z)>0,boolean)*(&&nby_&z-1)) /**One row per class level**/ 
                                                    %end;
                                                    columns=%superq(_ndisplay_class_&z) /**Determined by number of statistics called in the DISPLAY parameter**/                        
                                                    border=false valign=top halign=center;
                                                    /**Creates the headers for the class-level gridded block**/
                                                    %local pclass_check;
                                                    %do i = 1 %to %superq(_ndisplay_class_&z);
                                                        %macro _underline_(y=5); 
                                                            /**Draws lines above below headers**/
                                                            %if &&draw_underlines&z=1 and %sysevalf(%superq(by&z)=,boolean) %then %do;
                                                                drawline x1=0 x2=100 y1=&y y2=&y / drawspace=layoutpercent lineattrs=(color=&&ul_color&z pattern=&&ul_pattern&z thickness=&&ul_size&z);
                                                            %end;
                                                        %mend;
                                                        %let pclass_check=0;
                                                        %if %superq(_display_class&i._&z)=LEGEND %then %do;
                                                            %if %sysevalf(%superq(by&z)^=,boolean) %then %let legendheader&z=%superq(bylabel&z);
                                                            /**Legend**/
                                                            layout gridded / columns=1 rows=%sysfunc(countw(%superq(legendheader&z),`,m)) border=false
                                                                halign=center valign=bottom;
                                                                %do k=1 %to %sysfunc(countw(%superq(legendheader&z),`,m));
                                                                    entry halign=left "%scan(%superq(legendheader&z),&k,`,m)" / 
                                                                        valign=bottom textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                %end;
                                                                %_underline_;
                                                            endlayout;
                                                            /**Class Levels**/
                                                            %if %sysevalf("%superq(classdesc&z)"^="",boolean) and %sysevalf(%superq(class&z)^=,boolean) %then %do;
                                                                layout gridded / columns=1 rows=%sysfunc(countw(%superq(classdesc&z),`,m)) border=false
                                                                    halign=center valign=bottom;
                                                                    %do k=1 %to %sysfunc(countw(%superq(classdesc&z),`,m));
                                                                        entry halign=center "%scan(%superq(classdesc&z),&k,`,m)" / valign=bottom 
                                                                            textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                    %end;
                                                                    %_underline_;
                                                                endlayout;
                                                            %end;
                                                            %else %if %sysevalf(%superq(class&z)^=,boolean) %then %do;
                                                                layout gridded / columns=1 rows=%sysfunc(countw(%superq(label&z),`,m)) border=false
                                                                    halign=center valign=bottom;
                                                                    %do k=1 %to %sysfunc(countw(%superq(label&z),`,m));
                                                                        entry halign=center "%scan(%superq(label&z),&k,`,m)" / valign=bottom 
                                                                            textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                    %end;
                                                                    %_underline_;                                                                    
                                                                endlayout;
                                                            %end;
                                                            %else %do;
                                                                layout gridded / columns=1 rows=1 border=false opaque=false halign=center valign=bottom;
                                                                    entry halign=center " " / valign=bottom 
                                                                        textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                    %_underline_;
                                                                endlayout;
                                                            %end;
                                                            %let i = %sysevalf(&i+1);
                                                            %let pclass_check=1;  
                                                        %end;/**End Legend**/
                                                        %if %superq(pclass_check)=0 %then %do j=1 %to %sysfunc(countw(&classcolumns,|));
                                                            /**Total|Events|Medians**/
                                                            %if %superq(_display_class&i._&z)=%scan(&classcolumns,&j,|) %then %do;
                                                                %local header;
                                                                %let header=%superq(%sysfunc(compress(%superq(_display_class&i._&z)header))&z);
                                                                layout gridded / columns=1 rows=%sysfunc(countw(%superq(header),`,m)) border=false
                                                                    halign=center valign=bottom;
                                                                    %do k=1 %to %sysfunc(countw(%superq(header),`,m));
                                                                        entry halign=center "%scan(&header,&k,`,m)" / valign=bottom
                                                                            textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                    %end;
                                                                    %_underline_;
                                                                endlayout;
                                                                %let pclass_check=1;  
                                                            %end;
                                                        %end;/**Ends Total/Events/Medians**/
                                                        %if %superq(pclass_check)=0 and %sysevalf(%superq(class&z)^=,boolean) %then %do j=1 %to %sysfunc(countw(&hrcolumns,|));
                                                            /**Hazard Ratios**/
                                                            %if %superq(_display_class&i._&z)=%scan(&hrcolumns,&j,|) %then %do;
                                                                %local header;
                                                                %let header=%superq(%sysfunc(compress(%superq(_display_class&i._&z)header))&z);
                                                                layout gridded / columns=1 rows=%sysfunc(countw(%superq(header),`,m)) border=false
                                                                    halign=center valign=bottom;
                                                                    %do k=1 %to %sysfunc(countw(%superq(header),`,m));
                                                                    entry halign=center "%scan(&header,&k,`,m)" / valign=bottom 
                                                                        textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                    %end;
                                                                    %_underline_;
                                                                endlayout;
                                                                %let pclass_check=1;  
                                                            %end;
                                                        %end;/**Ends Hazard Ratios**/
                                                        %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=TIMELIST %then %do;
                                                            /**Survival Time-point Estimates**/
                                                            %if %sysevalf(%superq(listtimepoints&z)=1,boolean) %then %do;
                                                                layout gridded / columns=1 rows=%sysfunc(countw(%superq(timelistheader&z),`,m)) border=false
                                                                    halign=center valign=bottom;
                                                                    %do k=1 %to %sysfunc(countw(%superq(timelistheader&z),`,m));
                                                                        entry halign=center "%scan(&&timelistheader&z,&k,`,m)" / valign=bottom 
                                                                            textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                    %end;
                                                                    %_underline_;
                                                                 endlayout;
                                                            %end;
                                                            layout gridded / columns=1 rows=%sysfunc(countw(%superq(kmestheader&z),`,m)) border=false
                                                                halign=center valign=bottom;
                                                                %do k=1 %to %sysfunc(countw(%superq(kmestheader&z),`,m));
                                                                    %if %sysevalf(%superq(kmestheader&z)=%str(KM Est %(95%% CI%)),boolean) and
                                                                        %sysevalf(%qupcase(%superq(method&z))=CIF,boolean) %then %let kmestheader&z=CIF Est (95% CI);
                                                                    %else %if %sysevalf(%superq(kmestheader&z)=%str(KM Est %(95%% CI%)),boolean) and
                                                                        %sysevalf(%qupcase(%superq(sreverse&z))=1,boolean) %then %let kmestheader&z=1-KM Est (95% CI);
                                                                    entry halign=center "%scan(&&kmestheader&z,&k,`,m)" / valign=bottom 
                                                                        textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                %end;
                                                                %_underline_;
                                                            endlayout;
                                                            %let i = %sysevalf(&i+%superq(listtimepoints&z));
                                                            %let pclass_check=1; 
                                                            %if  %qupcase(%superq(risklocation&z))= TIMELIST and %sysevalf(%superq(risklist&z)^=,boolean) %then %do;
                                                                /**Patients-at-Risk Table**/    
                                                                layout gridded / columns=1 rows=%sysfunc(countw(%superq(risktableheader&z),`,m)) border=false
                                                                    halign=center valign=bottom;
                                                                    %do k=1 %to %sysfunc(countw(%superq(risktableheader&z),`,m));
                                                                        entry halign=center "%scan(&&risktableheader&z,&k,`,m)" / valign=bottom 
                                                                            textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                    %end;
                                                                    %_underline_;
                                                                endlayout;                                
                                                                %let i = %sysevalf(&i+1);
                                                            %end;
                                                        %end; 
                                                        

                                                        %local text;
                                                        %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=PVAL %then %do;
                                                            /**P-value**/
                                                            %if %qupcase(%superq(plotpval&z))=LR %then %let plotpval&z=Likelihood-Ratio; 
                                                            %else %if %qupcase(%superq(plotpval&z))=GRAY %then %let plotpval&z=Gray K-Sample Test; 
                                                            %else %if %qupcase(%superq(plotpval&z))=LOGRANK_ONESIDED %then %let plotpval&z=One-sided logrank; 
                                                            %if %qupcase(%superq(plotpvalmv&z))=LR %then %let plotpvalmv&z=Likelihood-Ratio; 
                                                            
                                                            %if %sysevalf(%superq(pvalheader&z)=,boolean)=0 %then %let text=%unquote(&&pvalheader&z);
                                                            %else %if %sysevalf(%superq(strata&z)^=,boolean) %then %let text=Stratified`%sysfunc(propcase(%superq(plotpval&z))) P-value;
                                                            %else %let text=%sysfunc(propcase(%superq(plotpval&z))) P-value;
                                                            layout gridded / columns=1 rows=%sysfunc(countw(%superq(text),`,m)) border=false
                                                                halign=center valign=bottom;
                                                                %do k=1 %to %sysfunc(countw(%superq(text),`,m));
                                                                    entry halign=center "%scan(&text,&k,`,m)" / valign=bottom 
                                                                        textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                %end;
                                                                %_underline_;
                                                            endlayout;
                                                            %let pclass_check=1; 
                                                        %end;
                                                        %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=PVALMV %then %do;
                                                            /**P-value**/
                                                            %if %qupcase(%superq(plotpvalmv&z))=LR %then %let plotpvalmv&z=Likelihood-Ratio; 
                                                            
                                                            %if %sysevalf(%superq(pvalmvheader&z)=,boolean)=0 %then %let text=%unquote(&&pvalmvheader&z);
                                                            %else %if %sysevalf(%superq(strata&z)^=,boolean) %then %let text=Stratified Adjusted`%sysfunc(propcase(%superq(plotpvalmv&z))) P-value;
                                                            %else %let text=Adjusted`%sysfunc(propcase(%superq(plotpvalmv&z))) P-value;
                                                            
                                                            layout gridded / columns=1 rows=%sysfunc(countw(%superq(text),`,m)) border=false
                                                                halign=center valign=bottom;
                                                                %do k=1 %to %sysfunc(countw(%superq(text),`,m));
                                                                    entry halign=center "%scan(&text,&k,`,m)" / valign=bottom 
                                                                        textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                %end;
                                                                %_underline_;
                                                            endlayout;
                                                            %let pclass_check=1; 
                                                        %end;
                                                        %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=PVAL_INTER %then %do;
                                                            /**P-value**/
                                                            %if %qupcase(%superq(pval_inter&z))=LR %then %let pval_inter&z=Likelihood-Ratio; 
                                                            
                                                            %if %sysevalf(%superq(pval_interheader&z)=,boolean)=0 %then %let text=%unquote(&&pval_interheader&z);
                                                            %else %if %sysevalf(%superq(strata&z)^=,boolean) %then %let text=Stratified`%sysfunc(propcase(%superq(pval_inter&z)))`Interaction P-value;
                                                            %else %let text=%sysfunc(propcase(%superq(pval_inter&z)))`Interaction P-value;
                                                            layout gridded / columns=1 rows=%sysfunc(countw(%superq(text),`,m)) border=false
                                                                halign=center valign=bottom;
                                                                %do k=1 %to %sysfunc(countw(%superq(text),`,m));
                                                                    entry halign=center "%scan(&text,&k,`,m)" / valign=bottom 
                                                                        textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                %end;
                                                                %_underline_;
                                                            endlayout;
                                                            %let pclass_check=1; 
                                                        %end;
                                                        %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=PVAL_INTERMV %then %do;
                                                            /**P-value**/
                                                            %if %qupcase(%superq(pval_intermv&z))=LR %then %let pval_intermv&z=Likelihood-Ratio; 
                                                            
                                                            %if %sysevalf(%superq(pval_intermvheader&z)=,boolean)=0 %then %let text=%unquote(&&pval_intermvheader&z);
                                                            %else %if %sysevalf(%superq(strata&z)^=,boolean) %then %let text=Stratified`Adjusted %sysfunc(propcase(%superq(pval_intermv&z)))`Interaction P-value;
                                                            %else %let text=Adjusted %sysfunc(propcase(%superq(pval_intermv&z)))`Interaction P-value;
                                                            layout gridded / columns=1 rows=%sysfunc(countw(%superq(text),`,m)) border=false
                                                                halign=center valign=bottom;
                                                                %do k=1 %to %sysfunc(countw(%superq(text),`,m));
                                                                    entry halign=center "%scan(&text,&k,`,m)" / valign=bottom 
                                                                        textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                %end;
                                                                %_underline_;
                                                            endlayout;
                                                            %let pclass_check=1; 
                                                        %end;
                                                        %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=RMST_PVAL %then %do;                                                            
                                                            %if %sysevalf(&&rmst_pvalheader&z=,boolean)=0 %then %let text=&&rmst_pvalheader&z;
                                                            %else %if %sysevalf(%superq(strata&z)^=,boolean) %then %let text=Stratified`RMST Chi-square`P-value;
                                                            %else %let text=RMST Chi-square`P-value;
                                                            /**Restricted Means Survival Time P-value comparison**/
                                                            layout gridded / columns=1 rows=%sysfunc(countw(%superq(text),`,m)) border=false
                                                                halign=center valign=bottom;
                                                                %do k=1 %to %sysfunc(countw(%superq(text),`,m));
                                                                    entry halign=center "%scan(&text,&k,`,m)" / valign=bottom 
                                                                        textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                %end;
                                                                %_underline_;
                                                            endlayout;
                                                            %let pclass_check=1; 
                                                        %end;
                                                        %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=CINDEX %then %do;                                                       
                                                            %if %sysevalf(&&cindexheader&z=,boolean)=0 %then %let text=&&cindexheader&z;
                                                            %else %let text=C-index (%sysevalf(100-100*&&alpha&z)% CI);
                                                            /**Univariate c-index**/
                                                            layout gridded / columns=1 rows=%sysfunc(countw(%superq(text),`,m)) border=false
                                                                halign=center valign=bottom;
                                                                %do k=1 %to %sysfunc(countw(%superq(text),`,m));
                                                                    entry halign=center "%scan(&text,&k,`,m)" / valign=bottom 
                                                                        textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                %end;
                                                                %_underline_;
                                                            endlayout;
                                                            %let pclass_check=1; 
                                                        %end;
                                                        %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=CINDEXMV %then %do;                                          
                                                            %if %sysevalf(&&cindexmvheader&z=,boolean)=0 %then %let text=&&cindexmvheader&z;
                                                            %else %let text=Multivariate`C-index (%sysevalf(100-100*&&alpha&z)% CI);
                                                            /**Multivariate c-index**/
                                                            layout gridded / columns=1 rows=%sysfunc(countw(%superq(text),`,m)) border=false
                                                                halign=center valign=bottom;
                                                                %do k=1 %to %sysfunc(countw(%superq(text),`,m));
                                                                    entry halign=center "%scan(&text,&k,`,m)" / valign=bottom 
                                                                        textattrs=(color=&fontcolor weight=bold size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                %end;
                                                                %_underline_;
                                                            endlayout;
                                                            %let pclass_check=1; 
                                                        %end;
                                                    %end;
                                                    
                                                    /**Creates the values for the class-level gridded block**/
                                                    %local c adj sfx sfx2;
                                                    %if %sysevalf(%superq(ntl_&z)=,boolean) %then %let ntl_&z = 1;
                                                    %do adj = (&&plot_unadjust&z=0) %to (%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT);
                                                        %if &adj=1 %then %do;
                                                            %let sfx=mv;
                                                            %if &&plot_unadjust&z=1 %then %let sfx2=%str((Adjusted));
                                                            %else %let sfx2=;
                                                        %end;
                                                        %else %do;
                                                            %let sfx=;
                                                            %let sfx2=;
                                                        %end;
                                                        
                                                        %local c_tally;
                                                        %let b=1;%let c_tally=0;
                                                        %do c = 1 %to %superq(nclass_&z);  
                                                            %let c_tally=%sysevalf(&c_tally + 1);   
                                                            %if %sysevalf(%superq(by&z)^=,boolean) and ((&b=1 and &c=1) or (&c_tally>%superq(nclass_&b._&z.))) %then %do;
                                                                %if &c_tally>%superq(nclass_&b._&z.) %then %do;
                                                                    %local _mstat_notc;
                                                                    %if &&_ipval_&z=0 %then %let _mstat_notc=%sysevalf(%superq(_ndisplay_mstats_&z)
                                                                          -1*(%sysfunc(find(%superq(display&z),tablecomments,i))>0)*%sysfunc(countw(%superq(tablecomments&z),`,m)));
                                                                    %else %let _mstat_notc=&&_ndisplay_mstats_&z;

                                                                    %if %sysevalf(%superq(_ndisplay_mstats_&z)>0,boolean) and &c^=%superq(nclass_&z)  and 
                                                                        %sysevalf(&_mstat_notc>0,boolean) %then %do j = 1 %to %superq(_ndisplay_class_&z);
                                                                        %if &j^=2 %then %do;
                                                                            layout gridded / rows=1 columns=1 opaque=false border=false;
                                                                            entry " "/ valign=top 
                                                                                textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=white);
                                                                            endlayout;
                                                                        %end;
                                                                        %else %if &j=2 %then %do;
                                                                            /**Print text for P-values and user entered Table Comments**/                           
                                                                            layout gridded / columns=1 
                                                                                %if %qupcase(&&method&z)=KM or %qupcase(&&method&z)=CIF %then %do;
                                                                                    rows=%superq(_ndisplay_mstats_&z)                                         
                                                                                %end;
                                                                                %else %do;
                                                                                    rows=&_mstat_notc 
                                                                                %end;
                                                                                valign=top halign=left
                                                                                opaque=false border=false;
                                                                                %do k = 1 %to &_mstat_notc;
                                                                                    layout gridded / columns=1 rows=1 border=false opaque=false; 
                                                                                    entry " "/ valign=top 
                                                                                        textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=white);
                                                                                    %if %superq(_display_mstats&k._&z)=PVAL %then %do;
                                                                                        /**P-value**/
                                                                                        %if %qupcase(%superq(plotpval&z))=LR %then %let plotpval&z=Likelihood-Ratio; 
                                                                                        %else %if %qupcase(%superq(plotpval&z))=GRAY %then %let plotpval&z=Gray K-Sample Test; 
                                                                                        %else %if %qupcase(%superq(plotpval&z))=LOGRANK_ONESIDED %then %let plotpval&z=One-sided logrank; 
                                                                                        %if %qupcase(%superq(plotpvalmv&z))=LR %then %let plotpvalmv&z=Likelihood-Ratio; 
                                                                                        drawtext textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=black ) 
                                                                                            %if %sysevalf(%superq(pval&sfx.header&z)=,boolean)=0 %then %do;
                                                                                                "%unquote(&&pval&sfx.header&z) %sysfunc(strip(%qscan(&&pval&sfx.&z,&b,|,m)))"
                                                                                            %end;
                                                                                            %else %if %sysevalf(%superq(strata&z)^=,boolean) %then %do;
                                                                                                "Stratified %sysfunc(propcase(%superq(plotpval&sfx.&z))) P-value: %sysfunc(strip(%qscan(&&pval&sfx.&z,&b,|,m)))"
                                                                                            %end;
                                                                                            %else %do;
                                                                                                "&sfx2.%sysfunc(propcase(%superq(plotpval&sfx.&z))) P-value: %sysfunc(strip(%qscan(&&pval&sfx.&z,&b,|,m)))"
                                                                                            %end;/ 
                                                                                            drawspace=layoutpercent justify=left y=50 x=0 width=10000 widthunit=percent anchor=left;
                                                                                    %end;
                                                                                    %if %superq(_display_mstats&k._&z)=PVALMV %then %do;
                                                                                        /**P-value**/
                                                                                        %if %qupcase(%superq(plotpvalmv&z))=LR %then %let plotpvalmv&z=Likelihood-Ratio; 
                                                                                        drawtext textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=black ) 
                                                                                            %if %sysevalf(&&pvalmvheader&z=,boolean)=0 %then %do;
                                                                                            "&&pvalmvheader&z %sysfunc(strip(%qscan(%superq(pvalmv&z),&b,|,m)))"
                                                                                            %end;
                                                                                            %else %do;
                                                                                                "Adjusted %sysfunc(propcase(%superq(plotpvalmv&z))) P-value: %sysfunc(strip(%qscan(%superq(pvalmv&z),&b,|,m)))"
                                                                                            %end;/ 
                                                                                            drawspace=layoutpercent justify=left y=50 x=0 width=10000 widthunit=percent anchor=left;
                                                                                    %end;
                                                                                    %if %superq(_display_mstats&k._&z)=RMST_PVAL %then %do;
                                                                                        /**Restricted Means Survival Time P-value comparison**/
                                                                                        drawtext textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=black ) 
                                                                                            %if %sysevalf(&&rmst_pvalheader&z=,boolean)=0 %then %do;
                                                                                                "&&rmst_pvalheader&z %sysfunc(strip(%qscan(%superq(rmst_pval&z),&b,|,m)))"
                                                                                            %end;
                                                                                            %else %if %sysevalf(%superq(strata&z)^=,boolean) %then %do;
                                                                                                "Stratified RMST Chi-square P-value: %sysfunc(strip(%qscan(%superq(rmst_pval&z),&b,|,m)))"
                                                                                            %end;
                                                                                            %else %do;
                                                                                                "RMST Chi-square P-value: %sysfunc(strip(%qscan(%superq(rmst_pval&z),&b,|,m)))"
                                                                                            %end;/ 
                                                                                            drawspace=layoutpercent justify=left y=50 x=0 width=10000 widthunit=percent anchor=left;
                                                                                    %end;
                                                                                    %if %superq(_display_mstats&k._&z)=CINDEX %then %do;
                                                                                        /**Univariate c-index**/
                                                                                        drawtext textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=black ) 
                                                                                            %if %sysevalf(&&cindexheader&z=,boolean)=0 %then %do;
                                                                                                "&&cindexheader&z %sysfunc(strip(%qscan(%superq(cindex&z),&b,|,m)))"
                                                                                            %end;
                                                                                            %else %do;
                                                                                                "C-index (%sysevalf(100-100*&&alpha&z)% CI): %sysfunc(strip(%qscan(%superq(cindex&z),&b,|,m)))"
                                                                                            %end;/ 
                                                                                            drawspace=layoutpercent justify=left y=50 x=0 width=10000 widthunit=percent anchor=left;
                                                                                    %end;
                                                                                    %if %superq(_display_mstats&k._&z)=CINDEXMV %then %do;
                                                                                        /**Multivariate c-index**/
                                                                                        drawtext textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=black ) 
                                                                                            %if %sysevalf(&&cindexmvheader&z=,boolean)=0 %then %do;
                                                                                                "&&cindexmvheader&z %sysfunc(strip(%qscan(%superq(cindexmv&z),&b,|,m)))"
                                                                                            %end;
                                                                                            %else %do;
                                                                                                "Multivariate C-index (%sysevalf(100-100*&&alpha&z)% CI): %sysfunc(strip(%qscan(%superq(cindexmv&z),&b,|,m)))"
                                                                                            %end;/ 
                                                                                            drawspace=layoutpercent justify=left y=50 x=0 width=10000 widthunit=percent anchor=left;
                                                                                    %end;
                                                                                endlayout;
                                                                                %end;
                                                                            endlayout;
                                                                        
                                                                            
                                                                        %end;
                                                                    %end;
                                                                    %let b=%eval(&b+1);
                                                                    %let c_tally=1;
                                                                %end;
                                                                %if %sysevalf(%superq(by&z)^=,boolean) %then %do i = 1 %to %superq(_ndisplay_class_&z);
                                                                    layout gridded / rows=1 columns=1 border=false opaque=false ;
                                                                        /**Draws lines above each BY Group**/
                                                                        %if &&draw_underlines&z=1 %then %do;
                                                                            drawline x1=0 x2=100 y1=100 y2=100 / drawspace=layoutpercent lineattrs=(color=&&ul_color&z pattern=&&ul_pattern&z thickness=&&ul_size&z);
                                                                        %end;
                                                                        %if %superq(_display_class&i._&z)=PVAL_INTER %then %do;
                                                                            /**P-value**/
                                                                            entry halign=center 
                                                                                %if &b=1 and (&c=1 or &c_tally=0) %then %do;
                                                                                    %if %sysevalf(%qscan(&&ipval&sfx.&z,&b,|,m)^=,boolean) %then %do;
                                                                                        "%sysfunc(strip(%qscan(&&ipval&sfx.&z,&b,|,m)))"
                                                                                    %end;
                                                                                    %else %do;
                                                                                        " "
                                                                                    %end;
                                                                                %end;
                                                                                %else %do;
                                                                                    "--"
                                                                                %end;
                                                                                /  valign=top
                                                                                    textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                                    %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                        color=%superq(color&z)
                                                                                    %end;
                                                                                    %else %if %superq(statcolor&z)=1 %then %do;
                                                                                        color=%scan(%superq(color&z), &c)
                                                                                    %end;
                                                                                    %else %do;
                                                                                        color=&fontcolor 
                                                                                    %end;);
                                                                            %let pclass_check=1; 
                                                                            
                                                                        %end;
                                                                        %else %if %superq(_display_class&i._&z)=PVAL_INTERMV %then %do;
                                                                            /**Multivariate P-value**/
                                                                            entry halign=center 
                                                                                %if &b=1 and (&c=1 or &c_tally=0) %then %do;
                                                                                    %if %sysevalf(%qscan(&&ipvalmv&sfx.&z,&b,|,m)^=,boolean) %then %do;
                                                                                        "%sysfunc(strip(%qscan(&&ipvalmv&sfx.&z,&b,|,m)))"
                                                                                    %end;
                                                                                    %else %do;
                                                                                        " "
                                                                                    %end;
                                                                                %end;
                                                                                %else %do;
                                                                                    "--"
                                                                                %end;
                                                                                /  valign=top
                                                                                    textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                                    %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                        color=%superq(color&z)
                                                                                    %end;
                                                                                    %else %if %superq(statcolor&z)=1 %then %do;
                                                                                        color=%scan(%superq(color&z), &c)
                                                                                    %end;
                                                                                    %else %do;
                                                                                        color=&fontcolor 
                                                                                    %end;);
                                                                            %let pclass_check=1; 
                                                                        %end;
                                                                        %else %do;
                                                                            entry " "/ valign=top 
                                                                                textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=white);
                                                                        %end;
                                                                        %if &i=1 %then %do;
                                                                            drawtext textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=black weight=bold) "%superq(by_&z._&b)"/ 
                                                                                drawspace=layoutpercent justify=left y=50 x=0 width=10000 widthunit=percent anchor=left;
                                                                        %end;
                                                                    endlayout;
                                                                %end;
                                                            %end;
                                                            %local pclass_check;
                                                            %do i = 1 %to %superq(_ndisplay_class_&z);
                                                                %let pclass_check=0;
                                                                %if %superq(_display_class&i._&z)=LEGEND %then %do;
                                                                    /**Legend**/
                                                                    /**Each DISCRETELEGEND statement only contains one class level of the plot**/
                                                                    discretelegend "plot&c.&sfx._&z" / opaque=false
                                                                        %if %sysevalf(%superq(legendlinelength&z)^=,boolean) %then %do;
                                                                            itemsize=(linelength=%superq(legendlinelength&z))
                                                                        %end;
                                                                        across=1 down=1 border=false  valign=top halign=center displayclipped=true
                                                                        valueattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                        %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                            color=%superq(color&z)
                                                                        %end;
                                                                        %else %if %superq(statcolor&z)=1 %then %do;
                                                                            color=%scan(%superq(color&z), &c)
                                                                        %end;
                                                                        %else %do;
                                                                            color=&fontcolor 
                                                                        %end;);
                                                                    /**Class Levels**/
                                                                    entry halign=%superq(classvalalign&z) %unquote("%superq(class_&z._&c) &sfx2")/ valign=top 
                                                                        textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                        %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                           color=%superq(color&z)
                                                                       %end;
                                                                       %else %if %superq(statcolor&z)=1 %then %do;
                                                                           color=%scan(%superq(color&z), &c)
                                                                       %end;
                                                                        %else %do;
                                                                            color=&fontcolor 
                                                                        %end;);
                                                                    %let i = %sysevalf(&i+1);
                                                                    %let pclass_check=1;  
                                                                %end;
                                                                %if %superq(pclass_check)=0 %then %do j=1 %to %sysfunc(countw(&classcolumns,|));
                                                                    /**Total|Events|Medians**/
                                                                    %if %superq(_display_class&i._&z)=%scan(&classcolumns,&j,|) %then %do;
                                                                        entry halign=center "%scan(%superq(%scan(&classcolumns,&j,|)&sfx.&z),&c,|)" / valign=top
                                                                            textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                        %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                            color=%superq(color&z)
                                                                        %end;
                                                                        %else %if %superq(statcolor&z)=1 %then %do;
                                                                            color=%scan(%superq(color&z), &c)
                                                                        %end;
                                                                        %else %do;
                                                                            color=&fontcolor 
                                                                        %end;);
                                                                        %let pclass_check=1;  
                                                                    %end;
                                                                %end;
                                                                %if %superq(pclass_check)=0 and %sysevalf(%superq(class&z)^=,boolean) %then %do j=1 %to %sysfunc(countw(&hrcolumns,|));
                                                                    /**Hazard Ratios**/
                                                                    %if %superq(_display_class&i._&z)=%scan(&hrcolumns,&j,|) %then %do;
                                                                        entry halign=center "%scan(%superq(%scan(&hrcolumns,&j,|)&sfx.&z),&c,|)" / valign=top
                                                                            textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                            %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                color=%superq(color&z)
                                                                            %end;
                                                                            %else %if %superq(statcolor&z)=1 %then %do;
                                                                                color=%scan(%superq(color&z), &c)
                                                                            %end;
                                                                            %else %do;
                                                                                color=&fontcolor 
                                                                            %end;);
                                                                    %let pclass_check=1;  
                                                                    %end;
                                                                %end;
                                                                %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=PVAL %then %do;
                                                                    /**P-value**/
                                                                    entry halign=center 
                                                                        %if &c=1 or &c_tally=0 %then %do;
                                                                            "%sysfunc(strip(%qscan(&&pval&sfx.&z,&b,|,m)))"
                                                                        %end;
                                                                        %else %do;
                                                                            "--"
                                                                        %end;
                                                                        /  valign=top
                                                                            textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                            %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                color=%superq(color&z)
                                                                            %end;
                                                                            %else %if %superq(statcolor&z)=1 %then %do;
                                                                                color=%scan(%superq(color&z), &c)
                                                                            %end;
                                                                            %else %do;
                                                                                color=&fontcolor 
                                                                            %end;);
                                                                    %let pclass_check=1; 
                                                                    
                                                                %end;
                                                                %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=PVALMV %then %do;
                                                                    /**Multivariate P-value**/
                                                                    entry halign=center 
                                                                        %if &c=1 or &c_tally=0 %then %do;
                                                                            "%sysfunc(strip(%qscan(&&pvalmv&sfx.&z,&b,|,m)))"
                                                                        %end;
                                                                        %else %do;
                                                                            "--"
                                                                        %end;
                                                                        /  valign=top
                                                                            textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                            %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                color=%superq(color&z)
                                                                            %end;
                                                                            %else %if %superq(statcolor&z)=1 %then %do;
                                                                                color=%scan(%superq(color&z), &c)
                                                                            %end;
                                                                            %else %do;
                                                                                color=&fontcolor 
                                                                            %end;);
                                                                    %let pclass_check=1; 
                                                                %end;
                                                                %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=PVAL_INTER %then %do;
                                                                    /**P-value**/
                                                                    entry halign=center 
                                                                        /*%if &b=1 and (&c=1 or &c_tally=0) %then %do;
                                                                            "%sysfunc(strip(%qscan(&&ipval&sfx.&z,&b,|,m)))"
                                                                        %end;
                                                                        %else %do;
                                                                            "--"
                                                                        %end;*/
                                                                        "  "
                                                                        /  valign=top
                                                                            textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                            %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                color=%superq(color&z)
                                                                            %end;
                                                                            %else %if %superq(statcolor&z)=1 %then %do;
                                                                                color=%scan(%superq(color&z), &c)
                                                                            %end;
                                                                            %else %do;
                                                                                color=&fontcolor 
                                                                            %end;);
                                                                    %let pclass_check=1; 
                                                                    
                                                                %end;
                                                                %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=PVAL_INTERMV %then %do;
                                                                    /**Multivariate P-value**/
                                                                    entry halign=center 
                                                                        /*%if &b=1 and (&c=1 or &c_tally=0) %then %do;
                                                                            "%sysfunc(strip(%qscan(&&ipvalmv&sfx.&z,&b,|,m)))"
                                                                        %end;
                                                                        %else %do;
                                                                            "--"
                                                                        %end;*/
                                                                        "  "
                                                                        /  valign=top
                                                                            textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                            %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                color=%superq(color&z)
                                                                            %end;
                                                                            %else %if %superq(statcolor&z)=1 %then %do;
                                                                                color=%scan(%superq(color&z), &c)
                                                                            %end;
                                                                            %else %do;
                                                                                color=&fontcolor 
                                                                            %end;);
                                                                    %let pclass_check=1; 
                                                                %end;
                                                                %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=RMST_PVAL %then %do;
                                                                    /**Restricted Means Survival Time P-value comparison**/
                                                                    entry halign=center 
                                                                        %if &c=1 or &c_tally=0 %then %do;
                                                                            "%sysfunc(strip(%qscan(%superq(rmst_pval&z),&b,|,m)))"
                                                                        %end;
                                                                        %else %do;
                                                                            "--"
                                                                        %end;
                                                                        /  valign=top
                                                                            textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                            %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                color=%superq(color&z)
                                                                            %end;
                                                                            %else %if %superq(statcolor&z)=1 %then %do;
                                                                                color=%scan(%superq(color&z), &c)
                                                                            %end;
                                                                            %else %do;
                                                                                color=&fontcolor 
                                                                            %end;);
                                                                    %let pclass_check=1; 
                                                                %end;
                                                                %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=CINDEX %then %do;
                                                                    /**Univariate c-index**/
                                                                    entry halign=center 
                                                                        %if &c=1 or &c_tally=0 %then %do;
                                                                            "%sysfunc(strip(%qscan(%superq(cindex&z),&b,|,m)))"
                                                                        %end;
                                                                        %else %do;
                                                                            "--"
                                                                        %end;
                                                                        /  valign=top
                                                                            textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                            %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                color=%superq(color&z)
                                                                            %end;
                                                                            %else %if %superq(statcolor&z)=1 %then %do;
                                                                                color=%scan(%superq(color&z), &c)
                                                                            %end;
                                                                            %else %do;
                                                                                color=&fontcolor 
                                                                            %end;);
                                                                    %let pclass_check=1; 
                                                                %end;
                                                                %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=CINDEXMV %then %do;
                                                                    /**Multivariate c-index**/
                                                                    entry halign=center 
                                                                        %if &c=1 or &c_tally=0 %then %do;
                                                                            "%sysfunc(strip(%qscan(%superq(cindexmv&z),&b,|,m)))"
                                                                        %end;
                                                                        %else %do;
                                                                            "--"
                                                                        %end;
                                                                        /  valign=top
                                                                            textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                            %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                color=%superq(color&z)
                                                                            %end;
                                                                            %else %if %superq(statcolor&z)=1 %then %do;
                                                                                color=%scan(%superq(color&z), &c)
                                                                            %end;
                                                                            %else %do;
                                                                                color=&fontcolor 
                                                                            %end;);
                                                                    %let pclass_check=1; 
                                                                %end;
                                                                %if %superq(pclass_check)=0 and %superq(_display_class&i._&z)=TIMELIST %then %do;
                                                                    /**Survival Time-point Estimates**/                                 
                                                                    %if %sysevalf(%superq(listtimepoints&z)=1,boolean) %then %do;
                                                                        /*Time-point labels*/
                                                                        layout gridded / columns=1 rows=%superq(ntl_&z) border=false opaque=false;
                                                                            %do k = 1 %to %superq(ntl_&z);
                                                                                entry halign=center "%scan(%superq(timelist&sfx.v&z),&k,|)" / valign=top
                                                                                    textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                                    %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                        color=%superq(color&z)
                                                                                    %end;
                                                                                    %else %if %superq(statcolor&z)=1 %then %do;
                                                                                        color=%scan(%superq(color&z), &c)
                                                                                    %end;
                                                                                    %else %do;
                                                                                        color=&fontcolor 
                                                                                    %end;);
                                                                            %end;
                                                                        endlayout;
                                                                    %end;
                                                                    /*Survival Estimates*/
                                                                    layout gridded / columns=1 rows=%superq(ntl_&z) border=false opaque=false;
                                                                        %do k = 1 %to %superq(ntl_&z);
                                                                            entry halign=center "%scan(%superq(timelist&sfx._&k._&z),&c,|)" / valign=top
                                                                                textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)"
                                                                                %if %superq(statcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                                    color=%superq(color&z)
                                                                                %end;
                                                                                %else %if %superq(statcolor&z)=1 %then %do;
                                                                                    color=%scan(%superq(color&z), &c)
                                                                                %end;
                                                                                %else %do;
                                                                                    color=&fontcolor 
                                                                                %end;);
                                                                        %end;
                                                                    endlayout;
                                                                    %let i = %sysevalf(&i+%superq(listtimepoints&z));
                                                                    %let pclass_check=1; 
                                                                %end;/**Ends Time-point estimate Section**/                               
                                                            %end;/**Ends column filling loop for class-level gridded block***/
                                                        %end;/**Ends row filling loop for class-level gridded block**/  
                                                     %end; /**End adjusting method potential block**/          
                                                endlayout;/**Ends class-level gridded block**/
                                            %end;/**Ends the class-level gridded block**/
                                            /**Creates model-level gridded block**/
                                            %if %superq(_ndisplay_model_&z)>0 %then %do;
                                                layout gridded / rows=1 columns=%superq(_ndisplay_model_&z) opaque=false
                                                    border=false valign=top halign=right;
                                                    %if %sysevalf(%superq(by&z)=,boolean) %then %_underline_(y=105);
                                                    %do i = 1 %to %superq(_ndisplay_model_&z);                            
                                                        /**Uses the plot that was colored white earlier to create the correct amount of space on the left**/
                                                        %if %superq(_display_model&i._&z)=LEGEND %then %do;
                                                            discretelegend "spacer" / opaque=false 
                                                                %if %sysevalf(%superq(legendlinelength&z)^=,boolean) %then %do;
                                                                    itemsize=(linelength=%superq(legendlinelength&z))
                                                                %end;
                                                                across=1 down=1 border=false valign=top halign=left displayclipped=true
                                                                valueattrs=(color=&fontcolor size=1pt family="%superq(ptabfamily&z)");
                                                        %end;
                                                        %if %superq(_display_model&i._&z)=STATS and %sysevalf(%superq(_ndisplay_mstats_&z)>0,boolean) %then %do;
                                                            /**Print text for P-values and user entered Table Comments**/                           
                                                            layout gridded / columns=1 
                                                                %if %qupcase(&&method&z)=KM or %qupcase(&&method&z)=CIF %then %do;
                                                                    rows=%superq(_ndisplay_mstats_&z)                                         
                                                                %end;
                                                                %else %do;
                                                                    rows=%sysevalf((1+&&plot_unadjust&z)*%superq(_ndisplay_mstats_&z)) /**One row per class level**/ 
                                                                %end;
                                                                valign=top halign=left
                                                                opaque=false border=false;
                                                                %local c adj sfx sfx2;
                                                                %if %sysevalf(%superq(ntl_&z)=,boolean) %then %let ntl_&z = 1;
                                                                %do adj = (&&plot_unadjust&z=0) %to (%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT);
                                                                    %if &adj=1 %then %do;
                                                                        %let sfx=mv;
                                                                        %if &&plot_unadjust&z=1 %then %let sfx2=%str(Adjusted );
                                                                        %else %let sfx2=;
                                                                    %end;
                                                                    %else %do;
                                                                        %let sfx=;
                                                                        %let sfx2=;
                                                                    %end;
                                                                    %do k = 1 %to %superq(_ndisplay_mstats_&z);
                                                                        %if %superq(_display_mstats&k._&z)=PVAL %then %do;
                                                                            /**P-value**/
                                                                            %if %qupcase(%superq(plotpval&z))=LR %then %let plotpval&z=Likelihood-Ratio; 
                                                                            %else %if %qupcase(%superq(plotpval&z))=GRAY %then %let plotpval&z=Gray K-Sample Test; 
                                                                            %else %if %qupcase(%superq(plotpval&z))=LOGRANK_ONESIDED %then %let plotpval&z=One-sided logrank; 
                                                                            %if %qupcase(%superq(plotpvalmv&z))=LR %then %let plotpvalmv&z=Likelihood-Ratio; 
                                                                            entry halign=left
                                                                                %if %sysevalf(%superq(pval&sfx.header&z)=,boolean)=0 %then %do;
                                                                                    "%unquote(&&pval&sfx.header&z) %sysfunc(strip(%qscan(&&pval&sfx.&z,&b,|,m)))"
                                                                                %end;
                                                                                %else %if %sysevalf(%superq(strata&z)^=,boolean) %then %do;
                                                                                    "Stratified %sysfunc(propcase(%superq(plotpval&sfx.&z))) P-value: %sysfunc(strip(%qscan(&&pval&sfx.&z,&b,|,m)))"
                                                                                %end;
                                                                                %else %do;
                                                                                    "&sfx2.%sysfunc(propcase(%superq(plotpval&sfx.&z))) P-value: %sysfunc(strip(%qscan(&&pval&sfx.&z,&b,|,m)))"
                                                                                %end; / valign=top textattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                        %end;
                                                                        %if %superq(_display_mstats&k._&z)=PVALMV %then %do;
                                                                            /**P-value**/
                                                                            %if %qupcase(%superq(plotpvalmv&z))=LR %then %let plotpvalmv&z=Likelihood-Ratio; 
                                                                            entry halign=left
                                                                            %if %sysevalf(&&pvalmvheader&z=,boolean)=0 %then %do;
                                                                                "&&pvalmvheader&z %sysfunc(strip(%qscan(&&pvalmv&sfx.&z,&b,|,m)))"
                                                                            %end;
                                                                            %else %do;
                                                                                "Adjusted %sysfunc(propcase(%superq(plotpvalmv&z))) P-value: %sysfunc(strip(%qscan(&&pvalmv&sfx.&z,&b,|,m)))"
                                                                            %end; / valign=top textattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                        %end;
                                                                        %if %superq(_display_mstats&k._&z)=RMST_PVAL %then %do;
                                                                            /**Restricted Means Survival Time P-value comparison**/
                                                                            entry halign=left
                                                                            %if %sysevalf(&&rmst_pvalheader&z=,boolean)=0 %then %do;
                                                                                "&&rmst_pvalheader&z %sysfunc(strip(%qscan(%superq(rmst_pval&z),&b,|,m)))"
                                                                            %end;
                                                                            %else %if %sysevalf(%superq(strata&z)^=,boolean) %then %do;
                                                                                "Stratified RMST Chi-square P-value: %sysfunc(strip(%qscan(%superq(rmst_pval&z),&b,|,m)))"
                                                                            %end;
                                                                            %else %do;
                                                                                "RMST Chi-square P-value: %sysfunc(strip(%qscan(%superq(rmst_pval&z),&b,|,m)))"
                                                                            %end; / valign=top textattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                        %end;
                                                                        %if %superq(_display_mstats&k._&z)=CINDEX %then %do;
                                                                            /**Univariate c-index**/
                                                                            entry halign=left
                                                                            %if %sysevalf(&&cindexheader&z=,boolean)=0 %then %do;
                                                                                "&&cindexheader&z %sysfunc(strip(%qscan(%superq(cindex&z),&b,|,m)))"
                                                                            %end;
                                                                            %else %do;
                                                                                "C-index (%sysevalf(100-100*&&alpha&z)% CI): %sysfunc(strip(%qscan(%superq(cindex&z),&b,|,m)))"
                                                                            %end; / valign=top textattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                        %end;
                                                                        %if %superq(_display_mstats&k._&z)=CINDEXMV %then %do;
                                                                            /**Multivariate c-index**/
                                                                            entry halign=left
                                                                            %if %sysevalf(&&cindexmvheader&z=,boolean)=0 %then %do;
                                                                                "&&cindexmvheader&z %sysfunc(strip(%qscan(%superq(cindexmv&z),&b,|,m)))"
                                                                            %end;
                                                                            %else %do;
                                                                                "Multivariate C-index (%sysevalf(100-100*&&alpha&z)% CI): %sysfunc(strip(%qscan(%superq(cindexmv&z),&b,|,m)))"
                                                                            %end; / valign=top textattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                        %end;
                                                                        %if %superq(_display_mstats&k._&z)=TABLECOMMENTS %then %do;
                                                                            /**User-provided Table Comments**/
                                                                            %do j=1 %to %sysfunc(countw(%superq(tablecomments&z),`,m));                   
                                                                                entry halign=left "%scan(&&tablecomments&z,&j,`,m)" / valign=top border=false
                                                                                    textattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                            %end;
                                                                            %let k=%sysevalf(&k+%sysfunc(countw(%superq(tablecomments&z),`,m))-1);
                                                                        %end;
                                                                    %end;
                                                                %end;
                                                            endlayout;
                                                        %end;
                                                        %if %superq(_display_model&i._&z)=CENSORS %then %do;
                                                            /**Create legend statement for censor values**/
                                                            discretelegend "cens&z" 
                                                                / border=false halign=right valign=top displayclipped=true opaque=false autoitemsize=TRUE
                                                                valueattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");   
                                                        %end;   
                                                    %end;
                                                endlayout;/**Ends model-level gridded block**/
                                            %end;/**Ends model-level gridded block**/
                                            /**Creates Interaction Test-level gridded block**/
                                            %if %superq(_ndisplay_interp_&z)>0 %then %do;
                                                layout gridded / rows=1 columns=%superq(_ndisplay_interp_&z) opaque=false
                                                    border=false valign=top halign=right;
                                                    /**Draws line above bottom section**/
                                                    %if &&draw_underlines&z=1 %then %do;
                                                        drawline x1=0 x2=100 y1=105 y2=105 / drawspace=layoutpercent lineattrs=(color=&&ul_color&z pattern=&&ul_pattern&z thickness=&&ul_size&z);
                                                    %end;
                                                    %do i = 1 %to %superq(_ndisplay_interp_&z);                            
                                                        /**Uses the plot that was colored white earlier to create the correct amount of space on the left**/
                                                        %if %superq(_display_interp&i._&z)=STATS and %sysevalf(%superq(_ndisplay_istats_&z)>0,boolean) %then %do;
                                                            /**Print text for P-values and user entered Table Comments**/                           
                                                            layout gridded / columns=1 
                                                                %if %qupcase(&&method&z)=KM or %qupcase(&&method&z)=CIF %then %do;
                                                                    rows=%superq(_ndisplay_istats_&z)                                         
                                                                %end;
                                                                %else %do;
                                                                    rows=%sysevalf((1+&&plot_unadjust&z)*%superq(_ndisplay_istats_&z)) /**One row per class level**/ 
                                                                %end;
                                                                valign=top halign=left
                                                                opaque=false border=false;
                                                                
                                                                %local c adj sfx sfx2;
                                                                %if %sysevalf(%superq(ntl_&z)=,boolean) %then %let ntl_&z = 1;
                                                                %do adj = (&&plot_unadjust&z=0) %to (%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT);
                                                                    %if &adj=1 %then %do;
                                                                        %let sfx=mv;
                                                                        %if &&plot_unadjust&z=1 %then %let sfx2=%str(Adjusted );
                                                                        %else %let sfx2=;
                                                                    %end;
                                                                    %else %do;
                                                                        %let sfx=;
                                                                        %let sfx2=;
                                                                    %end;
                                                                    %do k = 1 %to %superq(_ndisplay_istats_&z);
                                                                        %if %superq(_display_istats&k._&z)=PVAL_INTER %then %do;
                                                                            /**P-value**/
                                                                            %if %qupcase(%superq(pval_inter&sfx.&z))=LR %then %let pval_inter&sfx.&z=Likelihood-Ratio; 
                                                                            entry halign=left
                                                                                %if %sysevalf(%superq(pval_inter&sfx.header&z)=,boolean)=0 %then %do;
                                                                                    "%unquote(&&pval_inter&sfx.header&z) &&ipval&sfx.&z"
                                                                                %end;
                                                                                %else %if %sysevalf(%superq(strata&z)^=,boolean) %then %do;
                                                                                    "Stratified %sysfunc(propcase(%superq(pval_inter&sfx.&z))) Interaction P-value: &&ipval&sfx.&z"
                                                                                %end;
                                                                                %else %do;
                                                                                    "&sfx2.%sysfunc(propcase(%superq(pval_inter&sfx.&z))) Interaction P-value: &&ipval&sfx.&z"
                                                                                %end; / valign=top textattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                        %end;
                                                                        %if %superq(_display_istats&k._&z)=PVAL_INTERMV %then %do;
                                                                            /**P-value**/
                                                                            %if %qupcase(%superq(pval_intermv&z))=LR %then %let pval_intermv&z=Likelihood-Ratio; 
                                                                            entry halign=left
                                                                                %if %sysevalf(%superq(pval_intermv&sfx.header&z)=,boolean)=0 %then %do;
                                                                                    "%unquote(&&pval_intermv&sfx.header&z) &&ipvalmv&sfx.&z"
                                                                                %end;
                                                                                %else %if %sysevalf(%superq(strata&z)^=,boolean) %then %do;
                                                                                    "Adjusted Stratified %sysfunc(propcase(%superq(pval_intermv&sfx.&z))) interaction P-value: &&ipvalmv&sfx.&z"
                                                                                %end;
                                                                                %else %do;
                                                                                    "Adjusted &sfx2.%sysfunc(propcase(%superq(pval_intermv&sfx.&z))) interaction P-value: &&ipvalmv&sfx.&z"
                                                                                %end; / valign=top textattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                        %end;
                                                                        %if %superq(_display_istats&k._&z)=TABLECOMMENTS %then %do;
                                                                            /**User-provided Table Comments**/
                                                                            %do j=1 %to %sysfunc(countw(%superq(tablecomments&z),`,m));                   
                                                                                entry halign=left "%scan(&&tablecomments&z,&j,`,m)" / valign=top border=false
                                                                                    textattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");
                                                                            %end;
                                                                            %let k=%sysevalf(&k+%sysfunc(countw(%superq(tablecomments&z),`,m))-1);
                                                                        %end;
                                                                    %end;
                                                                %end;
                                                            endlayout;
                                                        %end;
                                                        %if %superq(_display_interp&i._&z)=CENSORS %then %do;                      
                                                            layout gridded / columns=2 rows=1
                                                                valign=top halign=right
                                                                opaque=false border=false;
                                                                entry halign=left " " / valign=top textattrs=(color=&fontcolor size=1pt family="%superq(ptabfamily&z)");
                                                                /**Create legend statement for censor values**/
                                                                discretelegend "cens&z" 
                                                                    / border=false halign=right valign=top displayclipped=true opaque=false autoitemsize=TRUE
                                                                    valueattrs=(color=&fontcolor size=%superq(ptabsize&z) family="%superq(ptabfamily&z)");   
                                                            endlayout;
                                                        %end;   
                                                    %end;
                                                endlayout;/**Ends model-level gridded block**/
                                            %end;/**Ends model-level gridded block**/
                                        endlayout;/**Ends the outer gridded block**/ 
                                    %end;
                                endlayout;/**Closes the LAYOUT OVERLAY**/
                                /**Creates the patients-at-risk block**/
                                /**Makes one block plot per class level**/
                                %if %sysevalf(%superq(risklist&z)=,boolean)=0 and %qupcase(%superq(risklocation&z))=BOTTOM %then %do k=1 %to %sysfunc(countw(%superq(pardisplay&z),%str( )));
                                /**Makes a header for the Patients-at-Risk table**/
                                    %if (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) and %sysevalf(%superq(parheader&z)^=,boolean)) or
                                     (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) and %sysevalf(%superq(ncensheader&z)^=,boolean)) or
                                     (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) and %sysevalf(%superq(neventsheader&z)^=,boolean)) or
                                     (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) and 
                                        (%sysevalf(%superq(parheader&z)^=,boolean) and %sysevalf(%superq(ncensheader&z)^=,boolean))) or
                                     (%sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) and 
                                        (%sysevalf(%superq(parheader&z)^=,boolean) and %sysevalf(%superq(neventsheader&z)^=,boolean))) %then %do;
                                        layout overlay / border=false walldisplay=none
                                            xaxisopts=(display=none type=linear                 
                                                /**Offset creates space at the top or bottom of the window that the plot cannot use, proportion from 0 to 1**/
                                                %if %sysevalf(%superq(xmaxoffset&z)^=,boolean) %then %do; offsetmax=%superq(xmaxoffset&z)%end;
                                                %if %sysevalf(%superq(xminoffset&z)^=,boolean) %then %do; offsetmin=%superq(xminoffset&z)%end;
                                                %else %if %sysevalf(%superq(risklist&z)=,boolean) or %sysevalf(%qupcase(&&risklocation&z)=BOTTOM,boolean) %then %do; offsetmin=0 %end;
                                                /**VIEWMAX and VIEWMIN are also required to show the desired range**/
                                                linearopts=(viewmin=%superq(xmin&z) viewmax=%superq(xmax&z)));
                                            blockplot x=eval(ifn(^missing(partitle_&z),%superq(xmin&z),.)) 
                                                %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) %then %do;
                                                    block=eval(tranwrd(repeat('a',length(partitle_&z )),'a','A0'x))
                                                %end;
                                                %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) %then %do;
                                                    block=eval(tranwrd(repeat('a',length(ncenstitle_&z  )),'a','A0'x))
                                                %end;
                                                %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) %then %do;
                                                    block=eval(tranwrd(repeat('a',length(neventstitle_&z )),'a','A0'x))
                                                %end;
                                                %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) %then %do;
                                                    block=eval(tranwrd(repeat('a',length(strip(partitle_&z)||' ('||strip(ncenstitle_&z)||')')),'a','A0'x))
                                                %end; 
                                                %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) %then %do;
                                                    block=eval(tranwrd(repeat('a',length(strip(partitle_&z)||' ('||strip(neventstitle_&z)||')')),'a','A0'x))
                                                %end;
                                                /
                                                %if %sysevalf(%qupcase(%superq(paralign&z))=LABELS,boolean) %then %do;
                                                    display=(values) 
                                                    %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) %then %do;
                                                        label="%superq(parheader&z)" 
                                                    %end;
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) %then %do;
                                                        label="%superq(ncensheader&z)"
                                                    %end;
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) %then %do;
                                                        label="%superq(neventsheader&z)"
                                                    %end;
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) %then %do;
                                                        label="%superq(parheader&z) (%superq(ncensheader&z))"
                                                    %end; 
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) %then %do;
                                                        label="%superq(parheader&z) (%superq(neventsheader&z))"
                                                    %end;
                                                    labelattrs=(color=&fontcolor size=%superq(parsize&z) weight=%superq(parweight&z) family="%superq(parfamily&z)")
                                                %end;
                                                %else %do;
                                                    display=(values)
                                                    valuehalign=%superq(paralign&z) 
                                                %end;
                                                valueattrs=(size=%superq(parsize&z) color=&background weight=%superq(parweight&z) family="%superq(parfamily&z)");
                                            
                                            %if %sysevalf(%qupcase(%superq(paralign&z))^=LABELS,boolean) %then %do;
                                                /*Draws Labels*/
                                                drawtext textattrs=(size=%superq(parsize&z) weight=%superq(parweight&z) family="%superq(parfamily&z)" color=&fontcolor)
                                                    %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) %then %do;
                                                        "%superq(parheader&z)"
                                                    %end;
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) %then %do;
                                                        "%superq(ncensheader&z)"
                                                    %end;
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) %then %do;
                                                        "%superq(neventsheader&z)"
                                                    %end;
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) %then %do;
                                                        "%superq(parheader&z) (%superq(ncensheader&z))"
                                                    %end; 
                                                    %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) %then %do;
                                                        "%superq(parheader&z) (%superq(neventsheader&z))"
                                                    %end;/
                                                    y=50  
                                                    %if %sysevalf(%qupcase(%superq(paralign&z))=LEFT,boolean) %then %do;
                                                        x=0 anchor=left justify=left
                                                    %end;
                                                    %else %if %sysevalf(%qupcase(%superq(paralign&z))=RIGHT,boolean) %then %do;
                                                        x=100 anchor=right justify=right
                                                    %end;
                                                    %else %if %sysevalf(%qupcase(%superq(paralign&z))=CENTER,boolean) %then %do;
                                                        x=50 anchor=center justify=center
                                                    %end; xspace=wallpercent yspace=layoutpercent width=10000 ;
                                            %end;
                                        endlayout;
                                    %end;
                                    %local adj;
                                    %do adj=(&&plot_unadjust&z=0) %to %sysevalf((%qupcase(&&method&z)=INVWTS or %qupcase(&&method&z)=DIRECT));
                                        %local sfx sfx2;
                                        %if &adj=0 %then %do;
                                            %let sfx=;%let sfx2=;
                                        %end;
                                        %else %do;
                                            %let sfx=_adj;
                                            %if &&plot_unadjust&z=1 %then %let sfx2=%str((Adjusted));
                                            %else %let sfx2=;
                                        %end;

                                        %let b=1;%let c_tally=0;
                                        %do i = 1 %to %superq(nclass_&z);     
                                            %let c_tally=%sysevalf(&c_tally + 1);
                                            %if %sysevalf(%superq(by&z)^=,boolean) and ((&b=1 and &i=1) or (&c_tally>%superq(nclass_&b._&z.))) %then %do;
                                                %if &c_tally>%superq(nclass_&b._&z.) %then %do;
                                                    %let b=%sysevalf(&b+1);
                                                    %let c_tally=0;
                                                %end;
                                                layout gridded / rows=1 columns=1 border=false opaque=false ;
                                                    /**Draws line for risk list when set to bottom**/
                                                    %if &&draw_underlines&z=1 %then %do;
                                                        drawline x1=0 x2=100 y1=100 y2=100 / drawspace=layoutpercent lineattrs=(color=&&ul_color&z pattern=&&ul_pattern&z thickness=&&ul_size&z);
                                                    %end;
                                                    entry " "/ valign=top 
                                                        textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=white);
                                                endlayout;
                                            %end;
                                        
                                            %if %qupcase(%superq(risklabellocation&z))=ABOVE and
                                                (%sysevalf(%superq(class&z)^=,boolean) or %sysevalf(%superq(classdesc&z)^=,boolean)) %then %do;
                                                layout overlay / border=false walldisplay=none
                                                    yaxisopts=(display=none)
                                                    xaxisopts=(display=none  type=linear                
                                                        /**Offset creates space at the top or bottom of the window that the plot cannot use, proportion from 0 to 1**/
                                                        %if %sysevalf(%superq(xmaxoffset&z)^=,boolean) %then %do; offsetmax=%superq(xmaxoffset&z)%end;
                                                        %if %sysevalf(%superq(xminoffset&z)^=,boolean) %then %do; offsetmin=%superq(xminoffset&z)%end;
                                                        %else %if %sysevalf(%superq(risklist&z)=,boolean) or %sysevalf(%qupcase(&&risklocation&z)=BOTTOM,boolean) %then %do; offsetmin=0 %end;
                                                        /**VIEWMAX and VIEWMIN are also required to show the desired range**/
                                                        linearopts=(viewmin=%superq(xmin&z) viewmax=%superq(xmax&z)));
                                                    scatterplot x=eval(t_&z*0) y=eval(t_&z*0) / markerattrs=(size=0pt);
                                                    /*Draws Labels*/
                                                    drawtext textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" 
                                                        %if %superq(riskcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                            color=%superq(color&z)
                                                        %end;
                                                        %else %if %superq(riskcolor&z)=1 %then %do;
                                                            color=%scan(%superq(color&z), &i)
                                                        %end;
                                                        %else %do;
                                                             color=&fontcolor
                                                        %end;) 
                                                        %unquote("%superq(class_&z._&i) &sfx2.") /
                                                        y=50  
                                                        %if %sysevalf(%qupcase(%superq(risklabelalign&z))=LEFT,boolean) %then %do;
                                                            x=0 anchor=left justify=left
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%superq(risklabelalign&z))=RIGHT,boolean) %then %do;
                                                            x=100 anchor=right justify=right
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%superq(risklabelalign&z))=CENTER,boolean) %then %do;
                                                            x=50 anchor=center justify=center
                                                        %end; xspace=datapercent yspace=layoutpercent width=10000 ;
                                                endlayout;
                                            %end;
                                            layout overlay / border=false walldisplay=none  
                                                yaxisopts=(display=none)
                                                xaxisopts=(display=none  type=linear  
                                                    /**Offset creates space at the top or bottom of the window that the plot cannot use, proportion from 0 to 1**/
                                                    %if %sysevalf(%superq(xmaxoffset&z)^=,boolean) %then %do; offsetmax=%superq(xmaxoffset&z)%end;
                                                    %if %sysevalf(%superq(xminoffset&z)^=,boolean) %then %do; offsetmin=%superq(xminoffset&z)%end;
                                                    %else %if %sysevalf(%superq(risklist&z)=,boolean) or %sysevalf(%qupcase(&&risklocation&z)=BOTTOM,boolean) %then %do; offsetmin=0 %end;
                                                    /**VIEWMAX and VIEWMIN are also required to show the desired range**/
                                                    linearopts=(tickvaluesequence=(start=%superq(xmin&z) end=%superq(xmax&z) increment=%superq(xincrement&z)) viewmin=%superq(xmin&z) viewmax=%superq(xmax&z)));
                                                %if &i=1 and %sysevalf(%superq(by&z)=,boolean) %then %_underline_(y=105);
                                                    blockplot x=eval(ifn(time&i.&sfx._&z ge %superq(xmin&z),time&i.&sfx._&z,.)) 
                                                        %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) %then %do;
                                                            block=eval(tranwrd(repeat('a',length(atrisk&i.&sfx._&z)),'a','A0'x))
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) %then %do;
                                                            block=eval(tranwrd(repeat('a',length(ncens&i.&sfx._&z)),'a','A0'x))
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) %then %do;
                                                            block=eval(tranwrd(repeat('a',length(nevent&i.&sfx._&z)),'a','A0'x))
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) %then %do;
                                                            block=eval(tranwrd(repeat('a',length(strip(atrisk&i.&sfx._&z)||' ('||strip(ncens&i.&sfx._&z)||')')),'a','A0'x))
                                                        %end;
                                                        %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) %then %do;
                                                            block=eval(tranwrd(repeat('a',length(strip(atrisk&i.&sfx._&z)||' ('||strip(ncens&i.&sfx._&z)||')')),'a','A0'x))
                                                        %end;
                                                        /repeatedvalues=TRUE display=(values )
                                                        valuehalign=start
                                                        valueattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=&background)
                                                        labelattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" color=&fontcolor);
        
                                                    /*Draws patients at risk*/
                                                    %if %sysevalf(%superq(risklist_v_&i.&sfx._&z)^=,boolean) %then %do j = 1 %to %sysfunc(countw(%superq(risklist_v_&i.&sfx._&z),|,m));
                                                        drawtext textattrs=(size=%superq(ptabsize&z) family="%superq(ptabfamily&z)" 
                                                            %if %superq(riskcolor&z)=1 and %sysfunc(countw(%superq(color&z))) = 1 %then %do;
                                                                color=%superq(color&z)
                                                            %end;
                                                            %else %if %superq(riskcolor&z)=1 %then %do;
                                                                color=%scan(%superq(color&z), &i)
                                                            %end;
                                                            %else %do;
                                                                 color=&fontcolor
                                                            %end;) 
                                                            %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR,boolean) %then %do;
                                                                "%scan(%superq(risklist_v_&i.&sfx._&z),&j,|,m)"
                                                            %end;
                                                            %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NCENS,boolean) %then %do;
                                                                "%scan(%superq(risklist_c_&i.&sfx._&z),&j,|,m)"
                                                            %end;
                                                            %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=NEVENTS,boolean) %then %do;
                                                                "%scan(%superq(risklist_e_&i.&sfx._&z),&j,|,m)"
                                                            %end;
                                                            %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NCENS,boolean) %then %do;
                                                                "%scan(%superq(risklist_v_&i.&sfx._&z),&j,|,m) (%scan(%superq(risklist_c_&i.&sfx._&z),&j,|,m))"
                                                            %end;
                                                            %else %if %sysevalf(%qupcase(%scan(%superq(pardisplay&z),&k,%str( )))=PAR_NEVENTS,boolean) %then %do;
                                                                "%scan(%superq(risklist_v_&i.&sfx._&z),&j,|,m) (%scan(%superq(risklist_e_&i.&sfx._&z),&j,|,m))"
                                                            %end;
                                                            /
                                                            x=%scan(%superq(risklist_t_&i.&sfx._&z),&j,|,m) y=50  xspace=datavalue yspace=layoutpercent width=10000 justify=center anchor=center;
                                                    %end;
                                            endlayout;
                                        %end;
                                    %end;     
                                %end;
                            endlayout;/**Ends inner lattice block**/    
                        endlayout; 
                    %end;/**Ends Model-by-model loop**/
                endlayout;/**Ends outer lattice block**/  
            endgraph;
        end;
    run;   
        
    /**Turn Results and ODS back on**/
    ods select all;
    ods results;                  
    /**Creates document to save**/
    %if %sysevalf(%superq(outdoc)=,boolean)=0 %then %do;
        ods escapechar='^';
        /**Sets up DPI and ODS generated file**/
        ods &destination 
            %if %qupcase(&destination)=RTF %then %do; 
                file="&outdoc"
                image_dpi=&dpi startpage=NO 
            %end;
            %else %if %qupcase(&destination)=HTML %then %do; 
                image_dpi=&dpi 
                %if %upcase(&sysscpl)=LINUX or %upcase(&sysscpl)=UNIX %then %do;
                    path="%substr(&outdoc,1,%sysfunc(find(&outdoc,/,-%sysfunc(length(&outdoc)))))"
                    file="%scan(&outdoc,1,/,b)"
                %end;
                %else %do;
                    path="%substr(&outdoc,1,%sysfunc(find(&outdoc,\,-%sysfunc(length(&outdoc)))))"
                    file="%scan(&outdoc,1,\,b)"
                %end;
                %if %sysevalf(%superq(gpath)=,boolean)=0 %then %do;
                    gpath="&gpath" (url=none)
                %end;
            %end;
            %else %if %qupcase(&destination)=PDF %then %do; 
                dpi=&dpi startpage=NO bookmarkgen=off notoc
                file="&outdoc"
            %end;
            %else %if %qupcase(&destination)=EXCEL %then %do; 
                file="&outdoc"
                dpi=&dpi options(sheet_interval='none') 
            %end;
            %else %if %qupcase(&destination)=POWERPOINT %then %do; 
                file="&outdoc"
                dpi=&dpi 
            %end;;
    %end;
    %else %if &_listing^=1 %then %do;
        ods listing close image_dpi=&dpi;
    %end;
    %else %if &_listing=1 %then %do;
        ods listing image_dpi=&dpi;
    %end;
    proc sql noprint;
        %local _ppt _other _destinations _styles k pfoot_list;
        select max(ifn(upcase(destination) ^in('LISTING' 'OUTPUT' 'POWERPOINT' 'RTF'),1,0)),
            max(ifn(upcase(destination) in('POWERPOINT'),1,0)),
            max(ifn(upcase(destination) in('RTF'),1,0))
            into :_other separated by '',:_ppt separated by '',:_rtf separated by '' from sashelp.vdest;
        select upcase(destination),upcase(style) into :_destinations separated by '|',:_styles separated by '|'
            from sashelp.vdest
            where upcase(destination)^in('OUTPUT');
    quit;
    /**Create plot if flagged**/
    %if &plot = 1 %then %do;
        /**Names and formats the image**/
        %if %sysevalf(%superq(plottype)^=,boolean) %then %do; 
            %if %qupcase(&plottype)=EMF or (&svg=1 and %qupcase(&destination)=RTF and %qupcase(&plottype)^=TIFF and %qupcase(&plottype)^=TIF)  
                or (&svg=1 and %qupcase(&destination)=EXCEL)
                or (&svg=1 and %qupcase(&destination)=POWERPOINT) %then %do;
                %local _any_trans;
                %let _any_trans=0;
                %do i =1 %to &nmodels;
                    %if (%sysevalf(%superq(class&i)=,boolean) and %sysevalf(%superq(plotci&i)^=0,boolean) and %sysevalf(%superq(plotcifill&i)=1,boolean)) or 
                        (%sysevalf(%superq(plotci&i)=1,boolean) and %sysevalf(%superq(plotcifill&i)=1,boolean)) %then %let _any_trans=1;
                %end;
                options printerpath='emf';
                ods graphics / imagefmt=&plottype;  
                %if &transparent=0 and &_any_trans=0 %then %do;
                    /**Modifies temporary registry keys to create better EMF image in 9.4**/
                    /**Taken from SAS Technical Support Martin Mincey**/
                    %local workdir;
                    %let workdir=%trim(%sysfunc(pathname(work))); 
                    /**Creates the new keys**/
                    data _null_;
                    %if %qupcase(&sysscp)=WIN %then %do; 
                        file "&workdir.\_newsurv_emf94.sasxreg";
                    %end;
                    %else %do;
                        file "&workdir./_newsurv_emf94.sasxreg";
                    %end;
                    put '[CORE\PRINTING\PRINTERS\EMF\ADVANCED]';
                    put '"Description"="Enhanced Metafile Format"';
                    put '"Metafile Type"="EMF"';
                    put '"Vector Alpha"=int:0';
                    put '"Image 32"=int:1';
                    run;    
                    %if %qupcase(&sysscp)=WIN %then %do; 
                        proc registry export="&workdir.\_newsurv_preexisting.sasxreg";/* Exports current SASUSER Keys */
                        proc registry import="&workdir.\_newsurv_emf94.sasxreg"; /* Import the new keys */
                        run;
                    %end;
                    %else %do;
                        proc registry export="&workdir./_newsurv_preexisting.sasxreg";/* Exports current SASUSER Keys */
                        proc registry import="&workdir./_newsurv_emf94.sasxreg"; /* Import the new keys */
                        run;
                    %end;
                %end;
                %else %do;
                    ods graphics / imagefmt=&plottype;  
                %end;
            %end;
            %else %if %qupcase(&plottype)=TIFF or %qupcase(&plottype)=TIF %then %do;
                ods graphics / imagefmt=png;    
            %end;
            %else %do;
                ods graphics / imagefmt=&plottype;  
            %end;          
        %end;
        %if %sysevalf(%superq(plotname)^=,boolean) %then %do; 
            ods graphics / reset=index imagename="&plotname";
        %end;  
        /**Turns on Scalable-Vector-Graphics**/
        %if &svg = 1 %then %do;
            %if %qupcase(&destination) = RTF or %qupcase(&destination) = EXCEL or %qupcase(&destination) = POWERPOINT %then %do;
                ods graphics / OUTPUTFMT=EMF;
            %end;
            %else %if %qupcase(&destination) = HTML %then %do;
                ods graphics / OUTPUTFMT=SVG;
            %end;
            %else %do;
                ods graphics / OUTPUTFMT=STATIC;
            %end;
        %end;
            
        proc template;
            %do i = 1 %to %sysfunc(countw(%superq(_destinations),|));
                define style styles.newsurv_axes&i;
                    parent=styles.%scan(%superq(_styles),&i,|);
                    class GraphAxisLines /
                        ContrastColor=&axiscolor
                        color=&axiscolor;
                    class GraphBorderLines /
                        contrastcolor=&axiscolor
                        color=&axiscolor;
                    class GraphWalls /
                        contrastcolor=&axiscolor
                        color=&axiscolor;
                    End;
            %end;
            %if &_listing^=1 and %sysevalf(%superq(gpath)^=,boolean) %then %do;
                define style styles.newsurv_axes_list;
                    parent=styles.listing;
                    class GraphAxisLines /
                        ContrastColor=&axiscolor
                        color=&axiscolor;
                    class GraphBorderLines /
                        contrastcolor=&axiscolor
                        color=&axiscolor;
                    class GraphWalls /
                        contrastcolor=&axiscolor
                        color=&axiscolor;
                    End;
            %end;
                
        run;
        
        %do i = 1 %to %sysfunc(countw(%superq(_destinations),|));
            ods %scan(%superq(_destinations),&i,|) style=newsurv_axes&i
                %if %sysevalf(%qupcase(%scan(%superq(_destinations),&i,|))=LISTING,boolean) and %sysevalf(%superq(gpath)^=,boolean) %then %do;
                    gpath="&gpath"
                %end;;
        %end;
        
        /**Save image to specified location**/
        %if %sysevalf(%superq(gpath)^=,boolean) and &_listing^=1 %then %do;
            ods listing gpath="&gpath" style=newsurv_axes_list;
        %end;
        /**Sets plot options**/
        ods graphics /  antialias antialiasmax=&antialiasmax scale=off width=&width height=&height 
            %if &sysvlong >= 9.04.01M5P091317 %then %do; LINEPATTERNOBSMAX=10000000 %end;;
        /**Generates the Plot**/
        options notes;
        proc sgrender data=_plot template=_km;
        run;
        options nonotes;
        /**Changes Potential Registry Changes back**/
        %if %qupcase(&plottype)=EMF or (&svg=1 and %qupcase(&destination)=RTF and %qupcase(&plottype)^=TIFF and %qupcase(&plottype)^=TIF)
            or (&svg=1 and %qupcase(&destination)=EXCEL)
            or (&svg=1 and %qupcase(&destination)=POWERPOINT) %then %do;
            %if &transparent=0 and &_any_trans=0 %then %do;
                proc registry clearsasuser; /* Deletes the SASUSER directory */
                proc registry import="&workdir./_newsurv_preexisting.sasxreg";/* Imports starting SASUSER Keys */
                run;
            %end;
        %end;
        /**Creates the TIFF file from the PNG file created earlier**/
        %else %if %qupcase(&plottype)=TIFF or %qupcase(&plottype)=TIF %then %do;
            %local _fncheck _fncheck2;
            options nonotes;
            %if %sysevalf(%superq(gpath)=,boolean) %then %do;
                filename nsurvpng "./&plotname..png"; 
                filename nsurvtif "./&plotname..tiff";
                data _null_;
                    x=fexist('nsurvpng');
                    x2=fdelete('nsurvtif');
                    call symput('_fncheck',strip(put(x,12.)));
                    call symput('_fncheck2',strip(put(x2,12.)));
                run;
                %if %sysevalf(%superq(_fncheck)^=1,boolean) %then %do;
                    filename nsurvpng "./&plotname.1.png"; 
                %end;
            %end;
            %else %do;
                filename nsurvpng "%sysfunc(tranwrd(&gpath./&plotname..png,//,/))"; 
                filename nsurvtif "%sysfunc(tranwrd(&gpath./&plotname..tiff,//,/))"; 
                data _null_;
                    x=fexist('nsurvpng');
                    x2=fdelete('nsurvtif');
                    call symput('_fncheck',strip(put(x,12.)));
                    call symput('_fncheck2',strip(put(x2,12.)));
                run;
                %if %sysevalf(%superq(_fncheck)^=1,boolean) %then %do;
                    filename nsurvpng "%sysfunc(tranwrd(&gpath./&plotname.1.png,//,/))"; 
                %end;
            %end;
            options notes;
            goptions device=&tiffdevice gsfname=nsurvtif 
                xmax=&width ymax=&height 
                xpixels=%sysevalf(%sysfunc(compress(&width,abcdefghijklmnopqrstuvwxyz,i))*&dpi) 
                ypixels=%sysevalf(%sysfunc(compress(&height,abcdefghijklmnopqrstuvwxyz,i))*&dpi)
                imagestyle=fit iback=nsurvpng
                %if &border=1 %then %do; border %end;
                %else %if &border=0 %then %do; noborder %end;;
            proc gslide;
            run;
            quit; 
            data _null_;
                x=fdelete('nsurvpng');
            run;
            filename nsurvpng clear;
            filename nsurvtif clear;
        %end;
        options nonotes; 
        
        
    %end;
    
    /**Print out summary table**/
    %if &summary=1 %then %do;
        %local _multmethodcheck _multmethodlist _any_by;
        %do i = 1 %to %sysfunc(countw(%superq(_destinations),|));
            %if %sysevalf(%qupcase(%qscan(%superq(_destinations),&i,|))=EXCEL,boolean) %then %do;
                ods excel options(sheet_name='NEWSURV' frozen_headers="2") style=_newsurvtable;
            %end;
            %else %do;
                ods %scan(%superq(_destinations),&i,|) style=_newsurvtable;
            %end;
        %end;
        /**Check p-values for footnote purposes**/
        proc sql noprint;
            %if &tablemergepval=1 %then %do;
                update %if %sysevalf(%superq(out)=,boolean)=0 %then %do; &out %end;
                       %else %do; _summary %end;
                    set pval=strip(coalescec(covpval,pval_inter))
                    where cmiss(covpval,pval_inter)^=2 and missing(pval);
                update %if %sysevalf(%superq(out)=,boolean)=0 %then %do; &out %end;
                %else %do; _summary %end;
                    set pvalmv=strip(coalescec(covpvalmv,pval_intermv))
                    where cmiss(covpvalmv,pval_intermv)^=2 and missing(pvalmv);
            %end;
            select max(^missing(by_level)) into :_any_by separated by '' from 
                %if %sysevalf(%superq(out)=,boolean)=0 %then %do; &out %end;
                %else %do; _summary %end;;
            select count(distinct modeltype) into :_multmethodcheck
                from 
                    %if %sysevalf(%superq(out)=,boolean)=0 %then %do; &out %end;
                    %else %do; _summary %end;;
            select distinct modeltype into :_multmethodlist separated by '|'
                from 
                    %if %sysevalf(%superq(out)=,boolean)=0 %then %do; &out %end;
                    %else %do; _summary %end;;
            /**Print Summary with PROC REPORT**/
            /**Determine columns to be showin in report**/                    
            /*Model Classifications*/
            %local _tndisplay_model modelcolumns;
            %let _tndisplay_model=0;
            %let modelcolumns=TITLE|FOOTNOTE;            
            /*Statistics*/
            %local _tndisplay_stat statcolumns _display_current _med_check _hr_check _tl_check _tlmv_check _pval_check _pvalmv_check _pval_inter_check _pval_intermv_check _covpval_check _covpvalmv_check;
            %let _pval_check=0;
            %let _pvalmv_check=0;
            %let _covpval_check=0;
            %let _covpvalmv_check=0;
            %let _tndisplay_stat=0;        
            %let statcolumns=TOTAL|EVENT|MEDIAN|TIMELIST|EV_N|N_EV|TOTALMV|EVENTMV|EV_NMV|N_EVMV|HR|MEDIANMV|TIMELISTMV|HRMV|PVAL|PVALMV|COVPVAL|COVPVALMV|CINDEX|CINDEXMV|RMST|RMTL|RMST_PVAL|PVAL_INTER|PVAL_INTERMV;          
            /*Statistics*/
            /**Take only first entry if repeated entries are listed**/
            %let _tabledisplay=%qupcase(%scan(&tabledisplay,1,%str( )));
            %do i = 2 %to %sysfunc(countw(&tabledisplay,%str( )));
                %let _display_current=%qupcase(%scan(%superq(tabledisplay),&i,%str( )));
                %let _test=0;
                %do j = 1 %to %sysevalf(&i-1);
                    %if &_display_current=%qupcase(%scan(%superq(_tabledisplay),&j,|)) %then %let _test=1;
                %end;
                %if ^&_test %then %let _tabledisplay=&_tabledisplay|&_display_current;
            %end;           
            %do i = 1 %to %sysfunc(countw(&_tabledisplay,|));            
                %let _display_current=%qupcase(%scan(%superq(_tabledisplay),&i,|));
                %let _test=0;
                select ifn(count(*)>0,1,0) into :_test 
                    from %if %sysevalf(%superq(out)=,boolean)=0 %then %do; &out %end;
                         %else %do; _summary %end;
                    where ^missing(&_display_current);
                %if &_test=1 %then %do j = 1 %to %sysfunc(countw(&modelcolumns,|));
                    %if &_display_current=%scan(&modelcolumns,&j,|) %then %do;
                        %let _tndisplay_model=%sysevalf(&_tndisplay_model+1);
                        %let _tndisplay_model_&_tndisplay_model=&_display_current;
                    %end;
                %end;
                %if &_test=1 %then %do j = 1 %to %sysfunc(countw(&statcolumns,|));
                    %if &_display_current=%scan(&statcolumns,&j,|) %then %do;
                        %let _tndisplay_stat=%sysevalf(&_tndisplay_stat+1);
                        %let _tndisplay_stat_&_tndisplay_stat=&_display_current;
                        %if %sysevalf(%qupcase(&_display_current)=MEDIAN,boolean) %then %let _med_check=1;
                        %else %if %sysevalf(%qupcase(&_display_current)=TIMELIST,boolean) %then %let _tl_check=1;
                        %else %if %sysevalf(%qupcase(&_display_current)=PVAL,boolean) %then %let _pval_check=1;
                        %else %if %sysevalf(%qupcase(&_display_current)=COVPVAL,boolean) %then %let _covpval_check=1;
                        %else %if %sysevalf(%qupcase(&_display_current)=PVALMV,boolean) %then %let _pvalmv_check=1;
                        %else %if %sysevalf(%qupcase(&_display_current)=COVPVALMV,boolean) %then %let _covpvalmv_check=1;
                        %else %if %sysevalf(%qupcase(&_display_current)=PVAL_INTER,boolean) %then %let _pval_inter_check=1;
                        %else %if %sysevalf(%qupcase(&_display_current)=PVAL_INTERMV,boolean) %then %let _pval_intermv_check=1;
                        %else %if %sysfunc(find(&_display_current,HR,i))>0 %then %let _hr_check=1;
                    %end;
                %end;
            %end;
        quit;
        data _report;
            set %if %sysevalf(%superq(out)=,boolean)=0 %then %do; &out %end;
                %else %do; _summary %end; end=last;
             
             array hlist (4) $20. _temporary_;
             array dlist (&_tndisplay_stat) $50. (%do i = 1 %to &_tndisplay_stat; %if &i>1 %then %do; , %end; "%qupcase(&&_tndisplay_stat_&i)" %end;);
             retain dlist;
             if _n_=1 then do;
                do i = 1 to dim(dlist);
                    if dlist(i) in("TIMELIST","MEDIAN") then _univ_=1;
                    if dlist(i) in("TIMELISTMV","MEDIANMV") then _multiv_=1;
                    retain _univ_ _multiv_;
                    drop _univ_ _multiv_;
                    if dlist(i) in("TIMELIST","MEDIAN","TIMELISTMV","MEDIANMV") and &_multmethodcheck=1 then do j = 1 to dim(hlist);
                        if missing(hlist(j)) or hlist(j) = 
                            ifc(&_multmethodlist=0,'KM',
                                ifc(&_multmethodlist=1,'1-KM',
                                    ifc(&_multmethodlist=2,'CIF',
                                        ifc(&_multmethodlist=3,'INVWTS',
                                            ifc(&_multmethodlist=4,'DIRECT',''))))) then do;
                            hlist(j)=ifc(&_multmethodlist=0,'KM',
                                        ifc(&_multmethodlist=1,'1-KM',  
                                            ifc(&_multmethodlist=2,'CIF',
                                                ifc(&_multmethodlist=3,'INVWTS',
                                                    ifc(&_multmethodlist=4,'DIRECT','')))));
                            j=dim(hlist);
                        end;
                    end;
                    else if dlist(i) in("HR","HRMV") then do j = 1 to dim(hlist);
                        if missing(hlist(j)) or hlist(j) = 'Cox' then do;
                            hlist(j)='Cox';
                            j=dim(hlist);
                        end;
                    end;
                end;

                nfoot=0;
                if cmiss(of hlist(*))<dim(hlist) then do i=1 to dim(hlist)-cmiss(of hlist(*));
                    nfoot+1;
                    if hlist(i)='KM' then call symputx('foot'||strip(put(nfoot,12.0)),"Kaplan-Meier method;",'l');
                    else if hlist(i)='1-KM' then call symputx('foot'||strip(put(nfoot,12.0)),"1-Kaplan-Meier method;",'l');
                    else if hlist(i)='CIF' then call symputx('foot'||strip(put(nfoot,12.0)),"Cumulative incidence method;",'l');
                    else if hlist(i)='INVWTS' then do;
                        if _univ_=1 then do;
                            call symputx('foot'||strip(put(nfoot,12.0)),"Kaplan-Meier method;",'l');
                        end;
                        if _multiv_=1 then do;
                            if _univ_=1 then nfoot+1;
                            call symputx('foot'||strip(put(nfoot,12.0)),"Inverse weights adjusting method;",'l');
                        end;
                    end;
                    else if hlist(i)='DIRECT' then do;
                        if _univ_=1 then do;
                            call symputx('foot'||strip(put(nfoot,12.0)),"Kaplan-Meier method;",'l');
                        end;
                        if _multiv_=1 then do;
                            if _univ_=1 then nfoot+1;
                            call symputx('foot'||strip(put(nfoot,12.0)),"Direct adjusting method;",'l');
                        end;
                    end;
                    else if hlist(i)='Cox' then call symputx('foot'||strip(put(nfoot,12.0)),"Cox model;",'l');
                end;
            end;
            array plist (6) $20. _temporary_;
            do i = 1 to dim(dlist);
                if dlist(i) in('PVAL' 'PVALMV' 'COVPVAL' 'COVPVALMV' 'PVAL_INTER' 'PVAL_INTERMV') then do;
                    if ^missing(compress(vvaluex(dlist(i)),'- ')) then do j = 1 to dim(plist);
                           if missing(plist(j)) or compress(scan(vvaluex(dlist(i)),2,'^'),' {super}') = plist(j) then do;
                               if missing(plist(j)) then plist(j)=compress(scan(vvaluex(dlist(i)),2,'^'),' {super}');
                               if dlist(i)='PVAL' then pval=strip(scan(vvaluex(dlist(i)),1,'^'))||'^{super '||strip(put(j+nfoot,12.0))||'}';
                               else if dlist(i)='PVALMV' then pvalmv=strip(scan(vvaluex(dlist(i)),1,'^'))||'^{super '||strip(put(j+nfoot,12.0))||'}';
                               else if dlist(i)='COVPVAL' then covpval=strip(scan(vvaluex(dlist(i)),1,'^'))||'^{super '||strip(put(j+nfoot,12.0))||'}';
                               else if dlist(i)='COVPVALMV' then covpvalmv=strip(scan(vvaluex(dlist(i)),1,'^'))||'^{super '||strip(put(j+nfoot,12.0))||'}';
                               else if dlist(i)='PVAL_INTER' then pval_inter=strip(scan(vvaluex(dlist(i)),1,'^'))||'^{super '||strip(put(j+nfoot,12.0))||'}';
                               else if dlist(i)='PVAL_INTERMV' then pval_intermv=strip(scan(vvaluex(dlist(i)),1,'^'))||'^{super '||strip(put(j+nfoot,12.0))||'}';
                               j=dim(plist);
                           end;
                   end;
                end;
            end;
            if last then do;
                if cmiss(of plist(*))<dim(plist) then do i=1 to dim(plist)-cmiss(of plist(*));
                    nfoot+1;
                    if plist(i)='#' then call symputx('foot'||strip(put(nfoot,12.0)),"Likelihood-ratio test;",'l');
                    else if plist(i)='$' then call symputx('foot'||strip(put(nfoot,12.0)),"Score test;",'l');
                    else if plist(i)='*' then call symputx('foot'||strip(put(nfoot,12.0)),"Logrank test;",'l');
                    else if plist(i)='**' then call symputx('foot'||strip(put(nfoot,12.0)),"One-sided logrank test;",'l');
                    else if plist(i)='@' then call symputx('foot'||strip(put(nfoot,12.0)),"Wilcoxon test;",'l');
                    else if plist(i)='+' then call symputx('foot'||strip(put(nfoot,12.0)),"Wald Chi-Square test;",'l');
                    else if plist(i)='G' then call symputx('foot'||strip(put(nfoot,12.0)),"Gray's k-sample test for equality of cumulative incidence functions;",'l');
                end;
                call symputx('nfoot',strip(put(nfoot,12.0)),'l');
            end;
            keep modelnum 
               %do i = 1 %to &_tndisplay_model;
                   %superq(_tndisplay_model_&i)
               %end;
               modeltype by_label by_order by_level classlevel subind subtitle 
               %do i = 1 %to &_tndisplay_stat;
                   %superq(_tndisplay_stat_&i)
               %end;;
        run;
        proc sql noprint;
            %local _ppt _other _destinations _styles k;
            select max(ifn(upcase(destination) ^in('LISTING' 'OUTPUT' 'POWERPOINT'),1,0)),
                max(ifn(upcase(destination) in('POWERPOINT'),1,0))
                into :_other separated by '',:_ppt separated by '' from sashelp.vdest;
        quit;
        %if &_listing = 1 %then %do;
            %do i = 1 %to %sysfunc(countw(%superq(_destinations),|));
                ods %scan(%superq(_destinations),&i,|) select none;
            %end;
            ods listing;
            ods listing select all;
            data _out_listing;
                set _report;
                array _chars_ (*) $2000. _character_;   
                retain _tl;
                do i = 1 to dim(_chars_);
                    if upcase(vname(_chars_(i))) in('PVAL' 'PVALMV' 'COVPVAL' 'COVPVALMV' 'PVAL_INTER' 'PVAL_INTER') and _chars_(i)^='--' then do;
                        if ^missing(_chars_(i)) then _chars_(i)=strip(scan(_chars_(i),1,'^'))||repeat('*',input(compress(scan(_chars_(i),2,'^'),'{super }'),12.)-1);
                    end;
                    else if upcase(vname(_chars_(i)))='SUBTITLE' then _chars_(i)='A0A0A0'x||strip(_chars_(i));
                end;  

                if subind=0 and %sysevalf(%superq(title)^=,boolean) and %sysfunc(find(&tabledisplay,title,i))>0 then subtitle=repeat('A0A0'x,1)||strip(subtitle); 
                else if %sysevalf(%superq(title)^=,boolean) and %sysfunc(find(&tabledisplay,title,i))>0 then subtitle=repeat('A0A0'x,subind)||strip(subtitle);
                else if subind>0 then subtitle=repeat('A0A0'x,subind-1)||strip(subtitle);
                
            run;
            proc contents data=_out_listing noprint out=_outldict;
            run;
            
            proc sql noprint;
                %local _list_cvars;
                select upcase(name) into :_list_cvars separated by '|' from _outldict where type=2 and upcase(name) ^in('TITLE' 'FOOTNOTE' 'BY_LABEL' 'BY_LEVEL');
                %do i = 1 %to %sysfunc(countw(%superq(_list_cvars),|,m));
                    %local _list_%scan(%superq(_list_cvars),&i,|,m);
                %end;
                select %do i = 1 %to %sysfunc(countw(%superq(_list_cvars),|,m)); 
                          %if &i>1 %then %do; , %end;
                          max(length(strip(%scan(%superq(_list_cvars),&i,|,m))))+4
                       %end;
                       into %do i = 1 %to %sysfunc(countw(%superq(_list_cvars),|,m)); 
                                %if &i>1 %then %do; , %end;
                                :_list_%scan(%superq(_list_cvars),&i,|,m) separated by ''
                            %end;
                       from _out_listing;
                       
                %local _list_totlength headerlength datalength;
                %let _list_totlength=0;
                %do i = 1 %to %sysfunc(countw(%superq(_list_cvars),|,m));
                    %if %sysevalf(%qupcase(%scan(%superq(_list_cvars),&i,|,m))=SUBTITLE,boolean) %then %let headerlength=0;
                    %else %do;
                        %if %sysfunc(countw(%superq(t%scan(%superq(_list_cvars),&i,|,m)header),~,m))>1 %then %do j=1 %to %sysfunc(countw(%superq(t%scan(%superq(_list_cvars),&i,|,m)header),~,m));
                            %if &j=1 %then %let headerlength=%length(%qscan(%superq(t%scan(%superq(_list_cvars),&i,|,m)header),&j,~,m));
                            %else %let headerlength=%sysfunc(max(&headerlength,%length(%qscan(%superq(t%scan(%superq(_list_cvars),&i,|,m)header),&j,~,m))));
                        %end;
                        %else %let headerlength=%length(%superq(t%scan(%superq(_list_cvars),&i,|,m)header));
                    %end;
                    %let datalength=%superq(_list_%scan(%superq(_list_cvars),&i,|,m));
                    %let _list_%scan(%superq(_list_cvars),&i,|,m)=%sysfunc(max(&headerlength,&datalength));
                    %let _list_totlength=%sysevalf(&_list_totlength+%sysfunc(max(&headerlength,&datalength)));
                %end;
                alter table _out_listing
                    modify %do i = 1 %to %sysfunc(countw(%superq(_list_cvars),|,m)); 
                                %if &i>1 %then %do; , %end;
                                %scan(%superq(_list_cvars),&i,|,m) char(%superq(_list_%scan(%superq(_list_cvars),&i,|,m)))
                            %end;;  
            quit;
            options linesize=%sysfunc(max(64,%sysfunc(min(256,&_list_totlength)))) nocenter notes; 
            proc report data=_out_listing spacing=0 nowd split='~' missing;
                
                columns
                    ("%sysfunc(tranwrd(%superq(tabletitle),`,~))~%sysfunc(repeat(-,&_list_totlength-1))"
                    modelnum /**Used for sorting and distinguishing Models**/
                        /**Titles and Footnotes are listed first to be used in compute blocks later**/    
                        %do i = 1 %to &_tndisplay_model;
                            %superq(_tndisplay_model_&i)
                        %end;
                        modeltype /**Used to determine if KM or CIF**/
                        by_order by_level  classlevel subind subtitle /*These are always shown*/
                        /*Statistics*/ 
                        %do i = 1 %to &_tndisplay_stat;
                            %superq(_tndisplay_stat_&i)
                        %end;);
                
                define subind / display noprint; /**Not Printed but defined**/
                define modelnum / order noprint;/**Used to keep models in order**/
                define modeltype / order noprint;/**Used to keep models in order**/
                define by_order / order noprint; /**Not Printed but defined**/
                define by_level / order noprint; /**Not Printed but defined**/
                define classlevel / order noprint; /**Not Printed but defined**/
                
                %do i = 1 %to &_tndisplay_model;
                    define %superq(_tndisplay_model_&i) / order noprint;/**Used in compute blocks later**/
                    compute 
                        %if %sysevalf(%qupcase(%superq(_tndisplay_model_&i))=TITLE,boolean) %then %do;
                            before 
                        %end;
                        %else %do;
                            after 
                        %end;
                        %superq(_tndisplay_model_&i) / style={just=l};
                        line @4 %superq(_tndisplay_model_&i) $2000.;
                    endcomp;
                %end;    
                compute before modelnum;
                    count+1;
                    x="%sysfunc(repeat(-,&_list_totlength-1))";
                    if count=1 then len=0;
                    else len=length(x);
                    line @1 x $varying. len;
                endcomp;
                     
                /**If KM and CIF are both used, then add subtitle to each model**/
                %if &_multmethodcheck>1 %then %do;
                    compute before modeltype / style={just=l};
                        length text $150.;
                        if modeltype=0 then text='Kaplan-Meier methods';
                        else if modeltype=1 then text='(1-Kaplan-Meier) methods';
                        else if modeltype=2 then text='Cumulative incidence methods';
                        else if modeltype=3 then text='Inverse weights adjusting methods';
                        else if modeltype=4 then text='Direct adjusting methods';
                        line @4 text $150.;
                    endcomp;
                %end;
                /**Widths are set to 30 to avoid throwing line-size errors**/
                /**This Summary Table is not designed to be viewed in the output window**/
                define subtitle / display "~%sysfunc(repeat(-,&_list_subtitle-1))"  id; /**Class level descriptions**/
                %do i =1 %to &_tndisplay_stat;
                    %if %qupcase(&&_tndisplay_stat_&i)=TIMELIST %then %do;
                        define timelist / display  center
                            %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                %if %sysevalf(%superq(foot&j)=%nrstr(Kaplan-Meier method;),boolean) %then %do;
                                    "&ttimelistheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_timelist-1))"
                                    %let j=&nfoot;
                                %end;
                                %else %if %sysevalf(%superq(foot&j)=%nrstr(1-Kaplan-Meier method;),boolean) %then %do;
                                    "&ttimelistheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_timelist-1))"
                                    %let j=&nfoot;
                                %end;
                                %else %if %sysevalf(%superq(foot&j)=%nrstr(Cumulative incidence method;),boolean) %then %do;
                                    "&ttimelistheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_timelist-1))"
                                    %let j=&nfoot;
                                %end;
                                %else %if &j=&nfoot %then %do;
                                    "&ttimelistheader~%sysfunc(repeat(-,&_list_timelist-1))"
                                %end;
                            %end;
                            %else %do;
                                "&ttimelistheader~%sysfunc(repeat(-,&_list_timelist-1))"
                            %end;;
                    %end;
                    %else %if %qupcase(&&_tndisplay_stat_&i)=TIMELISTMV %then %do;
                        define timelistmv / display  center
                            %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                %if %sysevalf(%superq(foot&j)=%nrstr(Direct adjusting method;),boolean) %then %do;
                                    "&ttimelistmvheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_timelistmv-1))"
                                    %let j=&nfoot;
                                %end;
                                %else %if %sysevalf(%superq(foot&j)=%nrstr(Inverse weights adjusting method;),boolean) %then %do;
                                    "&ttimelistmvheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_timelistmv-1))"
                                    %let j=&nfoot;
                                %end;
                                %else %if &j=&nfoot %then %do;
                                    "&ttimelistmvheader~%sysfunc(repeat(-,&_list_timelistmv-1))"
                                %end;
                            %end;
                            %else %do;
                                "&ttimelistmvheader~%sysfunc(repeat(-,&_list_timelistmv-1))"
                            %end;;
                    %end;
                    %else %if %qupcase(&&_tndisplay_stat_&i)=HR %then %do;
                        define hr / display  
                            %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                %if %sysevalf(%superq(foot&j)=%nrstr(Cox model;),boolean) %then %do;
                                    "&thrheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_hr-1))"
                                %end;
                            %end; center;
                    %end;
                    %else %if %qupcase(&&_tndisplay_stat_&i)=HRMV %then %do;
                        define hrmv / display  
                            %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                %if %sysevalf(%superq(foot&j)=%nrstr(Cox model;),boolean) %then %do;
                                    "&thrmvheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_hrmv-1))"
                                %end;
                            %end; center;
                    %end;
                    %else %if %qupcase(&&_tndisplay_stat_&i)=MEDIAN %then %do;
                        define median / display center                      
                            %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                %if %sysevalf(%superq(foot&j)=%nrstr(Kaplan-Meier method;),boolean) %then %do;
                                    "&tmedianheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_median-1))"
                                    %let j=&nfoot;
                                %end;
                                %else %if %sysevalf(%superq(foot&j)=%nrstr(1-Kaplan-Meier method;),boolean) %then %do;
                                    "&tmedianheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_median-1))"
                                    %let j=&nfoot;
                                %end;
                                %else %if %sysevalf(%superq(foot&j)=%nrstr(Cumulative incidence method;),boolean) %then %do;
                                    "&tmedianheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_median-1))"
                                    %let j=&nfoot;
                                %end;
                                %else %if &j=&nfoot %then %do;
                                    "&tmedianheader~%sysfunc(repeat(-,&_list_median-1))"
                                %end;
                            %end;
                            %else %do;
                                "&tmedianheader~%sysfunc(repeat(-,&_list_median-1))"
                            %end;;
                    %end;
                    %else %if %qupcase(&&_tndisplay_stat_&i)=MEDIANMV %then %do;
                        define medianmv / display center                      
                            %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                %if %sysevalf(%superq(foot&j)=%nrstr(Direct adjusting method;),boolean) %then %do;
                                    "&tmedianmvheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_medianmv-1))"
                                    %let j=&nfoot;
                                %end;
                                %else %if %sysevalf(%superq(foot&j)=%nrstr(Inverse weights adjusting method;),boolean) %then %do;
                                    "&tmedianmvheader%sysfunc(repeat(*,&j-1))~%sysfunc(repeat(-,&_list_medianmv-1))"
                                    %let j=&nfoot;
                                %end;
                                %else %if &j=&nfoot %then %do;
                                    "&tmedianmvheader~%sysfunc(repeat(-,&_list_medianmv-1))"
                                %end;
                            %end;
                            %else %do;
                                "&tmedianmvheader~%sysfunc(repeat(-,&_list_medianmv-1))"
                            %end;;
                    %end;
                    %else %do;
                        define %qupcase(&&_tndisplay_stat_&i) / display 
                        "%superq(%sysfunc(compress(t&&_tndisplay_stat_&i..header)))~%sysfunc(repeat(-,%superq(%qupcase(_list_%superq(_tndisplay_stat_&i)))-1))" center;
                    %end;
                %end;
                /**Choose the furthest right column in the dataset**/
                /**This allows other columns to be manipulated as they are all to the left of this column**/
                /**PROC REPORT does not give values to columns to the right of the currently processed columns**/
                /**Creates the overall footnotes at the bottom of the table**/
                compute after /style={just=l};
                    line @1 "%sysfunc(repeat(-,&_list_totlength-1))";
                    /**Creates footnotes with symbols based on which columns are requested with TABLEDISPLAY**/
                    %if &nfoot > 0 %then %do;
                        line @4 
                            %do i = 1 %to &nfoot;
                                "%sysfunc(repeat(*,&i-1))&&foot&i "
                            %end;;
                    %end;
                    %if &nfoot=0 and %sysevalf(%superq(tablefootnote)=,boolean) %then %do;
                        line @4 " ";
                    %end;
                    /**Lists the table footnote**/
                    %if %sysevalf(%superq(tablefootnote)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(tablefootnote),`,m));
                        line @4 "%qscan(%superq(tablefootnote),&i,`,m)";
                    %end;
                endcomp;
            run;
            options nonotes &_center;
            proc datasets nolist nodetails;
                %if &debug=0 %then %do;
                    delete _outldict _out_listing ;
                %end;
            quit;
            ods select all;
        %end;
        %if &_other = 1 or &_ppt = 1 or &_rtf=1 or %sysevalf(%superq(outdoc)^=,boolean) %then %do;
            ODS LISTING CLOSE;
            %local _rloop;
            %if &_other=1 %then %let _rloop=OTHER;
            %if &_ppt=1 and %sysevalf(%superq(_rloop)=,boolean) %then %let _rloop=PPT;
            %else %if &_ppt=1 %then %let _rloop=&_rloop|PPT;
            %do rloop = 1 %to %sysfunc(countw(&_rloop,|));
                %if %sysevalf(%scan(%superq(_rloop),&rloop,|)=OTHER,boolean) %then %do;
                    %if &_ppt=1 %then %do;
                        ods POWERPOINT exclude all;
                    %end;
                %end;
                %else %do;
                    ods POWERPOINT select all;
                    ods POWERPOINT  style=_newsurvtableppt;
                    %do k = 1 %to %sysfunc(countw(%superq(_destinations),|));
                        %if %sysevalf(%qupcase(%qscan(%superq(_destinations),&k,|))^=POWERPOINT,boolean) %then %do;
                            ods %scan(%superq(_destinations),&k,|) exclude all;
                        %end;
                    %end;
                %end;
                options notes;
                proc report data=_report 
                    nowd split='~' missing;
                    
                    columns
                        (modelnum /**Used for sorting and distinguishing Models**/
                            /**Titles and Footnotes are listed first to be used in compute blocks later**/    
                            %do i = 1 %to &_tndisplay_model;
                                %superq(_tndisplay_model_&i)
                            %end;
                            modeltype /**Used to determine if KM or CIF**/
                            by_order by_level classlevel subind subtitle /*These are always shown*/
                            /*Statistics*/ 
                            %do i = 1 %to &_tndisplay_stat;
                                %superq(_tndisplay_stat_&i)
                            %end;);
                    
                    define subind / order noprint; /**Not Printed but defined**/
                    define by_order / order noprint; /**Not Printed but defined**/
                    define by_level / order noprint; /**Not Printed but defined**/
                   
                    define classlevel / order noprint; /**Not Printed but defined**/
                    define modelnum / order noprint;/**Used to keep models in order**/
                    define modeltype / order noprint;/**Used to keep models in order**/
                    
                    %do i = 1 %to &_tndisplay_model;
                        define %superq(_tndisplay_model_&i) / order noprint;/**Used in compute blocks later**/
                        compute 
                            %if %sysevalf(%qupcase(%superq(_tndisplay_model_&i))=TITLE,boolean) %then %do;
                                before 
                            %end;
                            %else %do;
                                after 
                            %end;
                            %superq(_tndisplay_model_&i) / 
                            %if %sysevalf(%qupcase(%superq(_tndisplay_model_&i))=TITLE,boolean) %then %do;
                                style={fontweight=bold just=l bordertopstyle=solid bordertopcolor=black bordertopwidth=0.1}
                            %end;
                            %else %if %sysfunc(find(&tabledisplay,title,i))=0 %then %do;
                                style={fontweight=bold just=l borderbottomstyle=solid borderbottomcolor=black borderbottomwidth=0.1}
                            %end;
                            %else %do;
                                style={fontweight=bold just=l}
                            %end;;
                            line %superq(_tndisplay_model_&i) $2000.;
                        endcomp;
                    %end;       
                    /**If KM and CIF are both used, then add subtitle to each model**/
                    %if &_multmethodcheck>1 %then %do;
                        compute before modeltype / 
                            %if %sysfunc(find(&tabledisplay,title,i))=0 and
                                %sysfunc(find(&tabledisplay,footnote,i))=0 %then %do;
                                style={fontweight=bold just=l bordertopstyle=solid bordertopcolor=black bordertopwidth=0.1}
                            %end;
                            %else %do;
                                style={fontweight=bold just=l} 
                            %end;;
                            length text $150.;
                            if modeltype=0 then text='Kaplan-Meier methods';
                            else if modeltype=1 then text='(1-Kaplan-Meier) methods';
                            else if modeltype=2 then text='Cumulative incidence methods';
                            else if modeltype=3 then text='Inverse weights adjusting methods';
                            else if modeltype=4 then text='Direct adjusting methods';
                            line @1 text $150.;
                        endcomp;
                    %end;
                    /**Widths are set to 30 to avoid throwing line-size errors**/
                    /**This Summary Table is not designed to be viewed in the output window**/
                    define subtitle / order ''  id style={cellwidth=&tsubtitlewidth}; /**Class level descriptions**/
                    %do i =1 %to &_tndisplay_stat;
                        %if %qupcase(&&_tndisplay_stat_&i)=TIMELIST %then %do;
                            define timelist / display  
                                %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                    %if %sysevalf(%superq(foot&j)=%nrstr(Kaplan-Meier method;),boolean) %then %do;
                                        "&ttimelistheader^{super &j}"
                                    %end;
                                    %else %if %sysevalf(%superq(foot&j)=%nrstr(1-Kaplan-Meier method;),boolean) %then %do;
                                        "&ttimelistheader^{super &j}"
                                    %end;
                                    %else %if %sysevalf(%superq(foot&j)=%nrstr(Cumulative incidence method;),boolean) %then %do;
                                        "&ttimelistheader^{super &j}"
                                    %end;
                                %end;
                                %else %do;
                                    "&ttimelistheader"
                                %end;
                            style={cellwidth=&ttimelistwidth};
                        %end;
                        %else %if %qupcase(&&_tndisplay_stat_&i)=TIMELISTMV %then %do;
                            define timelistmv / display  
                                %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                    %if %sysevalf(%superq(foot&j)=%nrstr(Direct adjusting method;),boolean) %then %do;
                                        "&ttimelistmvheader^{super &j}"
                                    %end;
                                    %else %if %sysevalf(%superq(foot&j)=%nrstr(Inverse weights adjusting method;),boolean) %then %do;
                                        "&ttimelistmvheader^{super &j}"
                                    %end;
                                %end;
                                %else %do;
                                    "&ttimelistmvheader"
                                %end;
                            style={cellwidth=&ttimelistmvwidth};
                        %end;
                        %else %if %qupcase(&&_tndisplay_stat_&i)=HR %then %do;
                            define hr / display 
                                %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                    %if %sysevalf(%superq(foot&j)=%nrstr(Cox model;),boolean) %then %do;
                                        "&thrheader^{super &j}"
                                    %end;
                                %end; center style={cellwidth=&thrwidth};
                        %end;
                        %else %if %qupcase(&&_tndisplay_stat_&i)=HRMV %then %do;
                            define hrmv / display   
                                %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                    %if %sysevalf(%superq(foot&j)=%nrstr(Cox model;),boolean) %then %do;
                                        "&thrmvheader^{super &j}"
                                    %end;
                                %end; center
                            style={cellwidth=&thrmvwidth};
                        %end;
                        %else %if %qupcase(&&_tndisplay_stat_&i)=MEDIAN %then %do;
                            define median / display                          
                                %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                    %if %sysevalf(%superq(foot&j)=%nrstr(Kaplan-Meier method;),boolean) %then %do;
                                        "&tmedianheader^{super &j}"
                                    %end;
                                    %else %if %sysevalf(%superq(foot&j)=%nrstr(1-Kaplan-Meier method;),boolean) %then %do;
                                        "&tmedianheader^{super &j}"
                                    %end;
                                    %else %if %sysevalf(%superq(foot&j)=%nrstr(Cumulative incidence method;),boolean) %then %do;
                                        "&tmedianheader^{super &j}"
                                    %end;
                                %end;
                                %else %do;
                                    "&tmedianheader"
                                %end;
                            center  style={cellwidth=&tmedianwidth};
                        %end;
                        %else %if %qupcase(&&_tndisplay_stat_&i)=MEDIANMV %then %do;
                            define medianmv / display                          
                                %if &nfoot>0 %then %do j = 1 %to &nfoot;
                                    %if %sysevalf(%superq(foot&j)=%nrstr(Direct adjusting method;),boolean) %then %do;
                                        "&tmedianmvheader^{super &j}"
                                    %end;
                                    %else %if %sysevalf(%superq(foot&j)=%nrstr(Inverse weights adjusting method;),boolean) %then %do;
                                        "&tmedianmvheader^{super &j}"
                                    %end;
                                %end;
                                %else %do;
                                    "&tmedianmvheader"
                                %end;
                            center  style={cellwidth=&tmedianmvwidth};
                        %end;
                        %else %do;
                            define %qupcase(&&_tndisplay_stat_&i) / display  "%superq(%sysfunc(compress(t&&_tndisplay_stat_&i..header)))" center 
                             style={cellwidth=%superq(%sysfunc(compress(t%superq(_tndisplay_stat_&i)width)))};
                        %end;
                    %end;
                    /*Print title before table*/
                    compute before _page_/ 
                        style={leftmargin=0.06in bordertopstyle=none borderbottomstyle=solid borderbottomwidth=0.1 borderbottomcolor=black
                               vjust=bottom just=left color=black background=white};
                        %do i = 1 %to %sysfunc(max(1,%sysfunc(countw(%superq(tabletitle),`,m))));
                            line @1 "%scan(%superq(tabletitle),&i,`,m)";
                        %end;
                    endcomp;
                    /**Choose the furthest right column in the dataset**/
                    /**This allows other columns to be manipulated as they are all to the left of this column**/
                    /**PROC REPORT does not give values to columns to the right of the currently processed columns**/
                    compute &&_tndisplay_stat_&_tndisplay_stat;
                        %if &tableshading=1 %then %do;
                            /**Creates alternating-shading using modulo arithmatic**/
                            if classlevel=0 then do;
                                call define('subtitle','style/merge','style={fontweight=bold}');
                                shade=1;
                            end;
                            else shade+1;
                            if mod(shade,2)=0 then call define(_row_, 'style/merge','style={background=GREYEF');
                        %end;
                        %do i =1 %to &_tndisplay_stat;
                            %if %qupcase(&&_tndisplay_stat_&i)=TIMELIST %then %do;
                                if find(timelist,':')=0 then call define('timelist','style/merge','style={just=c}');
                            %end;
                        %end;
                        %do i =1 %to &_tndisplay_stat;
                            %if %qupcase(&&_tndisplay_stat_&i)=TIMELISTMV %then %do;
                                if find(timelistmv,':')=0 then call define('timelistmv','style/merge','style={just=c}');
                            %end;
                        %end;
                        /**Creates an indented list of class levels using the subind variable**/
                        if subind=1 then call define('subtitle','style/merge','style={indent=0.12in}');
                        else if subind=2 then call define('subtitle','style/merge','style={indent=0.18in}');
                        else if subind=3 then call define('subtitle','style/merge','style={indent=0.24in}');
                        else if subind=4 then call define('subtitle','style/merge','style={indent=0.3in}');

                        %if %sysfunc(find(&tabledisplay,title,i))=0 and
                            %sysfunc(find(&tabledisplay,footnote,i))=0 and &_multmethodcheck=1 %then %do;
                            if classlevel=0 then call define(_row_,'style/merge','style={bordertopstyle=solid bordertopcolor=black bordertopwidth=0.1}');
                        %end;
                    endcomp;
                    /**Creates the overall footnotes at the bottom of the table**/
                    compute after /style={leftmargin=0.06in bordertopstyle=solid bordertopwidth=0.1 bordertopcolor=black vjust=top color=black just=l};
                        /**Creates footnotes with symbols based on which columns are requested with TABLEDISPLAY**/
                        %if &nfoot > 0 %then %do;
                            line @1 
                                %do i = 1 %to &nfoot;
                                    "^{super &i}&&foot&i "
                                %end;;
                        %end;
                        %if &nfoot=0 and %sysevalf(%superq(tablefootnote)=,boolean) %then %do;
                            line @1 " ";
                        %end;
                        /**Lists the table footnote**/
                        %if %sysevalf(%superq(tablefootnote)=,boolean) =0 %then %do i = 1 %to %sysfunc(countw(%superq(tablefootnote),`,m));
                            line @1 "%qscan(%superq(tablefootnote),&i,`,m)";
                        %end;
                    endcomp;
                run;
                options nonotes; 
            %end;
        %end;
    %end;   
    
    %if %sysevalf(%qupcase(%superq(axiscolor))^=BLACK,boolean) or &summary=1 %then %do;
        %do i = 1 %to %sysfunc(countw(%superq(_destinations),|));
            ods %scan(%superq(_destinations),&i,|) style=%scan(%superq(_styles),&i,|);
        %end;
    %end;
        
    /**Closes the ODS file**/
    %if %sysevalf(%superq(outdoc)=,boolean)=0 %then %do;
        ods &destination close;
    %end;
    /**Outputs Plot Dataset**/
    %if %sysevalf(%superq(outp)=,boolean)=0 %then %do;
        data &outp;
            set _plot;
        run;
    %end; 
        
    %errhandl:
          
    /**Reset Graphics Options**/
    ods graphics / reset=all;
    %if &_listing=1 %then %do;
        ods Listing;
    %end;
    %else %do;
        ods listing close;
    %end;
    ods select all;
    ods results;   
    /**Delete temporary datasets**/
    proc datasets nolist nodetails;
        %if &debug=0 %then %do;
            delete _temp _options _plot _report
                %if %sysevalf(%superq(out)=,boolean) %then %do;
                    _summary
                %end; 
                %do z = 1 %to &nmodels;
                    _plot_&z
                %end; ;
        %end;
    quit;  
    /**Reload previous Options**/ 
    ods path &_odspath;
    options mergenoby=&_mergenoby &_notes &_qlm linesize=&_linesize msglevel=&_msglevel &_mprint;
    goptions device=&_device gsfname=&_gsfname &_gborder
        xmax=&_xmax ymax=&_ymax xpixels=&_xpixels ypixels=&_ypixels imagestyle=&_imagestyle iback=&_iback;
    %put NEWSURV has finished processing, runtime: %sysfunc(putn(%sysevalf(%sysfunc(TIME())-&_starttime.),mmss8.4)); 
    %mend;

