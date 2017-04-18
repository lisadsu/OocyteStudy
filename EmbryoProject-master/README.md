# EmbryoProject
This repository contains all the code used to gather and analyze data and generate figures for the embryo project.

It contains the following:

1. Labview code (Automated Pressure Control v20.vi) used to operate the system hardware (perform micropipette aspiration and record video). Requires drivers for Firgelli linear actuator (included in repo), as well as Thorlabs DCC1545M USB camera.

2. Matlab code to analyze the videos of micropipette aspiration, perform (mostly) automated detection of aspiration depth over time and calculate embryo mechanical parameters. Briefly:

a. To *load/save* experimental viability / ground truth data: /OocyteStudies/loadDataOocytes.m (for new oocyte measurements), loadDataClinicalStudy.m (for clinical study human embryos), saveNewMorphology.m (for old oocyte/embryo/human research data)

b. To *extract* mechanical parameters and save data files: extractParamsClinicalStudy.m (for clinical study human embryos), depthDetectionAuto.m (for mouse oocytes, embryos and generic other measurements). There is also a GUI (measGUI.m) to measure these parameters which will automatically detect the pressure used for measurement (does not work on clinical system because it requires an accompanying pressure file).

c. To *plot and visualize* mechanical data: /OocyteStudies/oocyteFertAnalysis.m (for new oocyte measurements), plotParamsClinicalStudy.m (for clinical study human embryos), plotAllMouseEmbryosSoFar.m / plotAllHumanSoFar.m / plotAllMouseOocytesSoFar.m (for old oocyte/embryo/human research embryo measurements)

d. To analyze cortical granule images for Nat Comm paper: CGprocessingPipeline.m, analyzeMicroinjectionImages.m

e. To perform further analysis: analyzeMicroinjectionData.m, makeMicroinjectionFigs.m, plotEmbryoToOocyteChange.m

f. To make figures for papers: makePaperFigs.m (for NatComm paper), makeReviewPaperFigs.m (for review paper in MHR), /OocyteStudies/makeOocytePaperFigs.m (for unpublished oocyte data). Figs for clinical study paper are in plotParamsClinicalStudy.m.

g. To evaluate classification performance and make predictions: Classification/classifyEmbryos.m, update the data to be loaded in either loadSVMDataClinical.m, loadSVMdataOocyte.m, loadSVMdata.m, etc. The output of loadSVMdata should include a list of class labels, a matrix of parameter values for all observations, a list of indices (either cross-validation or prediction) and a vector of embryo label numbers (if making predictions).

3. Shell scripts to run various aspects of the RNA-seq analysis pipeline (FastX Clipper/Trimmer, STAR aligner, HTSeq-Count, SamTools, RNASeQC) and prepare data for processing in R. There are also R scripts to run edgeR as well as ComBat, and generate figures for RNA-seq data. 


