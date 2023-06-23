# tempBTNAD
code use to analyze data and generate figures 

Figure 1

*Toxicity screen analysis
Main script:
- drug_screen_data_analysis.m
Assist functions:
- venn.m
Data input:
- drugscreen_normOD.mat

*Non-antibiotic toxicity screen analysis
Main script:
- drugScreen176analysis.m
Assist functions:
- natsort.m
Data input:
- 20210603sortedAntihostDrugs.xlsx
- OD_ecoli_20210605_plate1.xlsx
- OD_ecoli_20210605_plate2.xlsx

*Chemical similarity
Main script: 
- similarityAnalysis.m
Data input: 
- drugNamesForQuery_simmat_upd100223_cid.xlsx
- strainResultsAggreared.mat
- drugMetadata.xlsx

FIGURE 2&3

*For networks and heatmap:
Main script: 
- networkAnalysis.m
Assist functions: 
- checkPwEnrichment.m
- extractPwAnnFromGMT.m
- findHitStrains.m
- getNodeInsidePolygon.m
- plotSingleHCluster.m
- plotSingleNetwork.m
Data input: 
- strainResultsAggreared.mat
- drugMetadata.xlsx
- strainMetadata.csv
- GO bp Entrez GMT.csv
- GO cc Entrez GMT.csv
- GO mf Entrez GMT.csv
- KEGG Entrez GMT.csv
- Ecocyc Pathways Entrez GMT.csv
- COG function Entrez GMT.csv
- prophageStrains.mat
- pumpStrains.mat

*For GO enrichment analysis:
Main script: 
- GOenrichmentPlotting.m
Assist functions: 
- extractPwAnnFromGMT.m
Data input: 
- htClusters.mat
- abClusters.mat
- RevigoOutputAb.xlsx
- RevigoOutputHt.xlsx

*For pump analysis:
Main script: 
- pumpsAnalysis.m
Assist functions: 
- calcEmpiricalPval.m
Data input: 
- strainResultsAggreared.mat
- drugMetadata.xlsx
- strainMetadata.csv
- prophageStrains.mat
- pumpStrains.mat

FIGURE 4

* Dose response curves
Main script: 
- plotIC50curvesEvolvedLines.m
Data input: 
- CrossResistanceABX.mat

* Cross-resistance bars
Main script: 
- plotCrossresistanceBar.m
Data input: 
- EvoLinesIC50Results.mat
- EvoLinesIC50Results2.mat

* non-antibiotic cross-resistance
Main script: 
- plotIC50nonAb.m
Data input: 
- PYR_TCBZ_LAMRoomTempIC50.mat
