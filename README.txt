Workflow for R analysis for XLMS data from XlinkX, Proteome Discoverer
1.	From XlinkX, export txt file for CSMs and Crosslinks.
2.	R script to read txt files, for each acquisition mode and database,
1.	Read replicate Crosslinks txt
2.	Add identifying replicate number
3.	Create dataframe from all replicates -> Output replicates txt
4.	Rearrange Cross-link positions so that Position A> Position B (pmax) (pmin)
5.	Filter out non-BSA cross-links (grepl) -> Output BSA crosslinks only txt
6.	Filter out BSA cross-links (!grepl) -> Output non-BSA crosslinks txt
7.	Keep non duplicates (!duplicated) -> Output BSA crosslinks no duplicates txt
8.	Convert BSA crosslinks no duplicates data into dataframe for xiNET export -> Output xiNET csv file with BSA crosslinks no duplicates
9.	For each individual replicate, filter out XlinkX score based on score cut-off
10.	For each individual replicate, rearrange cross-link positions so that Position A> Position B (pmax)(pmin)
11.	For each individual replicate, keep non duplicates (!duplicated). 
12.	Create new dataframe for all non-duplicated crosslinks for each replicate.
13.	For each individual replicate, filter out non-BSA cross-links (grepl)
14.	Create new dataframe for non-duplicated BSA crosslinks for each replicate.
15.	For each individual replicate, determine non-BSA crosslinks by subtracting BSA crosslinks from all crosslinks
16.	Create dataframe for number of total crosslinks, BSA crosslinks and non-BSA crosslinks.
17.	Create dataframe for crosslink numbers for all 4 acquisition modes, both simple and complex sample and both targeted or general database
18.	Output venn diagram for overlap in crosslinks between replicates, with duplicated crosslinks removed. (draw.triple.venn)
	i.	Area1= replicate 1 crosslink number
	ii.	Area2= replicate 2 crosslink number
	iii.	Area3= replicate 3 crosslink number
	iv.	Area12=replicate1+replicate2 overlap
	v.	Area23=replicate2+replicate3 overlap
	vi.	Area13=replicate1+replicate3 overlap
	vii.	Area123=total crosslinks – (Area1+Area2+Area3-Area12-Area23-Area13)

Files:
xlBSA_replicatesanalysis #No filter
xlBSA_replicateanalysis #1% Error filter applied


Author: Zheng Ser
Kentsis Lab