# TPS_Substrate_Sensing
Supporting Information for "Differential Substrate Sensing in Terpene Synthases from Plants and Microorganisms. Insights from Structural, Bioinformatic, and EnzyDock Analyses"

- MC Integration code for active site   
	mc_intergal/scripts/0README explains

- The complete multiple sequence analysis   
	MSA/clustalo-I20240311-103911-0791-80447995-p1m-aln.fasta

- Modeller script   
	protein_prep/model-single-ligand.py

- CHARMM soft potential script   
	protein_prep/prep_prot_2a_soft.inp

- Starting structures for EnzyDock   
	pdb_for_enzydock/...

- Carbocation parameters   
	initial pdb + topology + parameter + patch files for ligands used for docking with EnzyDock.   
	cc_ligands/…

- PDB files for selected poses from EnzyDock runs   
	all_mm_best/…

- CSV file with energetic data from docking   
	all_mm_best2_to_si/all_poses.csv

- Electrostatic CHARMM script   
	charmm_analyses/get_potential_1n21.inp

- Ligand-protein interaction CHARMM script   
	charmm_analyses/get_list_4kux_1.inp
