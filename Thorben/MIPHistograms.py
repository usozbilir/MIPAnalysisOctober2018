#Thorben Quast, 16th November 2018
#output: 1. spectra with and without DWC tracking
#2. beam occupancies with and without hit requirement for efficiency study

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptFit()
from copy import deepcopy

#root_numpy loads data from a ROOT TTree and stores it conveniently into a root numpy array
import root_numpy as rn
import numpy as np
import pdb

#CONFIGURATION
max_DWC_track_chi2=150
energy_quantity = "rechit_amplitudeHigh"
nBins_spectra = 120
spectra_min = 0
spectra_max = 120
max_dX = 0.5	#cm
max_dY = 0.5	#cm

max_dX_efficiency = 1.2	#cm
max_dY_efficiency = 1.2	#cm



from workflows.tbAnalysisOctober2018.tasks.reconstruction import ProduceNtuples
from workflows.tbAnalysisOctober2018.configurations.runConfiguration import modulesUnderTest

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--run', help='run to analyse', type=int, required=True)
parser.add_argument('--NLayers', help='number of layers', type=int, required=True)
parser.add_argument("--alignment_file", type=str, help="file path to alignment information", default="0", required=True)
parser.add_argument("--output_file", type=str, help="output file path", default="0", required=True)
args = parser.parse_args()
#usage: python MIPHistograms.py --run <run number> --alignment_file <path to the alignment file from MIPAlignment.py> --NLayers 40 --output_file <output file with the stored energy spectra>

run = args.run
alignmentfile_path = args.alignment_file
outfile_path = args.output_file

#load the translational alignment parameters
print "Reading alignment file"
beam_shifts = {}
align_layer, align_shiftx, align_shifty = np.genfromtxt(open(args.alignment_file, "r"), usecols=(0, 1, 2), unpack=True, skip_header=2)
for i, layer in enumerate(align_layer):
	beam_shifts[layer] = [align_shiftx[i], align_shifty[i]]


tree_rechits = ROOT.TChain('rechitntupler/hits','rechits')
tree_tracks = ROOT.TChain('trackimpactntupler/impactPoints','tracks')
print "Adding Run", run
tree_rechits.Add(ProduceNtuples(runNumber=run).output().path)	
tree_tracks.Add(ProduceNtuples(runNumber=run).output().path)


outFile = ROOT.TFile(outfile_path, "RECREATE")
outFile.Close()

#dictionaries for storage of all TH*F
energy_spectra_all = {}
energy_spectra_hit = {}
occupancy_reco = {}
h2_occupancy = {}
h2_occupancy_hit = {}

#loading the ntuples using root_numpy
print "loading the data"
data_rechits = rn.tree2array(tree_rechits, branches=["rechit_x", "rechit_y", "rechit_module", "rechit_layer", "rechit_chip", "rechit_channel", energy_quantity])
data_tracks = rn.tree2array(tree_tracks, branches=["dwcReferenceType", "trackChi2_X", "trackChi2_Y"]+["impactX_HGCal_layer_%s"%LAYER for LAYER in range(1, args.NLayers+1)]+["impactY_HGCal_layer_%s"%LAYER for LAYER in range(1, args.NLayers+1)])

#later we need to know the module-->layer mapping, create it based on the ntuples
print "making module to layer map"
module_to_layer = {}
for ev in data_rechits:
	for hit_index in range(len(ev["rechit_module"])):
		module = ev["rechit_module"][hit_index]
		layer = ev["rechit_layer"][hit_index]
		if not module in module_to_layer:
			module_to_layer[module] = layer


#initialising the spectra and occupancy graphics
print "Initialising the MIP distributions"
for MODULE in modulesUnderTest:
	h2_occupancy[MODULE] = ROOT.TH2F("occupancy_module%s"%MODULE, "occupancy_module%s"%MODULE, 140, -7, 7, 140, -7, 7)
	h2_occupancy_hit[MODULE] = ROOT.TH2F("occupancy_hit_module%s"%MODULE, "occupancy_hit_module%s"%MODULE, 140, -7, 7, 140, -7, 7)	

	for chip in range(4):
		for ch in range(64):
			if ch % 2:
				continue
			key = 1000*MODULE + 100*chip + ch
			energy_spectra_hit[key] = ROOT.TH1F("module%s_chip%s_ch%s_signal" % (MODULE, chip, ch), "module%s_chip%s_ch%s_signal" % (MODULE, chip, ch), nBins_spectra, spectra_min, spectra_max)
			energy_spectra_all[key] = ROOT.TH1F("module%s_chip%s_ch%s_all" % (MODULE, chip, ch), "module%s_chip%s_ch%s_all" % (MODULE, chip, ch), nBins_spectra, spectra_min, spectra_max)
			occupancy_reco[key] = ROOT.TH2F("module%s_chip%s_ch%s_occupancy_reco" % (MODULE, chip, ch), "module%s_chip%s_ch%s_occupancy" % (MODULE, chip, ch), 140, -7, 7, 140, -7, 7)



print "Filling the energy spectra"
NEvents = len(data_tracks)
#event loop
for event in range(NEvents):
	if not event%100:
		print "[MIP Spectra] Run:",args.run,"-Event",event,"/",NEvents
	
	#hit loop
	for hit_entry in range(len(data_rechits[event]["rechit_module"])):
		MODULE = data_rechits[event]["rechit_module"][hit_entry]
		LAYER = data_rechits[event]["rechit_layer"][hit_entry]
		hit_x = data_rechits[event]["rechit_x"][hit_entry]
		hit_x = data_rechits[event]["rechit_y"][hit_entry]
		chip = data_rechits[event]["rechit_chip"][hit_entry]
		ch = data_rechits[event]["rechit_channel"][hit_entry]
		hit_amplitude = data_rechits[event][energy_quantity][hit_entry]

		key = 1000*MODULE + 100*chip + ch
		#energy spectra are filled without any requirement
		energy_spectra_all[key].Fill(hit_amplitude)

		#dwc track quality requirements
		impactX = data_tracks[event]["impactX_HGCal_layer_%s"%LAYER]
		impactY = data_tracks[event]["impactY_HGCal_layer_%s"%LAYER]
		if impactX == -999.:	#require: DWC synchronised data
			continue
		if data_tracks[event]["dwcReferenceType"] < 13.:	#require: either all or all but DWC C contributing
			continue
		if data_tracks[event]["trackChi2_X"] > max_DWC_track_chi2:	#require: good track in x
			continue
		if data_tracks[event]["trackChi2_Y"] > max_DWC_track_chi2:	#require: good track in y
			continue
		occupancy_reco[key].Fill(impactX, impactY)

		#DWC-tracking selection: require hit within a certain distance to the DWC track
		if abs(data_rechits[event]["rechit_x"][hit_entry]+impactX-beam_shifts[LAYER][0])>max_dX:
			continue
		if abs(data_rechits[event]["rechit_y"][hit_entry]+impactY-beam_shifts[LAYER][1])>max_dY:
			continue

		energy_spectra_hit[key].Fill(hit_amplitude)
	#end if hit loop

	#efficiency estimation part
	if data_tracks[event]["impactX_HGCal_layer_%s"%1] == -999.:
		continue
	if data_tracks[event]["dwcReferenceType"] < 13.:
		continue
	if data_tracks[event]["trackChi2_X"] > max_DWC_track_chi2:
		continue
	if data_tracks[event]["trackChi2_Y"] > max_DWC_track_chi2:
		continue

	for MODULE in modulesUnderTest:
		if not MODULE in module_to_layer:	
			continue
		LAYER = module_to_layer[MODULE]
		impactX = data_tracks[event]["impactX_HGCal_layer_%s"%LAYER] 
		impactY = data_tracks[event]["impactY_HGCal_layer_%s"%LAYER] 

		h2_occupancy[MODULE].Fill(impactX, impactY)

	efficiency_increased_modules = []
	for hit_entry in range(len(data_rechits[event]["rechit_module"])):
		MODULE = data_rechits[event]["rechit_module"][hit_entry]
		LAYER = module_to_layer[MODULE]
		if MODULE in efficiency_increased_modules:		#no double or more counting for the efficiency denominator
			continue
		if abs(data_rechits[event]["rechit_x"][hit_entry]+impactX-beam_shifts[LAYER][0])>max_dX_efficiency:
			continue
		if abs(data_rechits[event]["rechit_y"][hit_entry]+impactY-beam_shifts[LAYER][1])>max_dY_efficiency:
			continue
		efficiency_increased_modules.append(MODULE)
		h2_occupancy_hit[MODULE].Fill(impactX, impactY)

#end of event loop

#write everything to file
print "Writing to the file..."
for MODULE in modulesUnderTest:
	dir_name = "module_%s"%MODULE
	outFile = ROOT.TFile(outfile_path, "UPDATE")
	outFile.mkdir(dir_name)
	outFile.cd(dir_name)
	h2_occupancy[MODULE].Write()
	h2_occupancy_hit[MODULE].Write()
	for chip in range(4):
		for ch in range(64):
			if ch % 2:
				continue
			key = 1000*MODULE + 100*chip + ch
			if key in energy_spectra_hit:
				energy_spectra_hit[key].Write()
				energy_spectra_all[key].Write()
				occupancy_reco[key].Write()
	outFile.Close()	

