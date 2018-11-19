#Thorben Quast, 14th November 2018
#output: HGCal-DWC alignment for each muon run

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptFit()
from copy import deepcopy

import root_numpy as rn


#I am using a target based workflow system (https://github.com/spotify/luigi)
from workflows.tbAnalysisOctober2018.tasks.reconstruction import ProduceNtuples #the ProduceNtuples object ultimately returns me the path to the ntuple file
from workflows.tbAnalysisOctober2018.configurations.runConfiguration import modulesUnderTest

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--run', help='run to analyse', required=True)
parser.add_argument('--NLayers', help='number of layers', type=int, required=False, default=40)
parser.add_argument("--output_file", type=str, help="output file path to the alignment txt file", default="0")
args = parser.parse_args()

run = args.run
outfile_root_path = args.output_file.replace(".txt", ".root")
outfile_txt_path = args.output_file
#usage: python MIPAlignment.py --run <run number> --NLayers 40 --output_file <path to the txt file with the DWC-HGCAL alignment files>

#implementation of ROOT fits of gaussian distributions to a histogram incl. presetting of fit parameters
def gausFit(histogram, gaus):
	gaus.SetRange(histogram.GetXaxis().GetXmin(), histogram.GetXaxis().GetXmax())

	if histogram.GetEntries()==0:
		return

	max_pos = histogram.GetXaxis().GetBinCenter(histogram.GetMaximumBin())
	gaus.SetRange(max_pos-1.0, max_pos+1.0)
	gaus.SetParameter(0, histogram.GetMaximum())
	gaus.SetParLimits(0, histogram.GetMaximum()*0.6, histogram.GetMaximum()*1.6)
	gaus.SetParameter(1, max_pos)
	gaus.SetParameter(2, histogram.GetStdDev())

  	histogram.Fit(gaus, "RQ")




#loading the ntuples
tree_rechits = ROOT.TChain('rechitntupler/hits','rechits')
tree_tracks = ROOT.TChain('trackimpactntupler/impactPoints','tracks')
print "Adding Run", run
tree_rechits.Add(ProduceNtuples(runNumber=run).output().path)	
tree_tracks.Add(ProduceNtuples(runNumber=run).output().path)
tree_rechits.AddFriend(tree_tracks)

outFile = ROOT.TFile(outfile_root_path, "RECREATE")
outFile.Close()
beam_shifts = []
#iterate over all layers
for LAYER in range(1, 1+args.NLayers):
	dir_name = "layer_%s"%LAYER
	
	outFile = ROOT.TFile(outfile_root_path, "UPDATE")
	outFile.mkdir(dir_name)
	outFile.Close()

	tree_rechits.SetAlias("ref_x", "impactX_HGCal_layer_%s" % (LAYER))
	tree_rechits.SetAlias("ref_y", "impactY_HGCal_layer_%s" % (LAYER))
	tree_rechits.SetAlias("dX", "rechit_x+ref_x")	#mind the sign in the coordinates, rechit-x ~ -1 x hit-x
	tree_rechits.SetAlias("dY", "rechit_y+ref_y")   #mind the sign in the coordinates, rechit-y ~ -1 x hit-y

	baseselection = "(dwcReferenceType>=13)&&(trackChi2_X<150)&&(trackChi2_Y<150)&&(rechit_layer==%s)"%LAYER


	gaus_x = ROOT.TF1("gaus_layer%s"%LAYER, "gaus")
	gaus_y = ROOT.TF1("gaus_layer%s"%LAYER, "gaus")
	
	#plots the difference of hit-x to dwc-x for this layer
	h1_dx = ROOT.TH1F("h1_dx_layer%s"%LAYER, "h1_dx_layer%s"%LAYER, 100, -5., 5.)
	tree_rechits.Project("h1_dx_layer%s"%LAYER, "dX", "(run==%s)&&!((rechit_channel==22)&&(rechit_chip==3))&&(rechit_amplitudeHigh>30)&&%s"%(run, baseselection))
	#alignment is translational only, shift by the mean difference of hit-x vs. dwc-x
	gausFit(h1_dx, gaus_x)
	mu_x = gaus_x.GetParameter(1) if (gaus_x.GetChisquare()/gaus_x.GetNDF()<100 or LAYER==1) else beam_shifts[LAYER-1][1] 
	
	#analogous for y
	h1_dy = ROOT.TH1F("h1_dy_layer%s"%LAYER, "h1_dy_layer%s"%LAYER, 100, -5., 5.)
	tree_rechits.Project("h1_dy_layer%s"%LAYER, "dY", "(run==%s)&&!((rechit_channel==22)&&(rechit_chip==3))&&(rechit_amplitudeHigh>30)&&%s"%(run, baseselection))
	gausFit(h1_dy, gaus_y)
	mu_y = gaus_y.GetParameter(1) if (gaus_y.GetChisquare()/gaus_y.GetNDF()<100 or LAYER==1) else beam_shifts[LAYER-1][2]
	
	h1_dx.SetTitle("Shift: %s"% mu_x)
	h1_dy.SetTitle("Shift: %s"% mu_y)
	beam_shifts.append([LAYER, mu_x, mu_y])
	print "Run",args.run,",   Layer: ",LAYER,"-->",beam_shifts[LAYER-1]

	#save the histograms for later usage
	outFile = ROOT.TFile(outfile_root_path, "UPDATE")
	outFile.cd(dir_name)
	h1_dx.Write()
	h1_dy.Write()
	outFile.Close()

#write out the translational alignment parameters to a txt file
with open(outfile_txt_path, "a") as alignmentfile:
	import tabulate	#you might have to install that package in order to run this script
	alignmentfile.write(tabulate.tabulate(beam_shifts, headers=["layer", "shift x", "shift y"]))