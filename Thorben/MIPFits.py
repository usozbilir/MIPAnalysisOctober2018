#Thorben Quast, 19 November 2018
#MIP Spectrum fitting for a fixed module

import ROOT
ROOT.gStyle.SetOptFit()
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("/afs/cern.ch/user/t/tquast/CMS_HGCal_Upgrade/styles/CLICStyle_C.so")
ROOT.CLICStyle()
from copy import deepcopy
import os



import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--module', help='module to analyse', required=True, type=int)
parser.add_argument("--input_file", type=str, help="input file path", default="0")
parser.add_argument("--output_file", type=str, help="output file path", default="0")
args = parser.parse_args()
#usage: python MIPFits.py --module <moduleNumber> --input_file <path to the root file with the content as defined in MIPHistograms.py> --output_file <output .txt file with the stored fit results>
MODULE = args.module
outputFilePath = args.output_file
outputFilePath_root = outputFilePath.replace(".txt", ".root")


#MIP signal fit model
def gausConvLandau(x, par):
	#implementation adapted from: https://root.cern.ch/root/html/tutorials/fit/langaus.C.html
	invsq2pi = 0.3989422804014
	mpshift = -0.22278298

	np = 100
	sc = 5.

	mpc = par[1]
	xlow = x[0] - sc * par[3]
	xupp = x[0] + sc * par[3]

	step = (xupp-xlow) / np

	_sum = 0.

	for i in range(0, np):
		xx = xlow + (i-0.5) * step
		flandau = ROOT.TMath.Landau(xx, mpc, par[0], True)
		_sum += flandau * ROOT.TMath.Gaus(x[0]-xx, 0, par[3])

	return (par[2] * _sum)


output_dir = "/".join(outputFilePath.split("/")[0:-1])
if not os.path.exists(output_dir):
	os.mkdir(output_dir)

h1_energy_spectra_all = {}
h1_energy_spectra_hit = {}

outfile = ROOT.TFile(outputFilePath_root, "RECREATE")
outfile.Close()

#efficiency estimate: nominator divided by denominator
infile = ROOT.TFile(args.input_file, "READ")
h2_occupancy = deepcopy(infile.Get("module_%s/occupancy_module%s" % (MODULE, MODULE)))				
h2_occupancy_hit = deepcopy(infile.Get("module_%s/occupancy_hit_module%s" % (MODULE, MODULE)))	
infile.Close()
h2_occupancy_hit.Divide(h2_occupancy)
h2_occupancy_hit.SetStats(False) 
h2_occupancy_hit.GetZaxis().SetRangeUser(0.05, 1.)
h2_occupancy_hit.SetTitle("MIP efficiency estimate for layer: %s"%MODULE)
h2_occupancy_hit.GetXaxis().SetTitle("Impact X_{DWC} [cm]")
h2_occupancy_hit.GetYaxis().SetTitle("Impact Y_{DWC} [cm]")
canvas = ROOT.TCanvas("c", "c", 1000, 800)
h2_occupancy_hit.Draw("COLZ")
canvas.Print("%s/efficiency_estimate_module%s.pdf"%(output_dir, MODULE))


#load full and cleaned spectra for each channel on this module individually
print "Performing the fits:"
for chip in range(4):
	for ch in range(64):
		if ch%2:
			continue
		key = chip*100+ch

		infile = ROOT.TFile(args.input_file, "READ")
		h1_energy_spectra_all[key] = deepcopy(infile.Get("module_%s/module%s_chip%s_ch%s_all" % (MODULE, MODULE, chip, ch)))
		h1_energy_spectra_hit[key] = deepcopy(infile.Get("module_%s/module%s_chip%s_ch%s_signal" % (MODULE, MODULE, chip, ch)))
		infile.Close()

		#Rebin to a reasonable number of bins
		h1_energy_spectra_hit[key].RebinX(1)
		h1_energy_spectra_all[key].RebinX(1)

		Nentries_signal = h1_energy_spectra_hit[key].GetEntries()	
		FitChi2 = -1
		FitNDF = -1
		MPV = -1
		LandauMPV = -1
		LandauMPV_err = -1
		Gaussian_Width = -1

		gausConvLandau_func=None
		if Nentries_signal > 200:
			print "Fitting for chip:",chip," channel: ",ch
			#define model and set constraints
			gausConvLandau_func = ROOT.TF1("gausConvLandau_func",gausConvLandau, 20., 120., 4)
			gausConvLandau_func.SetParName(0, "Landau width")
			gausConvLandau_func.SetParName(1, "Landau MPV")
			gausConvLandau_func.SetParName(2, "Const")
			gausConvLandau_func.SetParName(3, "Gaussian width")
			gausConvLandau_func.SetParameter(0, 7.)
			gausConvLandau_func.SetParameter(1, 50.)
			gausConvLandau_func.SetParameter(2, 1.)
			gausConvLandau_func.SetParameter(3, 5.)
			gausConvLandau_func.SetParameter(2, h1_energy_spectra_hit[key].Integral())
			gausConvLandau_func.SetParLimits(0, 1., 10.)
			gausConvLandau_func.SetParLimits(1, 0., 70.)
			gausConvLandau_func.SetParLimits(3, 1., 20.)
			gausConvLandau_func.SetParameter(2, h1_energy_spectra_hit[key].Integral())
			fit_status = h1_energy_spectra_hit[key].Fit(gausConvLandau_func, "0S")
		
			reject_fit = False
			if (gausConvLandau_func.GetNDF()==0) or (gausConvLandau_func.GetChisquare()/gausConvLandau_func.GetNDF()) > 5.:
				reject_fit = True
			if (gausConvLandau_func.GetParameter(0)<=1.1) or (gausConvLandau_func.GetParameter(0)>=9.9):
				reject_fit = True
			if (gausConvLandau_func.GetParameter(1)<=0.1) or (gausConvLandau_func.GetParameter(1)>=69.9):
				reject_fit = True
			if (gausConvLandau_func.GetParameter(3)<=1.1) or (gausConvLandau_func.GetParameter(3)>=19.1):
				reject_fit = True
			
			if not reject_fit:
				FitChi2 = gausConvLandau_func.GetChisquare()
				FitNDF = gausConvLandau_func.GetNDF()
				MPV = gausConvLandau_func.GetMaximumX()
				LandauMPV = gausConvLandau_func.GetParameter(1)
				LandauMPV_err = gausConvLandau_func.GetParError(1)
				Gaussian_Width = gausConvLandau_func.GetParameter(3)
						
		
		canvas_MIPsFitted = ROOT.TCanvas("canvas_MIPs", "canvas_MIPs", 2000, 1000)
		h1_energy_spectra_all[key].SetLineColor(ROOT.kBlack)
		h1_energy_spectra_hit[key].SetLineColor(ROOT.kBlue+1)
		h1_energy_spectra_hit[key].GetXaxis().SetTitle("Reconstructed amplitude HG [ADC]")
		h1_energy_spectra_hit[key].GetYaxis().SetTitle("N_{events}")
		h1_energy_spectra_hit[key].Draw()
		h1_energy_spectra_all[key].Draw("SAME")
		if gausConvLandau_func!=None:
			gausConvLandau_func.Draw("SAME")
		h1_energy_spectra_hit[key].GetYaxis().SetRangeUser(1, h1_energy_spectra_all[key].GetMaximum())
		if gausConvLandau_func!=None:
			legend_c1 = ROOT.TLegend(0.10, 0.80, 0.60, 0.90)
			legend_c1.AddEntry(gausConvLandau_func, "Maximum: %s [HG ADC] "%round(gausConvLandau_func.GetMaximumX(), 2))
			legend_c1.Draw()
		canvas_MIPsFitted.SetLogy()
		canvas_MIPsFitted.Print("%s/module%s_ski%s_ch%s.pdf" % (output_dir, MODULE, chip, ch))


		outfile = ROOT.TFile(outputFilePath_root, "UPDATE")
		h1_energy_spectra_hit[key].Write()
		h1_energy_spectra_all[key].Write()
		outfile.Close()

		#write out the fit results to a txt file
		with open(outputFilePath, "a") as outfile:
			outfile.write(str(MODULE)+", ")
			outfile.write(str(chip)+", ")
			outfile.write(str(ch)+", ")
			outfile.write(str(int(Nentries_signal))+", ")
			outfile.write(str(FitChi2)+", ")
			outfile.write(str(FitNDF)+", ")
			outfile.write(str(MPV)+", ")
			outfile.write(str(LandauMPV)+", ")
			outfile.write(str(LandauMPV_err)+", ")
			outfile.write(str(Gaussian_Width)+", ")
			outfile.write("\n")