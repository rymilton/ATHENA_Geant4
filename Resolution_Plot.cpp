#include <iostream>
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "TAxis.h"

void Resolution_Plot() {

    Double_t g4_beamE[7] = {1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0};
	Double_t g4_beamE_piplus[14] = {1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 40., 50., 60., 70., 80., 90., 100.};
	Double_t g4_res_elec_5deg[7] = {0.0889507, 0.0633381, 0.0514269, 0.0405578, 0.0289683, 0.0206064, 0.0170256};
	Double_t g4_res_piplus_5deg[14] = {0.25023, 0.213926, 0.191676, 0.161576, 0.117062, 0.0926133, 0.0813332, 0.0739425, 0.0688675, 0.0639871, 0.0601296,
	0.0598187, 0.0577296, 0.0575132};

	Double_t g4_beamE_piplus_newbirk[13] = {1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40., 50., 60., 70., 80., 90., 100.};
	Double_t g4_res_piplus_5deg_newbirk[13] = {0.248197, 0.205701, 0.15331, 0.110433, 0.0870408, 0.0753077, 0.0679485, 0.0630292, 0.0602589, 0.057515, 0.0557724,
	0.0540185, 0.052538};

	double beam_E[10] = {3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20, 25, 32};
	double x_err[10] = {0,0,0,0,0,0,0,0,0};
	
	double factor[10] = {1.8, 2.0, 2.0, 1.8, 1.6, 1.6, 1.6, 1.2, 1.2, 1.2};
	double factor_err[10] = {0,0,0,0,0,0,0,0,0,0};
	
	double mean[10] = {487.1, 692.3, 1029.6, 1310.6, 1647.4, 1919.8, 2158.2 + 402.5/*2542.1*/, 2709.6 + 402.5 , 3364.5 + 402.5, 4413.8 + 402.5};
	double mean_err[10] = {0,0,0,0,0,0,0,0,0};
	
	//constant prop factor
	double mean_constant_factor[10] = {470.211, 621.5, 977.6, 1282.3, 1652.7, 1956.9, 2222.56 + 402.5, 2825.3 + 402.5, 3570.22 + 402.5, 4617.58 + 402.5};
	double mean_err_constantfactor[10] = {0,0,0,0,0,0,0,0,0};

	double sigma_constant_factor[10]= {154.9, 215.0, 231.8, 277.2, 345.6, 379.9, 433.185, 546.691, 642.301, 757.404};
	double sigma_err_constant_factor[10] = {0,0,0,0,0,0,0,0,0,0};
	double reso_err_constant_factor[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	
	//sigmas
	double sigma[10] = {165.0, 238.1, 247.6, 265.2, 329.5, 349.0, 402.5, 439.4, 536.3, 623.7};
	double sigma_err[10] = {0,0,0,0,0,0,0,0,0,0};
	double reso_err[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double ecal_reso_err[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	
	double ecal_mean[10] = {478.1, 622.4, 895.7, 1194.3, 1479.5, 1798.2, 2344.0, 2871.1, 3555.6, 4544.0};
	double ecal_mean_err[10] = {0,0,0,0,0,0,0,0,0};
	
	double ecal_sigma[10] = {37.0, 41.7, 51.5, 65.4, 68.5, 81.5, 104.5, 118.1, 137.6, 144.3};
	
	double beam_spread[10] = {.027,.023,.023,.023,.023,.023,.023,.023,.023,.023};


	
	double reso[10];
	double reso_constant_factor[10];
	double ecal_reso[10];
	for(int i=0;i<10;i++) {
		reso[i] = TMath::Sqrt(sigma[i]*sigma[i]/(mean[i]*mean[i]) - beam_spread[i]*beam_spread[i]);
		reso_constant_factor[i] = TMath::Sqrt(sigma_constant_factor[i]*sigma_constant_factor[i]/(mean_constant_factor[i]*mean_constant_factor[i])- beam_spread[i]*beam_spread[i]);
		//ecal
		ecal_reso[i] = TMath::Sqrt(ecal_sigma[i]*ecal_sigma[i]/(ecal_mean[i]*ecal_mean[i]) - beam_spread[i]*beam_spread[i]);
	}

	gStyle->SetOptFit(0000);


	TGraphErrors *gECalResolution = new TGraphErrors(10, beam_E, ecal_reso, 0, 0);
	gECalResolution->SetMarkerStyle(29);
	gECalResolution->SetMarkerColor(kGreen+2);
	gECalResolution->SetMarkerSize(2.4);
	gECalResolution->SetLineStyle(2);
	gECalResolution->SetLineColor(kGreen+2);
	
	TF1 *fresofitecal = new TF1("fresofitecal", "[1]/sqrt(x) + [0]", 1, 32);
	fresofitecal->SetLineStyle(2);
	fresofitecal->SetLineColor(kGreen+2);

	
	gECalResolution->Fit(fresofitecal, "R");
	gECalResolution->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gECalResolution->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gECalResolution->GetYaxis()->SetTitleOffset(1.4);

	
	TGraphErrors *gHCalResolution = new TGraphErrors(10, beam_E, reso, x_err, reso_err);
	gHCalResolution->SetMarkerStyle(33);
	gHCalResolution->SetMarkerColor(kRed);
	gHCalResolution->SetMarkerSize(2.4);
	gHCalResolution->SetLineStyle(2);
	gHCalResolution->SetLineColor(kRed);

	
	TGraphErrors *gSumFactor = new TGraphErrors(10, beam_E, factor, 0, factor_err);
	gSumFactor->SetMarkerStyle(20);
	gSumFactor->SetMarkerColor(kMagenta+2);
	gSumFactor->SetMarkerSize(2.1);

	TF1 *fresofit = new TF1("fresofit", "[1]/sqrt(x) + [0]", 1, 32);
	fresofit->SetLineStyle(2);
	fresofit->SetLineColor(kRed);



    TGraphErrors* gr_g4elec_res = new TGraphErrors(7, g4_beamE, g4_res_elec_5deg, 0, 0);
	gr_g4elec_res->SetMarkerStyle(kOpenStar);
	gr_g4elec_res->SetMarkerColor(kBlue);
	gr_g4elec_res->SetMarkerSize(2);
	gr_g4elec_res->SetLineStyle(2);
	gr_g4elec_res->SetLineColor(kBlue);
	
	TF1* f_g4elec_fit = new TF1("f_g4elec_fit", "[1]/sqrt(x) + [0]", 1, 32);
	f_g4elec_fit->SetParameter(0, .038);
	f_g4elec_fit->SetParameter(1, .43);
	f_g4elec_fit->SetLineStyle(2);
	f_g4elec_fit->SetLineColor(kBlue);

	gr_g4elec_res->Fit(f_g4elec_fit, "R");
	gr_g4elec_res->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gr_g4elec_res->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gr_g4elec_res->GetYaxis()->SetTitleOffset(1.4);

	TGraphErrors* gr_g4piplus_res = new TGraphErrors(14, g4_beamE_piplus, g4_res_piplus_5deg, 0, 0);
	gr_g4piplus_res->SetMarkerStyle(kOpenDiamond);
	gr_g4piplus_res->SetMarkerColor(kBlue);
	gr_g4piplus_res->SetMarkerSize(2);
	gr_g4piplus_res->SetLineStyle(2);
	gr_g4piplus_res->SetLineColor(kBlue);
	
	TF1* f_g4piplus_fit = new TF1("f_g4piplus_fit", "[1]/sqrt(x) + [0]", 3, 102);
	f_g4piplus_fit->SetParameter(0, .038);
	f_g4piplus_fit->SetParameter(1, .43);
	f_g4piplus_fit->SetLineStyle(2);
	f_g4piplus_fit->SetLineColor(kBlue);

	gr_g4piplus_res->Fit(f_g4piplus_fit, "R");
	gr_g4piplus_res->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gr_g4piplus_res->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gr_g4piplus_res->GetYaxis()->SetTitleOffset(1.4);


	TGraphErrors* gr_g4piplus_res_newbirk = new TGraphErrors(13, g4_beamE_piplus_newbirk, g4_res_piplus_5deg_newbirk, 0, 0);
	gr_g4piplus_res_newbirk->SetMarkerStyle(kOpenCircle);
	gr_g4piplus_res_newbirk->SetMarkerColor(kRed);
	gr_g4piplus_res_newbirk->SetMarkerSize(1.5);
	gr_g4piplus_res_newbirk->SetLineStyle(2);
	gr_g4piplus_res_newbirk->SetLineColor(kRed);
	
	TF1* f_g4piplus_fit_newbirk = new TF1("f_g4piplus_fit_newbirk", "[1]/sqrt(x) + [0]", 3, 102);
	f_g4piplus_fit_newbirk->SetParameter(0, .038);
	f_g4piplus_fit_newbirk->SetParameter(1, .43);
	f_g4piplus_fit_newbirk->SetLineStyle(2);
	f_g4piplus_fit_newbirk->SetLineColor(kRed);

	gr_g4piplus_res_newbirk->Fit(f_g4piplus_fit_newbirk, "R");
	gr_g4piplus_res_newbirk->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gr_g4piplus_res_newbirk->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gr_g4piplus_res_newbirk->GetYaxis()->SetTitleOffset(1.4);

	
 		
	//res from simulation
	TF1 *MC_hadron_resolution = new TF1("MC_hadron_resolution", "[1]/sqrt(x) + [0]", 1, 32);
	MC_hadron_resolution->SetParameter(0, .038);
	MC_hadron_resolution->SetParameter(1, .43);
	MC_hadron_resolution->SetLineStyle(3);
	MC_hadron_resolution->SetLineColor(kBlue);
	MC_hadron_resolution->SetLineWidth(2);
	
	TCanvas *cHCalResolution = new TCanvas("cHCalResolution", "HCal Resolution", 900, 700);
	cHCalResolution->cd();
	// gHCalResolution->Fit(fresofit, "R");

	gECalResolution->Draw("AP");
    gr_g4elec_res->Draw("sameP");
	gr_g4piplus_res->Draw("sameP");
	// gHCalResolution->Draw("sameP");
	gr_g4piplus_res_newbirk->Draw("sameP");


	gECalResolution->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gECalResolution->GetXaxis()->SetLimits(0.0, 105.0);
	gECalResolution->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gECalResolution->GetYaxis()->SetTitleOffset(1.4);
	gECalResolution->GetYaxis()->SetRangeUser(0.0, 0.5);
	gECalResolution->SetTitle("");
	
	// MC_hadron_resolution->Draw("sameL");
	
	TLegend *leg = new TLegend(.3, .5, .85, .85);
	leg->SetHeader("Resolution of ECAL + HCAL");
	// leg->AddEntry(gHCalResolution, "hadron,  #frac{#DeltaE}{E} = #frac{0.58 #pm 0.04}{#sqrt{E}} + 0.015 #pm 0.015", "LP");
	TString piplus_text;
	piplus_text.Form("G4 #pi^{+}, kB = .126, #frac{#DeltaE}{E} =  #frac{%0.3f}{#sqrt{E}} + %0.3f", f_g4piplus_fit->GetParameter(1), f_g4piplus_fit->GetParameter(0));
	leg->AddEntry(gr_g4piplus_res, piplus_text, "LP");

	TString piplus_text_newbirk;
	piplus_text_newbirk.Form("G4 #pi^{+}, kB = .078, #frac{#DeltaE}{E} =  #frac{%0.3f}{#sqrt{E}} + %0.3f", f_g4piplus_fit_newbirk->GetParameter(1), f_g4piplus_fit_newbirk->GetParameter(0));
	leg->AddEntry(gr_g4piplus_res_newbirk, piplus_text_newbirk, "LP");
	
	leg->AddEntry(gECalResolution, "electron,  #frac{#DeltaE}{E} =  #frac{0.11 #pm 0.007}{#sqrt{E}} + 0.007 #pm 0.002", "LP");

	TString electron_text;
	electron_text.Form("G4 electron, #frac{#DeltaE}{E} =  #frac{%0.2f}{#sqrt{E}} + %0.3f ", f_g4elec_fit->GetParameter(1), f_g4elec_fit->GetParameter(0));
	leg->AddEntry(gr_g4elec_res, electron_text, "LP");

	

	leg->SetTextSize(.028);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->Draw();
	
	cHCalResolution->SetRightMargin(.02);
	gStyle->SetStatX(.8);
	gStyle->SetStatY(.85);
	gStyle->SetStatW(.15);
	gStyle->SetStatH(.13);
	cHCalResolution->Update();

}