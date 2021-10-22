#include <iostream>
#include <string>
#include "TTree.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TGraph.h"
#include <map>
#include <vector>
#include <ROOT/TProcessExecutor.hxx>
#include <tuple>


std::map<Double_t, Double_t> Beam_MaxEnergy = {{1.0, 100}, {2.0, 100}, {5.0, 200}, {10.0, 400}, {20.0, 700}, {30.0, 1000}, {40.0, 1500}, {50.0, 2000}, {60., 2500}, {70., 3000}, {80., 3500}, {90., 4000}, {100., 4500}};

TH1D* ECalWeighting(TTree* tree, Double_t energy, Double_t weight)
{
    Double_t ECalEdep, HCalEdep;
    tree->SetBranchAddress("ECal_Edep_Total", &ECalEdep);
    tree->SetBranchAddress("HCal_Edep_Total", &HCalEdep);
    Int_t num_events = (Int_t) tree->GetEntries();
    char hist_name[100];
    sprintf(hist_name, "h_TotalEdep%f", weight);

    TH1D* h_TotalEdep = new TH1D(hist_name, "", 600, 0, Beam_MaxEnergy[energy]);
    for(Int_t i = 0; i < num_events; i++)
    {
        tree->GetEntry(i);
        if(ECalEdep < 0.183) ECalEdep = 0.;
        h_TotalEdep->Fill(ECalEdep/weight + HCalEdep);        
    }

    return h_TotalEdep;
}

Double_t ECalWeightingProcess(TTree* tree, Double_t energy)
{
    std::vector<Double_t> weights;
    Double_t w0 = 0.6;
    for(Int_t i = 0; i < 20; i++)
    {
        weights.push_back(w0 + i*0.05);
    }

    std::vector<TH1D*> WeightedHists;

    for(auto& weight : weights)
    {
        WeightedHists.push_back(ECalWeighting(tree, energy, weight));
    }

    TCanvas* c_WeightedHists = new TCanvas("c_WeightedHists", "", 1000, 1000);
    c_WeightedHists->Divide(4,5);

    Double_t minimum_resolution = 99999;
    Double_t optimal_weight = -1.;
    std::vector<Double_t> resolutions;

    for(Int_t i = 0; i < WeightedHists.size(); i++)
    {
        c_WeightedHists->cd(i+1);
        WeightedHists[i]->Draw("same");
        WeightedHists[i]->Fit("gaus", "q");
        Double_t mean = WeightedHists[i]->GetFunction("gaus")->GetParameter(1);
        Double_t sigma = WeightedHists[i]->GetFunction("gaus")->GetParameter(2);
        Double_t resolution;
        if(mean != 0) resolution = sigma/mean;
        else resolution = 99999.;
        resolutions.push_back(resolution);

        if( resolution < minimum_resolution)
        {
            minimum_resolution = resolution;
            optimal_weight = weights[i];
        } 
    }

    TCanvas* c_Weights = new TCanvas("c_Weights", "", 1000, 1000);
    TGraph* gr_Weights = new TGraph(weights.size(), weights.data(), resolutions.data());
    gr_Weights->Draw("ap");
    gr_Weights->SetMarkerStyle(kOpenCircle);
    gr_Weights->SetMarkerColor(kBlue);
    gr_Weights->GetXaxis()->SetTitle("Weight");
    gr_Weights->GetYaxis()->SetTitle("Resolution");
    gr_Weights->SetTitle("");

    std::cout<<"Optimal weight is "<<optimal_weight<<std::endl;
    return optimal_weight;
}
void Resolution(std::string particle = "e-", Double_t energy = 1.0)
{   

    // ROOT::EnableThreadSafety();
    // Int_t num_threads = 6;
    // ROOT::EnableImplicitMT(num_threads);

    
    Bool_t ECal_weight = kTRUE;
    if(particle == "e-" || particle == "pi0") ECal_weight = kFALSE;
    if(ECal_weight) std::cout<<"Using weighting for ECal."<<std::endl;
    else std::cout<<"Not using weighting for ECal."<<std::endl;

    TString file_name;
    //file_name.Form("test.root");
    // file_name.Form("%s5deg_Data/%s_%0.0fGeV.root", particle.c_str(), particle.c_str(), energy);
    file_name.Form("%s_%0.0fGeV.root", particle.c_str(), energy);
    std::cout<<"Opening "<<file_name<<std::endl;
    TFile* data_file = new TFile(file_name);

    TH1D* h_TotalEdep = new TH1D("h_TotalEdep", "", 600, 0, Beam_MaxEnergy[energy]);
    TTree* Total_tree = (TTree*) data_file->Get("EdepTotal");
    Int_t num_events = (Int_t) Total_tree->GetEntries();
    std::cout<<"Number of events: "<<num_events<<std::endl;

    Double_t optimal_weight = 1.0;
    if(ECal_weight) optimal_weight = ECalWeightingProcess(Total_tree, energy);

    Double_t ECalEdep, HCalEdep;
    Total_tree->SetBranchAddress("ECal_Edep_Total", &ECalEdep);
    Total_tree->SetBranchAddress("HCal_Edep_Total", &HCalEdep);

    for(Int_t i = 0; i < num_events; i++)
    {
        Total_tree->GetEntry(i);
        if(ECalEdep < 0.183) ECalEdep = 0.;
        Double_t ECal_energy = ECalEdep/optimal_weight;
        Double_t total_energy = ECal_energy + HCalEdep;
        h_TotalEdep->Fill(total_energy);
    }
    
    TCanvas* c_Resolution = new TCanvas("c_Resolution", "", 1000, 1000);
    h_TotalEdep->Draw();
    h_TotalEdep->GetXaxis()->SetTitle("Edep (MeV)");
    h_TotalEdep->GetYaxis()->SetTitle("Number of Events");

    TF1* f_gaus = new TF1("f_gaus", "gaus", 0, Beam_MaxEnergy[energy]);
    h_TotalEdep->Fit(f_gaus, "");

    Double_t mean = f_gaus->GetParameter(1);
    Double_t sigma = f_gaus->GetParameter(2);
    Double_t resolution;
    if(mean != 0) resolution = sigma/mean;
    else resolution = 0.;

    std::cout<<"Resolution is "<<resolution<<std::endl;

    TString info_text;
    info_text = Form("%s at %0.0f GeV", particle.c_str(), energy);

    h_TotalEdep->SetTitle(info_text);

    TString res_text;
    res_text = Form("Resolution = %0.5f", resolution);

    TLatex info_caption;
    info_caption.SetTextFont(62);
    info_caption.SetTextSize(.04);
    info_caption.SetNDC(kTRUE);
    info_caption.DrawLatex(.15, .85, info_text);

    TLatex res_caption;
    res_caption.SetTextFont(62);
    res_caption.SetTextSize(.04);
    res_caption.SetNDC(kTRUE);
    res_caption.DrawLatex(.15, .8, res_text);
}