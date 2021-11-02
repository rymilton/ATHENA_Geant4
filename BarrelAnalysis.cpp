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
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMultiGraph.h"


const Int_t num_total_layers = 51;
const Int_t num_ignored_layers = 11;
const Int_t num_net_layers = num_total_layers - num_ignored_layers;
std::map<Double_t, Double_t> Beam_MaxEnergy = {{1.0, 100}, {2.0, 100}, {5.0, 200}, {10.0, 400}, {20.0, 700}, {30.0, 1000}};
std::map<Double_t, Int_t> Energy_NumHits = {{1.0, 100}, {2.0, 200}, {3.0, 300}, {4.0, 400}, {5.0, 500}, {6.0, 600}, {7.0, 700}, {8.0, 800}, {9.0, 900}, {10.0, 1000}};

void BarrelAnalysis(Double_t energy = 1.0)
{
    TString file_name;
    file_name.Form("./neutron5deg_Data/neutron_%0.0fGeV.root", energy);
    std::cout<<"Opening "<<file_name<<std::endl;
    TFile* data_file = new TFile(file_name);

    TTree* Total_tree = (TTree*) data_file->Get("EdepTotal");
    TTree* HCal_tree = (TTree*) data_file->Get("HCalLayers");
    Int_t num_events = (Int_t) Total_tree->GetEntries();
    std::cout<<"Number of events: "<<num_events<<std::endl;

    Double_t ECalEdep;
    Int_t ECalHits;
    Int_t ECaleventID;
    Total_tree->SetBranchAddress("ECal_Edep_Total", &ECalEdep);
    Total_tree->SetBranchAddress("ECal_NumHits_Total", &ECalHits);
    Total_tree->SetBranchAddress("eventID", &ECaleventID);

    Int_t HCalTileHits;
    Double_t HCalTileEdep; 
    Int_t layerid;
    Int_t HCaleventID;
    HCal_tree->SetBranchAddress("HCal_Edep_Tile", &HCalTileEdep);
    HCal_tree->SetBranchAddress("HCal_NumHits_Tile", &HCalTileHits);
    HCal_tree->SetBranchAddress("HCal_Layerid", &layerid);
    HCal_tree->SetBranchAddress("eventID", &HCaleventID);

    TProfile* p_layerEdep = new TProfile("p_layerEdep", "", num_net_layers, 0.5, num_net_layers+.5);
    TProfile* p_layerHits = new TProfile("p_layerHits", "", num_net_layers, 0.5, num_net_layers+.5);
    TProfile* p_tileHits = new TProfile("p_tileHits", "", num_net_layers, 0.5, num_net_layers+.5);

    TH1I* h_EO1220Hits = new TH1I("h_EO1220Hits", "", Energy_NumHits[energy], 0, Energy_NumHits[energy]); // every other layer hits
    TH1I* h_all1220Hits = new TH1I("h_all1220Hits", "", Energy_NumHits[energy]*2, 0, Energy_NumHits[energy]*2);
    TH1I* h_EO1220Tile = new TH1I("h_EO1220Tile", "", 60, 0, 60); // every other layer tile hits
    TH1I* h_all1220Tiles = new TH1I("h_all1220Tiles", "", 100, 0, 100);
    TH1D* h_ECalEdep = new TH1D("h_ECalEdep", "", 200, 0, 100);

    TH1D* h_efficiency = new TH1D("h_efficiency", " ", 10, 0.5, 10.5);

    Double_t ECalEdep_array[num_events];
    Double_t HCalEdep_array[num_events];
    Int_t ECalHitsArray[num_events];

    for(Int_t i = 0; i< num_events; i++)
    {
        ECalEdep_array[i] = 0.;
        HCalEdep_array[i] = 0.;
        ECalHitsArray[i] = 0;
    }

    Double_t efficiencyWholeAll1220 = 0.;
    Double_t efficiencyWholeEO1220 = 0.;
    Double_t efficiencyWholeAllLayers = 0.;
    Double_t efficiencyWholeEO1232 = 0.;
    Double_t efficiencyWholeDead = 0.;

    Double_t efficiencyHCalAll1220 = 0.;
    Double_t efficiencyHCalEO1220 = 0.;
    Double_t efficiencyHCalAllLayers = 0.;
    Double_t efficiencyHCalEO1232 = 0.;

    Double_t efficiencyECal = 0.;
    Double_t avg_ECalEdep = 0.;

    Double_t ECalCut = 1.0;
    for(Int_t i = 0; i < num_events; i++)
    {
        Total_tree->GetEntry(i);
        
        h_ECalEdep->Fill(ECalEdep);
        if(ECalEdep < ECalCut)
        {
            ECalEdep = 0.;
            ECalHits = 0;
        }

        ECalEdep_array[ECaleventID] = ECalEdep;
        ECalHitsArray[ECaleventID] = ECalHits;
        avg_ECalEdep += ECalEdep;
        // h_ECalEdep->Fill(ECalEdep);
    }

    for(Int_t i = 0; i < num_events; i++)
    {
        Int_t EO1220Hits = 0;
        Int_t total1220Hits = 0;
        Int_t EO1232Hits = 0;
        Int_t totalHits = 0;
        Int_t totalwithDeadHits = 0;
        Double_t EO1220Tiles = 0;
        Double_t total1220Tiles = 0;
        Double_t HCalEventEdep = 0.;

        Double_t HCalLayerEdep[num_net_layers];
        Int_t HCalLayerHitsarray[num_net_layers];
        Int_t HCalTileHitsarray[num_net_layers];
        
        for(Int_t ilayer = 0; ilayer < num_net_layers; ilayer++)
        {
            HCalLayerEdep[ilayer] = 0.;
            HCalLayerHitsarray[ilayer] = 0;
            HCalTileHitsarray[ilayer] = 0;
        }

        for(Int_t itile = 0; itile < 36*num_total_layers; itile++)
        {
            HCal_tree->GetEntry(itile + i*36*num_total_layers);

            if(HCalTileEdep < 0.1)
            {
                HCalTileHits = 0;
                HCalTileEdep = 0.;
            }
            totalwithDeadHits += HCalTileHits;
            if(layerid < 11) continue;  

            totalHits += HCalTileHits;
            if(layerid == 21 || layerid == 23 || layerid == 25 || layerid == 27 || layerid == 29 || layerid == 31)
            {
                EO1232Hits += HCalTileHits;
            }
            else if(layerid == 11 || layerid == 13 || layerid == 15 || layerid == 17 || layerid == 19)
            {
                EO1220Hits += HCalTileHits;
                total1220Hits += HCalTileHits;
                EO1232Hits += HCalTileHits;
                if(HCalTileHits != 0)
                {
                    EO1220Tiles++;
                    total1220Tiles++;
                }
            }
            else if (layerid >= 11 && layerid <= 19)
            {
                total1220Hits += HCalTileHits;
                if(HCalTileHits != 0) total1220Tiles++;
            }
            HCalLayerEdep[layerid - num_ignored_layers] += HCalTileEdep;
            HCalLayerHitsarray[layerid - num_ignored_layers] += HCalTileHits;
            if(HCalTileHits != 0) HCalTileHitsarray[layerid - num_ignored_layers]++;
            HCalEventEdep += HCalTileEdep;
        }
        for(Int_t ilayer = 0; ilayer < num_net_layers; ilayer++)
        {
            p_layerEdep->Fill(ilayer+1, HCalLayerEdep[ilayer]);
            p_layerHits->Fill(ilayer+1, HCalLayerHitsarray[ilayer]); 
            p_tileHits->Fill(ilayer+1, HCalTileHitsarray[ilayer]);
        }
        
        HCalEdep_array[HCaleventID] += HCalEventEdep;
        Bool_t ECalHit;
        if(ECalHitsArray[HCaleventID] == 0) ECalHit = kFALSE;
        else ECalHit = kTRUE;

        if(totalwithDeadHits != 0 || ECalHit) efficiencyWholeDead++;
        if(totalHits != 0 || ECalHit) efficiencyWholeAllLayers++;
        if(EO1232Hits != 0 || ECalHit) efficiencyWholeEO1232++;
        if(total1220Hits != 0 || ECalHit) efficiencyWholeAll1220++;
        if(EO1220Hits != 0 || ECalHit) efficiencyWholeEO1220++;
        
        if(totalHits != 0) efficiencyHCalAllLayers++;
        if(EO1232Hits != 0) efficiencyHCalEO1232++;
        if(total1220Hits != 0 ) efficiencyHCalAll1220++;
        if(EO1220Hits != 0) efficiencyHCalEO1220++;
        
        if(ECalHit) efficiencyECal++;

        h_EO1220Hits->Fill(EO1220Hits);
        h_all1220Hits->Fill(total1220Hits);
        h_EO1220Tile->Fill(EO1220Tiles);
        h_all1220Tiles->Fill(total1220Tiles);
    }


    efficiencyWholeAll1220 /= num_events;
    efficiencyWholeAll1220 *= 100;
    efficiencyWholeAllLayers /= num_events;
    efficiencyWholeAllLayers *= 100;
    efficiencyWholeEO1220 /= num_events;
    efficiencyWholeEO1220 *= 100;
    efficiencyWholeEO1232 /= num_events;
    efficiencyWholeEO1232 *= 100;
    efficiencyHCalAll1220 /= num_events;
    efficiencyHCalAll1220 *= 100;
    efficiencyHCalAllLayers /= num_events;
    efficiencyHCalAllLayers *= 100;
    efficiencyHCalEO1220 /= num_events;
    efficiencyHCalEO1220 *= 100;
    efficiencyHCalEO1232 /= num_events;
    efficiencyHCalEO1232 *= 100;
    efficiencyECal /= num_events;
    efficiencyECal *= 100;
    efficiencyWholeDead /= num_events;
    efficiencyWholeDead *= 100;
    

    h_efficiency->SetBinContent(1, efficiencyWholeAllLayers);
    h_efficiency->SetBinContent(2, efficiencyWholeEO1232);
    h_efficiency->SetBinContent(3, efficiencyWholeAll1220);
    h_efficiency->SetBinContent(4, efficiencyWholeEO1220);
    h_efficiency->SetBinContent(5, efficiencyHCalAllLayers);
    h_efficiency->SetBinContent(6, efficiencyHCalEO1232);
    h_efficiency->SetBinContent(7, efficiencyHCalAll1220);
    h_efficiency->SetBinContent(8, efficiencyHCalEO1220);
    h_efficiency->SetBinContent(9, efficiencyECal);
    h_efficiency->SetBinContent(10, efficiencyWholeDead);


    TCanvas* c_LayerEdep = new TCanvas("c_LayerEdep", "", 1000, 1000);
    p_layerEdep->Draw("EX0");
    p_layerEdep->SetMarkerStyle(kFullCircle);
    p_layerEdep->SetMarkerColor(kBlack);
    p_layerEdep->SetTitle("");
    p_layerEdep->GetXaxis()->SetTitle("Layer number");
    p_layerEdep->GetYaxis()->SetTitle("Mean Energy deposit (MeV)");

    TCanvas* c_LayerHits = new TCanvas("c_LayerHits", "", 1000, 1000);
    c_LayerHits->Divide(1,2);
    c_LayerHits->cd(1);

    h_EO1220Hits->Draw();
    h_EO1220Hits->GetXaxis()->SetTitle("Number of hits, even layers in [12,20]");
    h_EO1220Hits->GetYaxis()->SetTitle("Number of Events");
    TString info_text;
    info_text = Form("Neutron, %0.0f GeV, 5#circ", energy);
    h_EO1220Hits->SetTitle(info_text);
    
    c_LayerHits->cd(2);
    h_all1220Hits->Draw("same");
    h_all1220Hits->GetXaxis()->SetTitle("Number of hits, all layers in [12,20]");
    h_all1220Hits->GetYaxis()->SetTitle("Number of Events");
    


    TCanvas* c_TileHits = new TCanvas("c_TileHits", "", 1000, 1000);
    c_TileHits->Divide(1,2);
    c_TileHits->cd(1);

    h_EO1220Tile->Draw();
    h_EO1220Tile->SetTitle(info_text);
    h_EO1220Tile->GetXaxis()->SetTitle("Number of tiles hit, even layers in [12,20]");
    h_EO1220Tile->GetYaxis()->SetTitle("Number of Events");

    c_TileHits->cd(2);
    h_all1220Tiles->Draw("same");
    h_all1220Tiles->GetXaxis()->SetTitle("Number of tiles hits, all layers in [12,20]");
    h_all1220Tiles->GetYaxis()->SetTitle("Number of Events");

    TCanvas* c_ECalEdep = new TCanvas("c_ECalEdep", "", 1200, 800);

    c_ECalEdep->SetLogy();
    h_ECalEdep->Draw();
    h_ECalEdep->GetXaxis()->SetTitle("ECal Edep (MeV)");
    h_ECalEdep->GetYaxis()->SetTitle("Number of Events (Log)");
    h_ECalEdep->SetTitle("");
    h_ECalEdep->SetAxisRange(0, 50, "x");
    h_ECalEdep->GetXaxis()->SetTickLength(.02);
    h_ECalEdep->GetXaxis()->SetTitleOffset(1.3);
    h_ECalEdep->GetYaxis()->SetTickLength(-0.02);
    h_ECalEdep->GetYaxis()->SetLabelOffset(.02);
    h_ECalEdep->GetYaxis()->SetTitleOffset(1.5);
    gStyle->SetOptStat(0);
    gPad->SetTickx();

    TLatex info_caption;
    info_caption.SetNDC(kTRUE);

    info_caption.SetTextFont(62);
    info_caption.SetTextSize(.05);
    info_caption.DrawLatex(.4, .84, info_text);
    
    TLine* l_cut = new TLine(ECalCut, 0.0, ECalCut, 20000);
    l_cut->SetLineWidth(2);
    l_cut->SetLineStyle(2);
    l_cut->SetLineColor(kRed);
    l_cut->Draw("same");


    TCanvas* c_LayerAvgs = new TCanvas("c_LayerAvgs", "", 1200, 800);
    c_LayerAvgs->Divide(1,2);
    c_LayerAvgs->cd(1);

    p_layerHits->Draw("EX0");
    p_layerHits->SetMarkerStyle(kFullCircle);
    p_layerHits->SetMarkerColor(kBlack);
    p_layerHits->SetTitle("");
    p_layerHits->GetXaxis()->SetTitle("Layer number");
    p_layerHits->GetYaxis()->SetTitle("Mean Number of Hits");

    c_LayerAvgs->cd(2);
    p_tileHits->Draw("EX0");
    p_tileHits->SetMarkerStyle(kFullCircle);
    p_tileHits->SetMarkerColor(kBlack);
    p_tileHits->SetTitle("");
    p_tileHits->GetXaxis()->SetTitle("Layer number");
    p_tileHits->GetYaxis()->SetTitle("Mean Number of Tiles Hit");

    TString outname;
    outname.Form("./neutron5deg_Data/AnalyzedNeutron_%0.0fGeV.root", energy);
    TFile* outfile = new TFile(outname, "RECREATE");
    p_layerEdep->Write();
    p_tileHits->Write();
    p_layerHits->Write();
    h_EO1220Hits->Write();
    h_all1220Hits->Write();
    h_EO1220Tile->Write();
    h_all1220Tiles->Write();
    h_efficiency->Write();

    
    std::cout<<"Efficiency for detector, all layers: "<<efficiencyWholeAllLayers<<std::endl;
    std::cout<<"Efficiency for detector, EO layers in [12,32]: "<<efficiencyWholeEO1232<<std::endl;
    std::cout<<"Efficiency for detector, all layers in [12,20]: "<<efficiencyWholeAll1220<<std::endl;
    std::cout<<"Efficiency for detector, EO layers in [12,20]: "<<efficiencyWholeEO1220<<std::endl;
    std::cout<<"Efficiency for HCal, EO layers in [12,20]: "<<efficiencyHCalEO1220<<std::endl;
    std::cout<<"Efficiency for ECal: "<<efficiencyECal<<std::endl;
    std::cout<<"Average ECal Edep: "<<avg_ECalEdep/num_events<<std::endl;
}

void BarrelPlot()
{
    TFile* neutron1GeV_file = new TFile("./neutron5deg_Data/AnalyzedNeutron_1GeV.root");
    TFile* neutron2GeV_file = new TFile("./neutron5deg_Data/AnalyzedNeutron_2GeV.root");
    TFile* neutron3GeV_file = new TFile("./neutron5deg_Data/AnalyzedNeutron_3GeV.root");
    TFile* neutron4GeV_file = new TFile("./neutron5deg_Data/AnalyzedNeutron_4GeV.root");
    TFile* neutron5GeV_file = new TFile("./neutron5deg_Data/AnalyzedNeutron_5GeV.root");
    TFile* neutron6GeV_file = new TFile("./neutron5deg_Data/AnalyzedNeutron_6GeV.root");
    TFile* neutron7GeV_file = new TFile("./neutron5deg_Data/AnalyzedNeutron_7GeV.root");
    TFile* neutron8GeV_file = new TFile("./neutron5deg_Data/AnalyzedNeutron_8GeV.root");
    TFile* neutron9GeV_file = new TFile("./neutron5deg_Data/AnalyzedNeutron_9GeV.root");
    TFile* neutron10GeV_file = new TFile("./neutron5deg_Data/AnalyzedNeutron_10GeV.root");

    TProfile* neutron1GeV_edep;
    TProfile* neutron2GeV_edep;
    TProfile* neutron3GeV_edep;
    TProfile* neutron4GeV_edep;
    TProfile* neutron5GeV_edep;
    TProfile* neutron6GeV_edep;
    TProfile* neutron7GeV_edep;
    TProfile* neutron8GeV_edep;
    TProfile* neutron9GeV_edep;
    TProfile* neutron10GeV_edep;

    TProfile* neutron1GeV_TileHits;
    TProfile* neutron2GeV_TileHits;
    TProfile* neutron3GeV_TileHits;
    TProfile* neutron4GeV_TileHits;
    TProfile* neutron5GeV_TileHits;
    TProfile* neutron6GeV_TileHits;
    TProfile* neutron7GeV_TileHits;
    TProfile* neutron8GeV_TileHits;
    TProfile* neutron9GeV_TileHits;
    TProfile* neutron10GeV_TileHits;

    TProfile* neutron1GeV_numHits;
    TProfile* neutron2GeV_numHits;
    TProfile* neutron3GeV_numHits;
    TProfile* neutron4GeV_numHits;
    TProfile* neutron5GeV_numHits;
    TProfile* neutron6GeV_numHits;
    TProfile* neutron7GeV_numHits;
    TProfile* neutron8GeV_numHits;
    TProfile* neutron9GeV_numHits;
    TProfile* neutron10GeV_numHits;

    TH1D* neutron1GeV_efficiency;
    TH1D* neutron2GeV_efficiency;
    TH1D* neutron3GeV_efficiency;
    TH1D* neutron4GeV_efficiency;
    TH1D* neutron5GeV_efficiency;
    TH1D* neutron6GeV_efficiency;
    TH1D* neutron7GeV_efficiency;
    TH1D* neutron8GeV_efficiency;
    TH1D* neutron9GeV_efficiency;
    TH1D* neutron10GeV_efficiency;

    TH1I* h_neutron1GeV_TileHits;
    TH1I* h_neutron2GeV_TileHits;
    TH1I* h_neutron3GeV_TileHits;
    TH1I* h_neutron4GeV_TileHits;
    TH1I* h_neutron5GeV_TileHits;
    TH1I* h_neutron6GeV_TileHits;
    TH1I* h_neutron7GeV_TileHits;
    TH1I* h_neutron8GeV_TileHits;
    TH1I* h_neutron9GeV_TileHits;
    TH1I* h_neutron10GeV_TileHits;

    neutron1GeV_file->GetObject("p_layerEdep", neutron1GeV_edep);
    neutron2GeV_file->GetObject("p_layerEdep", neutron2GeV_edep);
    neutron3GeV_file->GetObject("p_layerEdep", neutron3GeV_edep);
    neutron4GeV_file->GetObject("p_layerEdep", neutron4GeV_edep);
    neutron5GeV_file->GetObject("p_layerEdep", neutron5GeV_edep);
    neutron6GeV_file->GetObject("p_layerEdep", neutron6GeV_edep);
    neutron7GeV_file->GetObject("p_layerEdep", neutron7GeV_edep);
    neutron8GeV_file->GetObject("p_layerEdep", neutron8GeV_edep);
    neutron9GeV_file->GetObject("p_layerEdep", neutron9GeV_edep);
    neutron10GeV_file->GetObject("p_layerEdep", neutron10GeV_edep);

    neutron1GeV_file->GetObject("p_layerHits", neutron1GeV_numHits);
    neutron2GeV_file->GetObject("p_layerHits", neutron2GeV_numHits);
    neutron3GeV_file->GetObject("p_layerHits", neutron3GeV_numHits);
    neutron4GeV_file->GetObject("p_layerHits", neutron4GeV_numHits);
    neutron5GeV_file->GetObject("p_layerHits", neutron5GeV_numHits);
    neutron6GeV_file->GetObject("p_layerHits", neutron6GeV_numHits);
    neutron7GeV_file->GetObject("p_layerHits", neutron7GeV_numHits);
    neutron8GeV_file->GetObject("p_layerHits", neutron8GeV_numHits);
    neutron9GeV_file->GetObject("p_layerHits", neutron9GeV_numHits);
    neutron10GeV_file->GetObject("p_layerHits", neutron10GeV_numHits);

    neutron1GeV_file->GetObject("p_tileHits", neutron1GeV_TileHits);
    neutron2GeV_file->GetObject("p_tileHits", neutron2GeV_TileHits);
    neutron3GeV_file->GetObject("p_tileHits", neutron3GeV_TileHits);
    neutron4GeV_file->GetObject("p_tileHits", neutron4GeV_TileHits);
    neutron5GeV_file->GetObject("p_tileHits", neutron5GeV_TileHits);
    neutron6GeV_file->GetObject("p_tileHits", neutron6GeV_TileHits);
    neutron7GeV_file->GetObject("p_tileHits", neutron7GeV_TileHits);
    neutron8GeV_file->GetObject("p_tileHits", neutron8GeV_TileHits);
    neutron9GeV_file->GetObject("p_tileHits", neutron9GeV_TileHits);
    neutron10GeV_file->GetObject("p_tileHits", neutron10GeV_TileHits);

    neutron1GeV_file->GetObject("h_efficiency", neutron1GeV_efficiency);
    neutron2GeV_file->GetObject("h_efficiency", neutron2GeV_efficiency);
    neutron3GeV_file->GetObject("h_efficiency", neutron3GeV_efficiency);
    neutron4GeV_file->GetObject("h_efficiency", neutron4GeV_efficiency);
    neutron5GeV_file->GetObject("h_efficiency", neutron5GeV_efficiency);
    neutron6GeV_file->GetObject("h_efficiency", neutron6GeV_efficiency);
    neutron7GeV_file->GetObject("h_efficiency", neutron7GeV_efficiency);
    neutron8GeV_file->GetObject("h_efficiency", neutron8GeV_efficiency);
    neutron9GeV_file->GetObject("h_efficiency", neutron9GeV_efficiency);
    neutron10GeV_file->GetObject("h_efficiency", neutron10GeV_efficiency);

    neutron1GeV_file->GetObject("h_EO1220Tile", h_neutron1GeV_TileHits);
    neutron2GeV_file->GetObject("h_EO1220Tile", h_neutron2GeV_TileHits);
    neutron3GeV_file->GetObject("h_EO1220Tile", h_neutron3GeV_TileHits);
    neutron4GeV_file->GetObject("h_EO1220Tile", h_neutron4GeV_TileHits);
    neutron5GeV_file->GetObject("h_EO1220Tile", h_neutron5GeV_TileHits);
    neutron6GeV_file->GetObject("h_EO1220Tile", h_neutron6GeV_TileHits);
    neutron7GeV_file->GetObject("h_EO1220Tile", h_neutron7GeV_TileHits);
    neutron8GeV_file->GetObject("h_EO1220Tile", h_neutron8GeV_TileHits);
    neutron9GeV_file->GetObject("h_EO1220Tile", h_neutron9GeV_TileHits);
    neutron10GeV_file->GetObject("h_EO1220Tile", h_neutron10GeV_TileHits);

    

    Double_t energy[10] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
    std::vector<Double_t> efficiencyWholeAll;
    std::vector<Double_t> efficiencyWholeEO1232;
    std::vector<Double_t> efficiencyWholeAll1220;
    std::vector<Double_t> efficiencyWholeEO1220;
    std::vector<Double_t> efficiencyHCalAll;
    std::vector<Double_t> efficiencyHCalEO1232;
    std::vector<Double_t> efficiencyHCalAll1220;
    std::vector<Double_t> efficiencyHCalEO1220;
    std::vector<Double_t> efficiencyECal;
    std::vector<Double_t> efficiencywithoutDead;


    std::map<Int_t, TH1D*> efficiencyMap = {{1, neutron1GeV_efficiency}, {2, neutron2GeV_efficiency}, {3, neutron3GeV_efficiency}, {4, neutron4GeV_efficiency}, {5, neutron5GeV_efficiency},
    {6, neutron6GeV_efficiency}, {7, neutron7GeV_efficiency}, {8, neutron8GeV_efficiency}, {9, neutron9GeV_efficiency}, {10, neutron10GeV_efficiency}};


    for(Int_t i = 1; i < 11; i++)
    {
        efficiencyWholeAll.push_back(efficiencyMap[i]->GetBinContent(1));
        efficiencyWholeEO1232.push_back(efficiencyMap[i]->GetBinContent(2));
        efficiencyWholeAll1220.push_back(efficiencyMap[i]->GetBinContent(3));
        efficiencyWholeEO1220.push_back(efficiencyMap[i]->GetBinContent(4));
        efficiencyHCalAll.push_back(efficiencyMap[i]->GetBinContent(5));
        efficiencyHCalEO1232.push_back(efficiencyMap[i]->GetBinContent(6));
        efficiencyHCalAll1220.push_back(efficiencyMap[i]->GetBinContent(7));
        efficiencyHCalEO1220.push_back(efficiencyMap[i]->GetBinContent(8));
        efficiencyECal.push_back(efficiencyMap[i]->GetBinContent(9));
        efficiencywithoutDead.push_back(efficiencyMap[i]->GetBinContent(10));
    }

    

    Double_t ECalEdep[10] = {10.1014, 18.8301, 27.069, 35.2279, 43.3898, 51.2522, 58.5576, 65.5018, 71.8679, 78.6714};
    TGraph* gr_ECalEdep = new TGraph(10, energy, ECalEdep);

    TMultiGraph* mg_efficiency = new TMultiGraph();

    TGraph* gr_efficiencyWholeAll = new TGraph(10, energy, efficiencyWholeAll.data());
    TGraph* gr_efficiencyWholeEO1232 = new TGraph(10, energy, efficiencyWholeEO1232.data());
    TGraph* gr_efficiencyWholeAll1220 = new TGraph(10, energy, efficiencyWholeAll1220.data());
    TGraph* gr_efficiencyWholeEO1220 = new TGraph(10, energy, efficiencyWholeEO1220.data());
    TGraph* gr_efficiencywithoutDead = new TGraph(10, energy, efficiencywithoutDead.data());
    


    

    TCanvas* c_ECalEdep = new TCanvas("c_EcalEdep", "", 1200, 800);
    gr_ECalEdep->Draw("ap");
    gr_ECalEdep->SetMarkerStyle(kFullCircle);
    gr_ECalEdep->SetMarkerColor(kBlack);
    gr_ECalEdep->SetTitle("");
    gr_ECalEdep->GetXaxis()->SetTitle("Beam Energy (GeV)");
    gr_ECalEdep->GetYaxis()->SetTitle("Mean ECal Edep (MeV)");
    gr_ECalEdep->SetMinimum(0.);
    gr_ECalEdep->SetMaximum(84.9);
    gr_ECalEdep->GetYaxis()->SetNdivisions(310);
    gr_ECalEdep->GetXaxis()->SetLimits(0, 11.);
    gr_ECalEdep->GetXaxis()->SetNdivisions(308);
    TString info_text;
    info_text = Form("ECal before cut");
    TLatex info_caption;
    info_caption.SetTextFont(62);
    info_caption.SetTextSize(.04);
    info_caption.SetNDC(kTRUE);
    info_caption.DrawLatex(.43, .84, info_text);
    gPad->SetTicks();


    mg_efficiency->Add(gr_efficiencyWholeAll);
    mg_efficiency->Add(gr_efficiencyWholeEO1232);
    mg_efficiency->Add(gr_efficiencyWholeAll1220);
    mg_efficiency->Add(gr_efficiencyWholeEO1220);
    mg_efficiency->Add(gr_efficiencywithoutDead);
    
    


    TCanvas* c_efficiency = new TCanvas("c_efficiency", "", 1200, 800);
    mg_efficiency->Draw("ap");
    mg_efficiency->SetTitle("");
    mg_efficiency->GetXaxis()->SetTitle("Beam Energy (GeV)");
    mg_efficiency->GetYaxis()->SetTitle("Efficiency (%)");
    mg_efficiency->SetMinimum(60.);
    mg_efficiency->SetMaximum(102.5);
    mg_efficiency->GetYaxis()->SetNdivisions(310);
    mg_efficiency->GetXaxis()->SetLimits(0, 11.);
    mg_efficiency->GetXaxis()->SetNdivisions(308);

    gr_efficiencywithoutDead->SetMarkerStyle(kFullCircle);
    gr_efficiencywithoutDead->SetMarkerColor(kBlack);

    gr_efficiencyWholeAll->SetMarkerStyle(kOpenCircle);
    gr_efficiencyWholeAll->SetMarkerColor(kCyan);

    gr_efficiencyWholeAll1220->SetMarkerStyle(kOpenCircle);
    gr_efficiencyWholeAll1220->SetMarkerColor(kBlue);

    gr_efficiencyWholeEO1232->SetMarkerStyle(kOpenCircle);
    gr_efficiencyWholeEO1232->SetMarkerColor(kRed);

    gr_efficiencyWholeEO1220->SetMarkerStyle(kOpenCircle);
    gr_efficiencyWholeEO1220->SetMarkerColor(kGreen+2);

    info_text = Form("ECal + HCal: ECal 1 MeV Cut, HCal 0.1 MeV Cut");
    info_caption.DrawLatex(.26, .14, info_text);
    gPad->SetTicks();


    TLegend* leg_efficiencies = new TLegend(0.45, 0.45, 0.65, 0.65);
    leg_efficiencies->AddEntry(gr_efficiencywithoutDead, "No dead layers", "p");
    leg_efficiencies->AddEntry(gr_efficiencyWholeAll, "HCal Layers #in [1,39]", "p");
    leg_efficiencies->AddEntry(gr_efficiencyWholeAll1220, "HCal Layers #in [1,9]", "p");
    leg_efficiencies->AddEntry(gr_efficiencyWholeEO1232, "Every Other HCal Layer #in [1, 21]", "p");
    leg_efficiencies->AddEntry(gr_efficiencyWholeEO1220, "Every Other HCal Layer #in [1, 9]", "p");
    leg_efficiencies->SetTextFont(43);
    leg_efficiencies->SetTextSize(29.76);
    leg_efficiencies->SetBorderSize(0);
    leg_efficiencies->Draw("same");

    TMultiGraph* mg_ECalHCalefficiency = new TMultiGraph();


    TGraph* gr_efficiencyHCalEO1220 = new TGraph(10, energy, efficiencyHCalEO1220.data());
    TGraph* gr_efficiencyECal = new TGraph(10, energy, efficiencyECal.data());

    mg_ECalHCalefficiency->Add(gr_efficiencyHCalEO1220);
    mg_ECalHCalefficiency->Add(gr_efficiencyECal);

    TCanvas* c_ECalHCalefficiency = new TCanvas("c_ECalHCalefficiency", "", 1200, 800);
    mg_ECalHCalefficiency->Draw("ap");
    mg_ECalHCalefficiency->SetTitle("");
    mg_ECalHCalefficiency->GetXaxis()->SetTitle("Beam Energy (GeV)");
    mg_ECalHCalefficiency->GetYaxis()->SetTitle("Efficiency (%)");
    mg_ECalHCalefficiency->SetMinimum(40.);
    mg_ECalHCalefficiency->SetMaximum(105);
    mg_ECalHCalefficiency->GetYaxis()->SetNdivisions(310);
    mg_ECalHCalefficiency->GetXaxis()->SetLimits(0, 11.);
    mg_ECalHCalefficiency->GetXaxis()->SetNdivisions(308);

    gr_efficiencyHCalEO1220->SetMarkerStyle(kOpenCircle);
    gr_efficiencyHCalEO1220->SetMarkerColor(kRed);

    gr_efficiencyECal->SetMarkerStyle(kOpenCircle);
    gr_efficiencyECal->SetMarkerColor(kBlue);
    
    info_text = Form("ECal 1 MeV Cut, HCal 0.1 MeV Cut");
    info_caption.DrawLatex(.32, .14, info_text);
    gPad->SetTicks();

    TLegend* leg_ECalHCalefficiency = new TLegend(0.15, 0.75, 0.25, 0.85);
    leg_ECalHCalefficiency->AddEntry(gr_efficiencyHCalEO1220, "HCal, Every Other Layer #in [1,9]", "p");
    leg_ECalHCalefficiency->AddEntry(gr_efficiencyECal, "ECal", "p");
    leg_ECalHCalefficiency->SetTextFont(43);
    leg_ECalHCalefficiency->SetTextSize(29.76);
    leg_ECalHCalefficiency->SetBorderSize(0);
    leg_ECalHCalefficiency->Draw("same");



    TCanvas* c_edep = new TCanvas("c_edep", "", 1200, 800);

    neutron1GeV_edep->Draw("EX0");
    neutron1GeV_edep->SetMarkerColor(kRed);
    neutron1GeV_edep->SetAxisRange(-0.2, 5., "y");

    neutron2GeV_edep->Draw("EX0 same");
    neutron2GeV_edep->SetMarkerColor(kBlue);
    
    neutron3GeV_edep->Draw("EX0 same");
    neutron3GeV_edep->SetMarkerColor(kGreen);

    neutron4GeV_edep->Draw("EX0 same");
    neutron4GeV_edep->SetMarkerColor(kOrange);

    neutron5GeV_edep->Draw("EX0 same");
    neutron5GeV_edep->SetMarkerColor(kCyan);

    neutron6GeV_edep->Draw("EX0 same");
    neutron6GeV_edep->SetMarkerColor(kMagenta);

    neutron7GeV_edep->Draw("EX0 same");
    neutron7GeV_edep->SetMarkerColor(kGreen+2);

    neutron8GeV_edep->Draw("EX0 same");
    neutron8GeV_edep->SetMarkerColor(kBlack);

    neutron9GeV_edep->Draw("EX0 same");
    neutron9GeV_edep->SetMarkerColor(kGray);

    neutron10GeV_edep->Draw("EX0 same");
    neutron10GeV_edep->SetMarkerColor(kViolet+1);

    gStyle->SetOptStat(0);

    TLine* l0 = new TLine(0.5, 0, 40.5, 0);
    l0->SetLineStyle(2);
    l0->Draw("same");

    TLegend* leg_energies = new TLegend(0.65, 0.45, 0.85, 0.85);
    leg_energies->AddEntry(neutron1GeV_edep, "1 GeV", "p");
    leg_energies->AddEntry(neutron2GeV_edep, "2 GeV", "p");
    leg_energies->AddEntry(neutron3GeV_edep, "3 GeV", "p");
    leg_energies->AddEntry(neutron4GeV_edep, "4 GeV", "p");
    leg_energies->AddEntry(neutron5GeV_edep, "5 GeV", "p");
    leg_energies->AddEntry(neutron6GeV_edep, "6 GeV", "p");
    leg_energies->AddEntry(neutron7GeV_edep, "7 GeV", "p");
    leg_energies->AddEntry(neutron8GeV_edep, "8 GeV", "p");
    leg_energies->AddEntry(neutron9GeV_edep, "9 GeV", "p");
    leg_energies->AddEntry(neutron10GeV_edep, "10 GeV", "p");
    leg_energies->SetTextFont(43);
    leg_energies->SetTextSize(29.76);

    leg_energies->SetBorderSize(0);
    leg_energies->Draw("same");


    TCanvas* c_TileHits = new TCanvas("c_TileHits", "", 1200, 800);

    neutron1GeV_TileHits->Draw("EX0");
    neutron1GeV_TileHits->SetMarkerColor(kRed);
    neutron1GeV_TileHits->SetAxisRange(-0.1, 3, "y");

    neutron2GeV_TileHits->Draw("EX0 same");
    neutron2GeV_TileHits->SetMarkerColor(kBlue);
    
    neutron3GeV_TileHits->Draw("EX0 same");
    neutron3GeV_TileHits->SetMarkerColor(kGreen);

    neutron4GeV_TileHits->Draw("EX0 same");
    neutron4GeV_TileHits->SetMarkerColor(kOrange);

    neutron5GeV_TileHits->Draw("EX0 same");
    neutron5GeV_TileHits->SetMarkerColor(kCyan);

    neutron6GeV_TileHits->Draw("EX0 same");
    neutron6GeV_TileHits->SetMarkerColor(kMagenta);

    neutron7GeV_TileHits->Draw("EX0 same");
    neutron7GeV_TileHits->SetMarkerColor(kGreen+2);

    neutron8GeV_TileHits->Draw("EX0 same");
    neutron8GeV_TileHits->SetMarkerColor(kBlack);

    neutron9GeV_TileHits->Draw("EX0 same");
    neutron9GeV_TileHits->SetMarkerColor(kGray);

    neutron10GeV_TileHits->Draw("EX0 same");
    neutron10GeV_TileHits->SetMarkerColor(kViolet+1);

    gStyle->SetOptStat(0);

    l0->Draw("same");

    TLegend* leg_TileHits = new TLegend(0.65, 0.45, 0.85, 0.85);
    leg_TileHits->AddEntry(neutron1GeV_TileHits, "1 GeV", "p");
    leg_TileHits->AddEntry(neutron2GeV_TileHits, "2 GeV", "p");
    leg_TileHits->AddEntry(neutron3GeV_TileHits, "3 GeV", "p");
    leg_TileHits->AddEntry(neutron4GeV_TileHits, "4 GeV", "p");
    leg_TileHits->AddEntry(neutron5GeV_TileHits, "5 GeV", "p");
    leg_TileHits->AddEntry(neutron6GeV_TileHits, "6 GeV", "p");
    leg_TileHits->AddEntry(neutron7GeV_TileHits, "7 GeV", "p");
    leg_TileHits->AddEntry(neutron8GeV_TileHits, "8 GeV", "p");
    leg_TileHits->AddEntry(neutron9GeV_TileHits, "9 GeV", "p");
    leg_TileHits->AddEntry(neutron10GeV_TileHits, "10 GeV", "p");
    leg_TileHits->SetTextFont(43);
    leg_TileHits->SetTextSize(29.76);
    leg_TileHits->SetBorderSize(0);
    leg_TileHits->Draw("same");


    TCanvas* c_numHits = new TCanvas("c_numHits", "", 1200, 800);

    neutron1GeV_numHits->Draw("EX0");
    neutron1GeV_numHits->SetMarkerColor(kRed);
    neutron1GeV_numHits->SetAxisRange(-1.9, 60, "y");

    neutron2GeV_numHits->Draw("EX0 same");
    neutron2GeV_numHits->SetMarkerColor(kBlue);
    
    neutron3GeV_numHits->Draw("EX0 same");
    neutron3GeV_numHits->SetMarkerColor(kGreen);

    neutron4GeV_numHits->Draw("EX0 same");
    neutron4GeV_numHits->SetMarkerColor(kOrange);

    neutron5GeV_numHits->Draw("EX0 same");
    neutron5GeV_numHits->SetMarkerColor(kCyan);

    neutron6GeV_numHits->Draw("EX0 same");
    neutron6GeV_numHits->SetMarkerColor(kMagenta);

    neutron7GeV_numHits->Draw("EX0 same");
    neutron7GeV_numHits->SetMarkerColor(kGreen+2);

    neutron8GeV_numHits->Draw("EX0 same");
    neutron8GeV_numHits->SetMarkerColor(kBlack);

    neutron9GeV_numHits->Draw("EX0 same");
    neutron9GeV_numHits->SetMarkerColor(kGray);

    neutron10GeV_numHits->Draw("EX0 same");
    neutron10GeV_numHits->SetMarkerColor(kViolet+1);

    gStyle->SetOptStat(0);

    l0->Draw("same");

    TLegend* leg_numHits = new TLegend(0.65, 0.45, 0.85, 0.85);
    leg_numHits->AddEntry(neutron1GeV_numHits, "1 GeV", "p");
    leg_numHits->AddEntry(neutron2GeV_numHits, "2 GeV", "p");
    leg_numHits->AddEntry(neutron3GeV_numHits, "3 GeV", "p");
    leg_numHits->AddEntry(neutron4GeV_numHits, "4 GeV", "p");
    leg_numHits->AddEntry(neutron5GeV_numHits, "5 GeV", "p");
    leg_numHits->AddEntry(neutron6GeV_numHits, "6 GeV", "p");
    leg_numHits->AddEntry(neutron7GeV_numHits, "7 GeV", "p");
    leg_numHits->AddEntry(neutron8GeV_numHits, "8 GeV", "p");
    leg_numHits->AddEntry(neutron9GeV_numHits, "9 GeV", "p");
    leg_numHits->AddEntry(neutron10GeV_numHits, "10 GeV", "p");
    leg_numHits->SetTextFont(43);
    leg_numHits->SetTextSize(29.76);
    leg_numHits->SetBorderSize(0);
    leg_numHits->Draw("same");



    TCanvas* c_TileHitsHists = new TCanvas("c_TileHitsHists", "", 1200, 800);
    c_TileHitsHists->SetLogy();
    h_neutron1GeV_TileHits->Draw();
    h_neutron1GeV_TileHits->SetLineColor(kRed);

    h_neutron1GeV_TileHits->GetXaxis()->SetTitle("Number of tiles hit, every other layer #in [1,9]");
    h_neutron1GeV_TileHits->GetYaxis()->SetTitle("Number of events (Log)");
    h_neutron1GeV_TileHits->SetTitle("");
    h_neutron1GeV_TileHits->SetAxisRange(0,50, "x");
    

    h_neutron2GeV_TileHits->Draw("same");
    h_neutron2GeV_TileHits->SetLineColor(kBlue);
    
    h_neutron3GeV_TileHits->Draw("same");
    h_neutron3GeV_TileHits->SetLineColor(kGreen);

    h_neutron4GeV_TileHits->Draw("same");
    h_neutron4GeV_TileHits->SetLineColor(kOrange);

    h_neutron5GeV_TileHits->Draw("same");
    h_neutron5GeV_TileHits->SetLineColor(kCyan);

    h_neutron6GeV_TileHits->Draw("same");
    h_neutron6GeV_TileHits->SetLineColor(kMagenta);

    h_neutron7GeV_TileHits->Draw("same");
    h_neutron7GeV_TileHits->SetLineColor(kGreen+2);

    h_neutron8GeV_TileHits->Draw("same");
    h_neutron8GeV_TileHits->SetLineColor(kBlack);

    h_neutron9GeV_TileHits->Draw("same");
    h_neutron9GeV_TileHits->SetLineColor(kGray);

    h_neutron10GeV_TileHits->Draw("same");
    h_neutron10GeV_TileHits->SetLineColor(kViolet+1);

    TLegend* leg_TileHitsHists = new TLegend(0.65, 0.45, 0.85, 0.85);
    leg_TileHitsHists->AddEntry(h_neutron1GeV_TileHits, "1 GeV", "l");
    leg_TileHitsHists->AddEntry(h_neutron2GeV_TileHits, "2 GeV", "l");
    leg_TileHitsHists->AddEntry(h_neutron3GeV_TileHits, "3 GeV", "l");
    leg_TileHitsHists->AddEntry(h_neutron4GeV_TileHits, "4 GeV", "l");
    leg_TileHitsHists->AddEntry(h_neutron5GeV_TileHits, "5 GeV", "l");
    leg_TileHitsHists->AddEntry(h_neutron6GeV_TileHits, "6 GeV", "l");
    leg_TileHitsHists->AddEntry(h_neutron7GeV_TileHits, "7 GeV", "l");
    leg_TileHitsHists->AddEntry(h_neutron8GeV_TileHits, "8 GeV", "l");
    leg_TileHitsHists->AddEntry(h_neutron9GeV_TileHits, "9 GeV", "l");
    leg_TileHitsHists->AddEntry(h_neutron10GeV_TileHits, "10 GeV", "l");
    leg_TileHitsHists->SetTextFont(43);
    leg_TileHitsHists->SetTextSize(29.76);
    leg_TileHitsHists->SetBorderSize(0);
    leg_TileHitsHists->Draw("same");


}