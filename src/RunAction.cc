
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(0);     

  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();

  // Create directories 
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetNtupleMerging(true);
  // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple

  analysisManager->CreateNtuple("EdepTotal", "EdepTotal");
  analysisManager->CreateNtupleDColumn("ECal_Edep_Active_Total");
  analysisManager->CreateNtupleDColumn("ECal_EdepPi0_Active_Total");
  analysisManager->CreateNtupleDColumn("HCal_Edep_Active_Total");
  analysisManager->CreateNtupleDColumn("HCal_EdepPi0_Active_Total");
  analysisManager->CreateNtupleDColumn("ECal_Edep_Absorber_Total");
  analysisManager->CreateNtupleDColumn("ECal_EdepPi0_Absorber_Total");
  analysisManager->CreateNtupleDColumn("HCal_Edep_Absorber_Total");
  analysisManager->CreateNtupleDColumn("HCal_EdepPi0_Absorber_Total");
  analysisManager->CreateNtupleIColumn("ECal_Num_Active_Pi0");
  analysisManager->CreateNtupleIColumn("HCal_Num_Active_Pi0");
  analysisManager->CreateNtupleIColumn("ECal_Num_Absorber_Pi0");
  analysisManager->CreateNtupleIColumn("HCal_Num_Absorber_Pi0");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("ECalBlocks", "ECalBlocks");
  analysisManager->CreateNtupleDColumn("ECal_Edep_Active_Block");
  analysisManager->CreateNtupleDColumn("ECal_EdepPi0_Active_Block");
  analysisManager->CreateNtupleDColumn("ECal_Edep_Absorber_Block");
  analysisManager->CreateNtupleDColumn("ECal_EdepPi0_Absorber_Block");
  analysisManager->CreateNtupleIColumn("ECal_Num_Active_Pi0");
  analysisManager->CreateNtupleIColumn("ECal_Num_Absorber_Pi0");
  analysisManager->CreateNtupleIColumn("ECal_BlockXid");
  analysisManager->CreateNtupleIColumn("ECal_BlockYid");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("HCalTowers", "HCalTowers");
  analysisManager->CreateNtupleDColumn("HCal_Edep_Active_Tower");
  analysisManager->CreateNtupleDColumn("HCal_EdepPi0_Active_Tower");
  analysisManager->CreateNtupleDColumn("HCal_Edep_Absorber_Tower");
  analysisManager->CreateNtupleDColumn("HCal_EdepPi0_Absorber_Tower");
  analysisManager->CreateNtupleIColumn("HCal_TowerXid");
  analysisManager->CreateNtupleIColumn("HCal_TowerYid");
  analysisManager->CreateNtupleIColumn("HCal_Num_Active_Pi0");
  analysisManager->CreateNtupleIColumn("HCal_Num_Absorber_Pi0");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("HCalTiles", "HCalTiles");
  analysisManager->CreateNtupleDColumn("HCal_Edep_Active_Tile");
  analysisManager->CreateNtupleDColumn("HCal_EdepPi0_Active_Tile");
  analysisManager->CreateNtupleDColumn("HCal_Edep_Absorber_Tile");
  analysisManager->CreateNtupleDColumn("HCal_EdepPi0_Absorber_Tile");
  analysisManager->CreateNtupleIColumn("HCal_Layerid");
  analysisManager->CreateNtupleIColumn("HCal_NumPi0_Active_Tile");
  analysisManager->CreateNtupleIColumn("HCal_NumPi0_Absorber_Tile");
  analysisManager->CreateNtupleIColumn("HCal_TowerXid");
  analysisManager->CreateNtupleIColumn("HCal_TowerYid");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("Pi0", "Pi0");
  analysisManager->CreateNtupleDColumn("Energy");
  analysisManager->CreateNtupleDColumn("PosX");
  analysisManager->CreateNtupleDColumn("PosY");
  analysisManager->CreateNtupleDColumn("PosZ");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  analysisManager->OpenFile(analysisManager->GetFileName()); // File name set via macro
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  auto analysisManager = G4AnalysisManager::Instance();

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
