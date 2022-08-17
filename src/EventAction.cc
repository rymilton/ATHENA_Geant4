
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "CalorimeterSD.hh"
#include "CalorHit.hh"
#include "Analysis.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "GlobalValues.hh"

using namespace GlobalValues;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
 : G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorHitsCollection* 
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<CalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::PrintEventStatistics(
                              G4double ECalEdep, G4double HCalEdep) const
{
  // print event statistics
  G4cout
     << "   ECal: total energy: " 
     << std::setw(7) << G4BestUnit(ECalEdep, "Energy")
     << G4endl
     << "        HCal: total energy: " 
     << std::setw(7) << G4BestUnit(HCalEdep, "Energy")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{  
  auto eventID = event->GetEventID();

  char nameHolder[200];
  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Getting and reading out HCal event data

  G4double hcal_active_edep = 0.; // Total Edep in HCal scintillators
  G4double hcal_active_edepPi0 = 0.; // Total Edep from pi0 in HCal scint.
  G4int hcal_active_num_Pi0 = 0; // Number of pi0 in HCal scint.

  G4double hcal_absorber_edep = 0.; // Total Edep in HCal absorbers (including WLS plates and steel plates)
  G4double hcal_absorber_edepPi0 = 0.; // Total Edep from pi0 in HCal absorbers
  G4int hcal_absorber_num_Pi0 = 0; // Number of pi0 in HCal absorbers

  // Getting HCal information.

  for(G4int i = 0; i < NumHCalTowers; i++)
  {
    for(G4int j = 0; j < NumHCalTowers; j++)
    {
      G4double hcal_active_tower_edep = 0.;
      G4double hcal_active_tower_edepPi0 = 0.;
      G4double hcal_absorber_tower_edep = 0.;
      G4double hcal_absorber_tower_edepPi0 = 0.;
      G4int layer_tracker = 0; // Tracks which layer the tiles are in
      G4int hcal_active_tower_numPi0 = 0;
      G4int hcal_absorber_tower_numPi0 = 0;

      // Getting scintillator info
      sprintf(nameHolder, "HCal_ActiveHitsCollection%d%d", i, j);
      G4int HCal_ActiveHCID = G4SDManager::GetSDMpointer()->GetCollectionID(nameHolder);
      auto HCal_ActiveHC = GetHitsCollection(HCal_ActiveHCID, event);
      auto HCal_ActiveHit = (*HCal_ActiveHC)[HCal_ActiveHC->entries()-1]; // entries()-1 kept track of information for whole tower segmentation
      hcal_active_edep += HCal_ActiveHit->GetEdep();
      hcal_active_edepPi0 += HCal_ActiveHit->GetEdepPi0();
      hcal_active_num_Pi0 += HCal_ActiveHit->GetNumPi0();
      hcal_active_tower_edep += HCal_ActiveHit->GetEdep();
      hcal_active_tower_edepPi0 += HCal_ActiveHit->GetEdepPi0();
      hcal_active_tower_numPi0 += HCal_ActiveHit->GetNumPi0();
      
      // Getting absorber info (just absorbers in towers here)
      sprintf(nameHolder, "HCal_AbsorberHitsCollection%d%d", i, j);
      G4int HCal_AbsorberHCID = G4SDManager::GetSDMpointer()->GetCollectionID(nameHolder);
      auto HCal_AbsorberHC = GetHitsCollection(HCal_AbsorberHCID, event);
      auto HCal_AbsorberHit = (*HCal_AbsorberHC)[HCal_AbsorberHC->entries()-1]; // entries()-1 kept track of information for whole tower segmentation
      hcal_absorber_edep += HCal_AbsorberHit->GetEdep();
      hcal_absorber_edepPi0 += HCal_AbsorberHit->GetEdepPi0();
      hcal_absorber_num_Pi0 += HCal_AbsorberHit->GetNumPi0();
      hcal_absorber_tower_edep += HCal_AbsorberHit->GetEdep();
      hcal_absorber_tower_edepPi0 += HCal_AbsorberHit->GetEdepPi0();
      hcal_absorber_tower_numPi0 += HCal_AbsorberHit->GetNumPi0();

      // Looping over HCal layers. Getting individual tile information
      for(G4int k = 0; k < NumHCalLayers; k++)
      { 
        // Ntuple with id 3 holds HCal tile information
        auto HCal_ActiveTileHit = (*HCal_ActiveHC)[k]; // Tile is each of scintillating plates in the HCal towers
        auto HCal_AbsorberTileHit = (*HCal_AbsorberHC)[k]; // Individual absorber in the HCal towers

        analysisManager->FillNtupleDColumn(3, 0,  HCal_ActiveTileHit->GetEdep());
        analysisManager->FillNtupleDColumn(3, 1,  HCal_ActiveTileHit->GetEdepPi0());
        analysisManager->FillNtupleDColumn(3, 2,  HCal_AbsorberTileHit->GetEdep());
        analysisManager->FillNtupleDColumn(3, 3,  HCal_AbsorberTileHit->GetEdepPi0());
        analysisManager->FillNtupleIColumn(3, 4, layer_tracker);
        analysisManager->FillNtupleIColumn(3, 5,  HCal_ActiveTileHit->GetNumPi0()); 
        analysisManager->FillNtupleIColumn(3, 6,  HCal_AbsorberTileHit->GetNumPi0()); 
        analysisManager->FillNtupleIColumn(3, 7, i);
        analysisManager->FillNtupleIColumn(3, 8, j);
        analysisManager->FillNtupleIColumn(3, 9, eventID);
        analysisManager->AddNtupleRow(3);
        layer_tracker++;
      }
      
      // Ntuple with id 2 holds HCal tower information
      analysisManager->FillNtupleDColumn(2, 0, hcal_active_tower_edep);
      analysisManager->FillNtupleDColumn(2, 1, hcal_active_tower_edepPi0);
      analysisManager->FillNtupleDColumn(2, 2, hcal_absorber_tower_edep);
      analysisManager->FillNtupleDColumn(2, 3, hcal_absorber_tower_edepPi0);
      analysisManager->FillNtupleIColumn(2, 4, i);
      analysisManager->FillNtupleIColumn(2, 5, j);
      analysisManager->FillNtupleIColumn(2, 6, hcal_active_tower_numPi0); 
      analysisManager->FillNtupleIColumn(2, 7, hcal_absorber_tower_numPi0); 
      analysisManager->FillNtupleIColumn(2, 8, eventID);
      analysisManager->AddNtupleRow(2);
    }
  }
  // Info from the steel plates and WLS plates in the HCal
  // Combining this info with absorber info (i.e. non-scintillating materials)
  G4int HCal_PlatesHCID = G4SDManager::GetSDMpointer()->GetCollectionID("HCal_PlatesHitCollection");
  auto HCal_PlatesHC = GetHitsCollection(HCal_PlatesHCID, event);
  auto HCal_PlatesHit = (*HCal_PlatesHC)[HCal_PlatesHC->entries()-1]; // entries()-1 kept track of information for whole tower segmentation
  hcal_absorber_edep += HCal_PlatesHit->GetEdep();
  hcal_absorber_edepPi0 += HCal_PlatesHit->GetEdepPi0();
  hcal_absorber_num_Pi0 += HCal_PlatesHit->GetNumPi0();


  // Getting and reading out ECal event data

  G4double ecal_fiber_active_edep = 0.; // Total Edep for ECal fiber cores
  G4double ecal_fiber_active_edepPi0 = 0.; // Total pi0 edep for ECal fiber cores
  G4int ecal_fiber_active_num_Pi0 = 0; // Total number pi0 hits for ECal fiber cores

  G4double ecal_absorber_edep = 0.; // Total Edep for ECal Absorber (including cladding & glue)
  G4double ecal_absorber_edepPi0 = 0.; // Total pi0 edep for ECal absorber
  G4int ecal_absorber_num_Pi0 = 0; // Total number pi0 hits for ECal fiber cores

  for(G4int i = 0; i < NumECalBlocks; i++)
  {
    for(G4int j = 0; j < NumECalBlocks; j++)
    {
      // Fiber core info
      sprintf(nameHolder, "ECal_FiberHitsCollection%d%d", i, j);
      G4int ECal_FiberHCID = G4SDManager::GetSDMpointer()->GetCollectionID(nameHolder);
      auto ECal_FiberHC = GetHitsCollection(ECal_FiberHCID, event);
      auto ECal_FiberHit = (*ECal_FiberHC)[ECal_FiberHC->entries()-1]; // entries()-1 kept track of information for whole block
      ecal_fiber_active_edep += ECal_FiberHit->GetEdep();
      ecal_fiber_active_edepPi0 += ECal_FiberHit->GetEdepPi0();
      ecal_fiber_active_num_Pi0 += ECal_FiberHit->GetNumPi0();

      // Absorber info (only tungsten powder and cladding here)
      sprintf(nameHolder, "ECal_AbsorberHitsCollection%d%d", i, j);
      G4int ECal_AbsHCID = G4SDManager::GetSDMpointer()->GetCollectionID(nameHolder);
      auto ECal_AbsHC = GetHitsCollection(ECal_AbsHCID, event);
      auto ECal_AbsHit = (*ECal_AbsHC)[ECal_AbsHC->entries()-1]; // entries()-1 kept track of information for whole block
      ecal_absorber_edep += ECal_AbsHit->GetEdep();
      ecal_absorber_edepPi0 += ECal_AbsHit->GetEdepPi0();
      ecal_absorber_num_Pi0 += ECal_AbsHit->GetNumPi0();

      // Ntuple with id 1 holds ECal information
      analysisManager->FillNtupleDColumn(1, 0, ECal_FiberHit->GetEdep());
      analysisManager->FillNtupleDColumn(1, 1, ECal_FiberHit->GetEdepPi0());
      analysisManager->FillNtupleDColumn(1, 2, ECal_AbsHit->GetEdep());
      analysisManager->FillNtupleDColumn(1, 3, ECal_AbsHit->GetEdepPi0());
      analysisManager->FillNtupleIColumn(1, 4, ECal_FiberHit->GetNumPi0());
      analysisManager->FillNtupleIColumn(1, 5, ECal_AbsHit->GetNumPi0());
      analysisManager->FillNtupleIColumn(1, 6, i);
      analysisManager->FillNtupleIColumn(1, 7, j);
      analysisManager->FillNtupleIColumn(1, 8, eventID);
      analysisManager->AddNtupleRow(1);
    }
  }

  // Info from the glue in the HCal
  // Combining this info with absorber info (i.e. non-scintillating materials) 
  G4int ECal_GlueHCID = G4SDManager::GetSDMpointer()->GetCollectionID("ECal_GlueHitCollection");
  auto ECal_GlueHC = GetHitsCollection(ECal_GlueHCID, event);
  auto ECal_GlueHit = (*ECal_GlueHC)[ECal_GlueHC->entries()-1]; // entries()-1 kept track of information for whole tower segmentation
  ecal_absorber_edep += ECal_GlueHit->GetEdep();
  ecal_absorber_edepPi0 += ECal_GlueHit->GetEdepPi0();
  ecal_absorber_num_Pi0 += ECal_GlueHit->GetNumPi0();


  // Ntuple with id 0 holds info for entire detector
  analysisManager->FillNtupleDColumn(0, 0, ecal_fiber_active_edep);
  analysisManager->FillNtupleDColumn(0, 1, ecal_fiber_active_edepPi0);
  analysisManager->FillNtupleDColumn(0, 2, hcal_active_edep);
  analysisManager->FillNtupleDColumn(0, 3, hcal_active_edepPi0);
  analysisManager->FillNtupleDColumn(0, 4, ecal_absorber_edep);
  analysisManager->FillNtupleDColumn(0, 5, ecal_absorber_edepPi0);
  analysisManager->FillNtupleDColumn(0, 6, hcal_absorber_edep);
  analysisManager->FillNtupleDColumn(0, 7, hcal_absorber_edepPi0);
  analysisManager->FillNtupleIColumn(0, 8, ecal_fiber_active_num_Pi0);
  analysisManager->FillNtupleIColumn(0, 9, hcal_active_num_Pi0);
  analysisManager->FillNtupleIColumn(0, 10, ecal_absorber_num_Pi0);
  analysisManager->FillNtupleIColumn(0, 11, hcal_absorber_num_Pi0);
  analysisManager->FillNtupleIColumn(0, 12, eventID);
  analysisManager->AddNtupleRow(0); 

  
  if(eventID % 1000 == 0) G4cout << "---> End of event: " << eventID << G4endl; 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
