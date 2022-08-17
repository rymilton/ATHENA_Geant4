
/// \file CalorimeterSD.cc
/// \brief Implementation of the CalorimeterSD class

#include "CalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "Analysis.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::CalorimeterSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName,
                            G4int nofCells)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr),
   fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::~CalorimeterSD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
    = new CalorHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  auto hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  // Create hits
  // fNofCells for cells + one more for total sums 
  for (G4int i=0; i<fNofCells+1; i++ ) {
    fHitsCollection->insert(new CalorHit());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CalorimeterSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory*)
{  
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  // step length
  // Using step length for Birk's formula, which applies only to charged particles
  G4double stepLength = (step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.) ? step->GetStepLength() : 0.;

  auto touchable = (step->GetPreStepPoint()->GetTouchable());
  
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  // Get calorimeter cell id
  auto layerNumber = touchable->GetReplicaNumber(1);

  // Get hit accounting data for this cell
  auto hit = (*fHitsCollection)[layerNumber];
  if ( ! hit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << layerNumber; 
    G4Exception("CalorimeterSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
  }
  
  // Get hit for total accounting
  auto hitTotal 
    = (*fHitsCollection)[fHitsCollection->entries()-1];

  // Adjusting the energy for the Birk's constant
  G4Material* mat = volume->GetLogicalVolume()->GetMaterial();
  G4double charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
  G4double birk = mat->GetIonisation()->GetBirksConstant();
  
  if (birk * edep * stepLength * charge != 0)
  {
    edep /= (1. + birk * edep / stepLength); // Done for charged particles in organic scintillators
  }

  // Tracking and saving pi0 information
  // Used during analysis to calculate the electromagnetic fraction
  G4double energyPi0 = 0.;
  G4int numPi0 = 0;
  if (step->GetTrack()->GetDefinition()->GetParticleName() == "pi0")
  {
    energyPi0 = step->GetTrack()->GetTotalEnergy();
    auto analysisManager = G4AnalysisManager::Instance();
    G4int event_number = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    analysisManager->FillNtupleDColumn(4, 0, energyPi0);
    analysisManager->FillNtupleDColumn(4, 1, step->GetPreStepPoint()->GetPosition().x()/cm);
    analysisManager->FillNtupleDColumn(4, 2, step->GetPreStepPoint()->GetPosition().y()/cm);
    analysisManager->FillNtupleDColumn(4, 3, step->GetPreStepPoint()->GetPosition().z()/cm + 8.5);
    analysisManager->FillNtupleIColumn(4, 4, event_number);
    analysisManager->AddNtupleRow(4);
    numPi0++;
  }


  // Add values
  hit->Add(edep, stepLength, energyPi0, numPi0);
  hitTotal->Add(edep, stepLength, energyPi0, numPi0);  
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     auto nofHits = fHitsCollection->entries();
     G4cout
       << G4endl 
       << "-------->Hits Collection: in this event they are " << nofHits 
       << " hits in the tracker chambers: " << G4endl;
     for ( std::size_t i=0; i<nofHits; ++i ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
