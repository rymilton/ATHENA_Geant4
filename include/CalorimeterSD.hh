
/// \file CalorimeterSD.hh
/// \brief Definition of the CalorimeterSD class

#ifndef CalorimeterSD_h
#define CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"

#include "CalorHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

/// Calorimeter sensitive detector class
///
/// In Initialize(), it creates one hit for each calorimeter layer and one more
/// hit for accounting the total quantities in all layers.
///
/// The values are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class CalorimeterSD : public G4VSensitiveDetector
{
  public:
    CalorimeterSD(const G4String& name, 
                     const G4String& hitsCollectionName, 
                     G4int nofCells);
    virtual ~CalorimeterSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    CalorHitsCollection* fHitsCollection;
    G4int  fNofCells;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

