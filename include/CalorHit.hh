
/// \file CalorHit.hh
/// \brief Definition of the CalorHit class

#ifndef CalorHit_h
#define CalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class CalorHit : public G4VHit
{
  public:
    CalorHit();
    CalorHit(const CalorHit&);
    virtual ~CalorHit();

    // operators
    const CalorHit& operator=(const CalorHit&);
    G4bool operator==(const CalorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double de, G4double dl, G4double dePi0, G4int nPi0);

    // get methods
    G4double GetEdep() const;
    G4double GetTrackLength() const;
    G4double GetEdepPi0() const;
    G4int GetNumPi0() const;

  private:
    G4double fEdep;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength; ///< Track length in the  sensitive volume
    G4double fEdepPi0;
    G4int fNumPi0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using CalorHitsCollection = G4THitsCollection<CalorHit>;

extern G4ThreadLocal G4Allocator<CalorHit>* CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* CalorHit::operator new(size_t)
{
  if (!CalorHitAllocator) {
    CalorHitAllocator = new G4Allocator<CalorHit>;
  }
  void *hit;
  hit = (void *) CalorHitAllocator->MallocSingle();
  return hit;
}

inline void CalorHit::operator delete(void *hit)
{
  if (!CalorHitAllocator) {
    CalorHitAllocator = new G4Allocator<CalorHit>;
  }
  CalorHitAllocator->FreeSingle((CalorHit*) hit);
}

inline void CalorHit::Add(G4double de, G4double dl, G4double dePi0, G4int nPi0) {
  fEdep += de; 
  fTrackLength += dl;
  fEdepPi0 += dePi0;
  fNumPi0 += nPi0;
}

inline G4double CalorHit::GetEdep() const { 
  return fEdep; 
}

inline G4double CalorHit::GetTrackLength() const { 
  return fTrackLength; 
}
inline G4double CalorHit::GetEdepPi0() const {
  return fEdepPi0;
}
inline G4int CalorHit::GetNumPi0() const {
  return fNumPi0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
