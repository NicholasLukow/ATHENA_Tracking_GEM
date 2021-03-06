// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICFRICH_STEPPINGACTION_H
#define EICFRICH_STEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class dRICH_Detector;

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;

class dRICH_SteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  dRICH_SteppingAction(dRICH_Detector*, const PHParameters* parameters);

  //! destructor
  virtual ~dRICH_SteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  dRICH_Detector* m_Detector;
  const PHParameters* m_Params;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainer;
  PHG4Hit* m_Hit;
  PHG4HitContainer* m_SaveHitContainer;
  G4VPhysicalVolume* m_SaveVolPre;
  G4VPhysicalVolume* m_SaveVolPost;

  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  int m_ActiveFlag;
  int m_BlackHoleFlag;
  double m_EdepSum;
  double m_EionSum;
};

#endif // EICFRICH_STEPPINGACTION_H
