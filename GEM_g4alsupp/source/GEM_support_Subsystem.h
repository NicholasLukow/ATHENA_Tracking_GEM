// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GEM_AL_SUPPORT_SUBSYSTEM_H
#define GEM_AL_SUPPORT_SUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

class PHCompositeNode;
class PHG4Detector;
class GEM_support_Detector;
class PHG4SteppingAction;

/**
   * \brief Detector Subsystem module
   *
   * The detector is constructed and registered via GEM_support_Detector
   *
   *
   * \see GEM_support_Detector
   * \see GEM_support_Subsystem
   *
   */
class GEM_support_Subsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  GEM_support_Subsystem(const std::string& name = "GEM_support_");

  //! destructor
  virtual ~GEM_support_Subsystem() {}

  /*!
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*) override;

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*) override;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector() const override;

  PHG4SteppingAction* GetSteppingAction() const override { return m_SteppingAction; }
  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

 protected:
  // \brief Set default parameter values
  void SetDefaultParameters() override;

 private:
  //! detector construction
  /*! derives from PHG4Detector */
  GEM_support_Detector  *m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction *m_SteppingAction;
};

#endif // ALLSI_AL_SUPPORT_SUBSYSTEM_H
