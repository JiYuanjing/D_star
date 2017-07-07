#ifndef StPicoDstarAnaHists__h
#define StPicoDstarAnaHists__h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoDstarEvent
 *  simultaneously and do analysis.
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "StAnaCuts.h"
#include "THnSparse.h"
#include "THn.h"

class TH1F;
class TH2F;
class TH3F;
class TFile;
class TString;
class StPicoPrescales;
class StPicoEvent;
class StPicoTrack;
class StKaonPion;
class TNtuple;
class StD0Pion;

class StPicoDstarAnaHists
{
  public:
   StPicoDstarAnaHists(TString fileBaseName,bool fillQaHists=true, bool fillBackgroundTrees=false);
   virtual ~StPicoDstarAnaHists();
   void addEvent(StPicoEvent const *);
   void addEventBeforeCut(StPicoEvent const *);
   void addCent(const double refmultCor, int centrality, const double reweight, const float vz);
   void addKaonPion(StKaonPion const*, bool unlike, bool tpc, bool tof, int centrality, const double reweight);
   //
//   void addD0Pion(StD0Pion const*, bool unlike, bool tpc, bool tof, int centrality, const double reweight);
void addD0SoftPion(StD0Pion const* const d0p, StKaonPion const* const kp, bool unlike,  int centrality, const double reweight);
void addSideBandBackground(StD0Pion const* const d0p, StKaonPion const* const k,bool unlike, int centrality, const double reweight);

//
   void addBackground(StKaonPion const*, StPicoTrack const* kaon, StPicoTrack const* pion, int ptBin, bool SB);
   void addTpcDenom1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX);
   void addHFTNumer1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX);
   void addQaNtuple(int, float, float, float, float, float, int, const double, float, int, int);
   void addDcaPtCent(float dca, float dcaXy, float  dcaZ, bool IsPion, bool IsKaon, bool IsProton, float pt,  int centrality, float Eta, float Phi, float Vz, float ZdcX);
   int getEtaIndexDca(float Eta) ;
   int getPhiIndexDca(float Phi) ;
   int getVzIndexDca(float Vz) ;
   // int getZdcxIndex(float ZdcX) ;
   int getEtaIndexRatio(float Eta) ;
   int getPhiIndexRatio(float Phi) ;
   int getVzIndexRatio(float Vz) ;
   void closeFile();
private:
   StPicoDstarAnaHists() {}

   bool mFillQaHists;
   bool mFillBackgroundTrees;
   StPicoPrescales* mPrescales;
   TFile* mOutFile;
   TH1F* mh1TotalEventsInRun;
   TH1F* mh1TotalEventsInRunBeforeCut;
   TH2F* mh2InvariantMassVsPt;
   TH2F* mh2InvariantMassVsPtLike;
   TH2F* mh2InvariantMassVsPtTof;
   TH2F* mh2InvariantMassVsPtTofLike;
   //centrality
   TH1F* mh1Cent;
   TH1F* mh1CentWg;
   TH1F* mh1gRefmultCor;
   TH1F* mh1gRefmultCorWg;
   TH2F* mh2CentVz;
   TH2F* mh2CentVzWg;
   TH3F* mh3InvariantMassVsPtVsCent;
   TH3F* mh3InvariantMassVsPtVsCentLike;
   TH3F* mh3InvariantMassVsPtVsCentTof;
   TH3F* mh3InvariantMassVsPtVsCentTofLike;
   //dstar
   TH3F* mh3dstarInvariantMassVsPtVsCent;
   TH3F* mh3dstarInvariantMassVsPtVsCentTof;
   //dstar bar
   TH3F* mh3dstarbarInvariantMassVsPtVsCent;
   TH3F* mh3dstarbarInvariantMassVsPtVsCentTof;
   //LikeSign Pos+Neg
   TH3F* mh3InvariantMassVsPtVsCentLikePos;
   TH3F* mh3InvariantMassVsPtVsCentTofLikePos;
   //LikeSign Pos+Neg
   TH3F* mh3InvariantMassVsPtVsCentLikeNeg;
   TH3F* mh3InvariantMassVsPtVsCentTofLikeNeg;
   //HFT Left
   TH3F* mh3InvariantMassVsPtVsCentLeft;
   TH3F* mh3InvariantMassVsPtVsCentLikeLeft;
   TH3F* mh3InvariantMassVsPtVsCentTofLeft;
   TH3F* mh3InvariantMassVsPtVsCentTofLikeLeft;
   //HFT Right
   TH3F* mh3InvariantMassVsPtVsCentRight;
   TH3F* mh3InvariantMassVsPtVsCentLikeRight;
   TH3F* mh3InvariantMassVsPtVsCentTofRight;
   TH3F* mh3InvariantMassVsPtVsCentTofLikeRight;
   //HFT Mixed
   TH3F* mh3InvariantMassVsPtVsCentMixed;
   TH3F* mh3InvariantMassVsPtVsCentLikeMixed;
   TH3F* mh3InvariantMassVsPtVsCentTofMixed;
   TH3F* mh3InvariantMassVsPtVsCentTofLikeMixed;
   //HFT ratio QA
   TH2F* mh2Tpc1PtCent;
   TH2F* mh2Tpc1PhiVz;
   TH2F* mh2HFT1PtCent;
   TH2F* mh2HFT1PhiVz;
   TH2F* mh2Tpc1PtCentPartEtaVzPhi[anaCuts::nParticles][anaCuts::nEtasRatio][anaCuts::nVzsRatio][anaCuts::nPhisRatio];
   TH2F* mh2HFT1PtCentPartEtaVzPhi[anaCuts::nParticles][anaCuts::nEtasRatio][anaCuts::nVzsRatio][anaCuts::nPhisRatio];

   //HFT Dca
   TH3F* mh3DcaXyZPtCentPartEtaVzPhi[anaCuts::nParticles][anaCuts::nEtasDca][anaCuts::nVzsDca][anaCuts::nCentsDca];

   TH3F* mh3DcaPtCent;
   TH3F* mh3DcaXyPtCent;
   TH3F* mh3DcaZPtCent;

   TNtuple* mNtDstarBackgroungSameSign[anaCuts::nPtBins];
   TNtuple* mNtDstarBackgroungSideBand[anaCuts::nPtBins];
//Dstar
   TH2F* mh2InvariantMassVsPtDstar; 
   TH2F* mh2InvariantMassVsPtSBDstar; 
   TH2F* mh2InvariantMassVsPtLikeDstar;
   TH3F* mh3InvariantMassVsPtVsCentDstar;
   TH3F* mh3InvariantMassVsPtVsCentSBDstar;
   TH3F* mh3InvariantMassVsPtVsCentLikeDstar;
   TH2F* mh2InvariantMassVsPtDstarD0;
   TH2F* mh2InvariantMassVsPtSBD0;

};

#endif
