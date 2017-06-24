#ifndef StPicoDstarAnaMaker_h
#define StPicoDstarAnaMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoDstarEvent
 *  simultaneously and do analysis.
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
#include "StThreeVectorF.hh"

class TString;
class TFile;
class TNtuple;
class StPicoEvent;
class StPicoD0Event;
class StKaonPion;
class StD0Pion;
class StPicoTrack;
class StPicoDstMaker;
class StPicoDstarAnaHists;
class StRefMultCorr;

class StPicoDstarAnaMaker : public StMaker
{
public:
   StPicoDstarAnaMaker(char const * name, TString const inputFilesList,
                    TString const outBaseName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil);
   virtual ~StPicoDstarAnaMaker();

   virtual Int_t Init();
   virtual Int_t Make();
   virtual Int_t Finish();

    int getEntries() const;
    void fillQaHistograms(bool b = true);
    void fillBackgroundTrees(bool b = true);

private:

   StPicoDstarAnaMaker() {}
   void readNextEvent();

    int  getD0PtIndex(StKaonPion const* ) const;
    bool isGoodTrigger(StPicoEvent const*) const;
    bool isGoodEvent(StPicoEvent const*, StThreeVectorF const& vtx) const;
    bool isGoodQaTrack(StPicoTrack const* ,StThreeVectorF const& momentum ,double dca) const;
    bool isGoodTrack(StPicoTrack const*, StThreeVectorF const&) const;
    bool isGoodSoftPionTrack(StPicoTrack const*, StThreeVectorF const&) const;
    bool isTpcPion(StPicoTrack const*) const;
    bool isTpcKaon(StPicoTrack const*) const;
    bool isTpcProton(StPicoTrack const*) const;
    bool isTofPion(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    bool isTofKaon(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    bool isTofProton(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    bool isGoodPair(StKaonPion const*) const;
    bool isSideBand(float m) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const& vtx) const;
    int trkHalf(StPicoTrack const*, StThreeVectorF const& vtx) const;
    bool isGoodD0(StKaonPion const*) const;
   StPicoDstMaker* mPicoDstMaker;
   StPicoD0Event* mPicoD0Event;
   StRefMultCorr* mGRefMultCorrUtil;

    TString mInputFilesList;
    TString mOutFileBaseName;
    TChain* mChain;
    int mEventCounter;
    bool mFillQaHists;
    bool mFillBackgroundTrees;

   // -------------- USER variables -------------------------
   // add your member variables here.
   // Remember that ntuples size can be really big, use histograms where appropriate
   StPicoDstarAnaHists* mHists;

   ClassDef(StPicoDstarAnaMaker, 1)
};

inline void StPicoDstarAnaMaker::fillQaHistograms(bool b) { mFillQaHists = b;}
inline void StPicoDstarAnaMaker::fillBackgroundTrees(bool b) { mFillBackgroundTrees = b;}

inline int StPicoDstarAnaMaker::getEntries() const
{
   return mChain ? mChain->GetEntries() : 0;
}

inline void StPicoDstarAnaMaker::readNextEvent()
{
   mChain->GetEntry(mEventCounter++);
}
#endif
