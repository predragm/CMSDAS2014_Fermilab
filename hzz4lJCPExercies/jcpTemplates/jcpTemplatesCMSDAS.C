//=================Usage=======================================================================================================================================//
//
// .L joetree.C
// .L jcpTemplatesCMSDAS.C
// jcpTemplatesCMSDAS(<hypothesis>, <process>, <processFile>, <finalState>)
//
// where
//   <hypothesis> : 0M, 0hP
//   <process>    : SM, ALT, qqZZ, Data
//   <processFile>: see Twiki
//   <finalState> : 2e2mu, 4mu, 4e, 4l
//
//=================Includes====================================================================================================================================//

// ROOT include
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TMath.h>
#include <TROOT.h>
#include "TRandom.h"

// C include
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdio.h>

// CMSDAS Fermilab include
#include "joetree.h"

using namespace std;


//=================Global definitions==========================================================================================================================//

// default cuts
const double CUT_ELPT = 7.;
const double CUT_MUPT = 5.;
const double CUT_ELETA = 2.5;
const double CUT_MUETA = 2.4;
const double CUT_MZ1LOW = 40.;
const double CUT_MZ1HIGH = 120.;
const double CUT_MZ2LOW = 12.;
const double CUT_MZ2HIGH = 120.;
const double CUT_M4LLOW = 121.;
const double CUT_M4LHIGH = 131.;
// directory names
const TString templatesDirJCP = "../jcpTemplates/templates1D";
// variable template binning
const int nbinsXps=21;
float binsXps[nbinsXps+1]={0.000, 0.030, 0.060, 0.100, 0.200, 0.300, 0.400, 0.500, 0.550, 0.600, 0.633, 0.666, 0.700, 0.733, 0.766, 0.800, 0.833, 0.866, 0.900, 0.933, 0.966, 1.000};
// print out settings
const int SILENT = true;
const int NEVT_PRINTOUT = 100000;
int printOutWidth = 12;
int printOutPrecision = 3;


//=================Implementations=============================================================================================================================//

//_______________________________________________________________________________________________________________________________________________
void analysisInit() {
    gErrorIgnoreLevel = kWarning;
    gErrorIgnoreLevel = kError;
}

//_______________________________________________________________________________________________________________________________________________
int normaliseHist(TH1D* &h1D, double norm = 1.){
    if (h1D->Integral()==0) return -1;
    h1D->Scale(norm/h1D->Integral());
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
void fillEmptyBinsHist1D(TH1* &h1D, double floor = .00001) {
    int nXbins=h1D->GetNbinsX();
    for(int i=1; i<=nXbins; i++){
        h1D->SetBinContent(i,h1D->GetBinContent(i)+floor);
    }
}

//_______________________________________________________________________________________________________________________________________________
void smoothAndNormaliseTemplate(TH1D* &h1D_KDm4l, double norm = 1.){
    cout.precision(2*printOutPrecision);
    int nXbins=h1D_KDm4l->GetNbinsX();
    // smooth
    TString smthAlg = "k5b";
    h1D_KDm4l->Smooth(1, smthAlg);
    // norm + floor + norm
    normaliseHist(h1D_KDm4l);
    double floor = .001/(nXbins);
    fillEmptyBinsHist1D(h1D_KDm4l,floor);
    normaliseHist(h1D_KDm4l, norm);
}

//_______________________________________________________________________________________________________________________________________________
void storeTemplateJCP(TH1D* h1D_KDm4l, TString histName, TString templateLocation, TString fOption){
    // store 1D template
    TFile ftJCP1D(templateLocation, fOption);
    ftJCP1D.cd();
    h1D_KDm4l->SetName(histName);
    h1D_KDm4l->Write();
    ftJCP1D.Close();
}

//_______________________________________________________________________________________________________________________________________________
// loop over the tree and fill the JCP template
int getHist1DJCP(TChain* tree, TString altModelTag, TH1* &h1D_kd, double ptElCut, double ptMuCut, double mZ2Cut, double mZ1Cut, double m4l_low, double m4l_high, bool scale = false){
    // counters
    int nEvtMassWindow = 0;

    // get the tree and number of tree entries
    h1D_kd->Sumw2();
    joetree jtree(tree);
    Long64_t nentries = jtree.fChain->GetEntries();

    // print-out
    if(!SILENT) {
        cout << setw(printOutWidth) << "nEntries: " << nentries << setw(printOutWidth) << endl;
        cout << "Event" << ":" << "Run" << ":" << "LumiSect" << ":" << "mass4l" << ":" << "massZ1" << ":" << "massZ2" << ":" << "weight" << ":" << "kd_ALT" << endl;
    }

    // loop over the tree and fill histogram
    for (Long64_t jentry=0; jentry<nentries;jentry++) { // tree loop
        Long64_t ientry = jtree.LoadTree(jentry);
        if (ientry < 0) break;
        jtree.fChain->GetEntry(jentry);

        // if events passes selection - fill hist with a proper weight
        double weight = (scale)? jtree.f_pu_weight*jtree.f_eff_weight: 1.;
        if (jtree.f_Z1mass<mZ1Cut || jtree.f_Z1mass>120 || jtree.f_Z2mass<mZ2Cut || jtree.f_Z2mass>120 || jtree.f_mass4l<m4l_low || jtree.f_mass4l>m4l_high) continue;
        // support only 0- and 0h+ models at the moment
        if (altModelTag=="0M")       {h1D_kd->Fill(jtree.f_KD_k0minus, weight);}
        else if (altModelTag=="0hP") {h1D_kd->Fill(jtree.f_KD_k0hplus, weight);}
        else return -1;

        // print-out
        if(!SILENT) {
            if (nEvtMassWindow<=NEVT_PRINTOUT) {
                cout << fixed;
                cout.precision(printOutPrecision);
                cout << jtree.f_event << ":" << jtree.f_run << ":" << jtree.f_lumi << ":" << jtree.f_mass4l << ":" << jtree.f_Z1mass << ":" << jtree.f_Z2mass << ":" << jtree.f_weight << ":" << jtree.f_KD << endl;
            }
        }
        // counter
        nEvtMassWindow++;
    }

    if(!SILENT) { // print-out
        cout  << setw(printOutWidth) << "nEventsInMassWindow: " << nEvtMassWindow << endl;
    }
    
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
// prepare JCP templates for the requested final states, hypotheses and processes
int getTemplateJCP(TString altModelTag, TString processNameTag, TString processFileName, TString sfinalState,
                   double elpT = CUT_ELPT, double mupT = CUT_MUPT, double mZ2_low = CUT_MZ2LOW, double mZ1_low = CUT_MZ1LOW, double m4l_low = CUT_M4LLOW, double m4l_high = CUT_M4LHIGH) {

    // prepare final state and event yields variables
    double eventYieldSig, eventYieldBkg, eventYieldDat; // expected event yields in +/- 5 GeV mass window, still TBC
    if      (sfinalState == "4mu")     {eventYieldSig = 7.1;  eventYieldBkg = 3.7; eventYieldDat = 9; }
    else if (sfinalState == "4e")      {eventYieldSig = 3.4;  eventYieldBkg = 1.7; eventYieldDat = 5; }
    else if (sfinalState == "2e2mu")   {eventYieldSig = 8.8;  eventYieldBkg = 4.6; eventYieldDat = 14; }
    else return -1;

    // get chains
    TChain* sigTree;
    sigTree = new TChain("HZZ4LeptonsAnalysisReduced");
    sigTree->Add(processFileName);

    // define histograms with variable binning
    TH1D *h1D_KDm4l = new TH1D(altModelTag, altModelTag + " - " + sfinalState, nbinsXps,binsXps);

    // loop over the tree and fill histograms
    getHist1DJCP(sigTree, altModelTag, h1D_KDm4l, elpT, mupT, mZ2_low, mZ1_low, m4l_low, m4l_high, true);

    // normalisation, with smoothing and filling of the empty bins
    if (processNameTag=="SM" ||
        processNameTag=="ALT")       {smoothAndNormaliseTemplate(h1D_KDm4l, eventYieldSig);}
    else if (processNameTag=="qqZZ") {smoothAndNormaliseTemplate(h1D_KDm4l, eventYieldBkg);}
    else if (processNameTag=="Data") {smoothAndNormaliseTemplate(h1D_KDm4l, eventYieldDat);} // test

    // store the 1D templates
    TString templateLocation = templatesDirJCP+"_D"+altModelTag+"/Hist_1D_"+processNameTag+"_8TeV_"+sfinalState+".root";
    storeTemplateJCP(h1D_KDm4l, "histRescaled", templateLocation, "recreate");

    return 0;
}

//_______________________________________________________________________________________________________________________________________________
// prepare JCP templates for the requested final state
void templatesJCP(TString altModelTag, TString processNameTag, TString processFileName, TString sfinalState = "4l") {
    cout << "[preparing 1D JCP templates for CMSDAS exercise, hypothesis: "+altModelTag+", process: "+processNameTag+", final state: "+sfinalState+"]" << endl;
    if (sfinalState == "4l") {
        getTemplateJCP(altModelTag, processNameTag, processFileName, "2e2mu");
        getTemplateJCP(altModelTag, processNameTag, processFileName, "4e");
        getTemplateJCP(altModelTag, processNameTag, processFileName, "4mu");
    } else {
        getTemplateJCP(altModelTag, processNameTag, processFileName, sfinalState);
    }
}

//_______________________________________________________________________________________________________________________________________________
// prepare JCP templates for given parameters
void jcpTemplatesCMSDAS(TString altModelTag = "0M", TString processNameTag = "SM", TString processFileName = "SMHiggsToZZTo4L_M-126.root", TString sfinalState = "4l"){
    analysisInit();
    templatesJCP(altModelTag, processNameTag, processFileName, sfinalState);
}
