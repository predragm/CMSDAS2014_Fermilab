#include <Riostream.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "Math/DistFunc.h"

#include "tdrstyle.C"

double Lumi_7TeV = 5.1;
double Lumi_8TeV = 19.6;

void extractJCPHypothesisTestResults(TString altHypothesis = "0M", bool quite = false){ 

  gErrorIgnoreLevel = kWarning;
  gErrorIgnoreLevel = kError;

  TString altModelName = "ALT";
  if (altHypothesis=="0M") {altModelName = "0^{-}";}
  else if (altHypothesis=="0hP") {altModelName = "0_{h}^{+}";}

  double Hist_Range_min = -30;
  double Hist_Range_max = +30;

  setTDRStyle( true );

  char fileName[128];
  sprintf(fileName,"qmu.root");
  TFile *fq=new TFile(fileName,"READ");
  TTree *t=(TTree*)fq->Get("q");

  float q,m,w;
  int type;
  t->SetBranchAddress("q",&q);
  t->SetBranchAddress("mh",&m);
  t->SetBranchAddress("weight",&w);
  t->SetBranchAddress("type",&type);

  TH1F *hSM=new TH1F("hSM","",4000,Hist_Range_min,Hist_Range_max);
  TH1F *hPS=new TH1F("hPS","",4000,Hist_Range_min,Hist_Range_max);
  TH1F *hObs=new TH1F("hObserved","",1000,Hist_Range_min,Hist_Range_max);
  if (!quite) cout<<"Start to lopp on tree in file "<<fileName<<endl;

  std::vector<float> v_SM, v_PS,v_Obs;

  for(int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if(i==0)
          if (!quite) cout<<"MASS in the TREE = "<<m<<endl<<endl;

    if(type<0){ //STD hypothesis
      hSM->Fill( -2*q );
      v_SM.push_back( -2*q );
    }
    else if(type>0){//ALT hypothesis (-> PS)
      hPS->Fill( -2*q );
      v_PS.push_back( -2*q );
    }
    else{
      hObs->Fill( -2*q );
      v_Obs.push_back( -2*q );
    }
  }//end loop on tree entries
  if (!quite) cout << "Actual data input: " << v_Obs[0] << endl;

  if (!quite) cout<<"Finished to loop, sorting vectors "<<v_SM.size()<<" "<<v_PS.size()<<" "<<v_Obs.size()<<endl;
  sort(v_SM.begin(),v_SM.end());//sort in ascending order
  sort(v_PS.begin(),v_PS.end()); 
  sort(v_Obs.begin(),v_Obs.end());
  int ntoysSM= hSM->GetEntries();
  int ntoysPS= hPS->GetEntries();

  //we assume that SM is on the right and PS on the left of zero
  if(v_PS.at(0)>v_SM.at(ntoysSM-1)){
    cout<<"Swapped distributions !!! The alternative model shouldstay on the negative side of the significance."<<endl;
    cout<<"Please edit the code and change the sign of q when filling histos and vectors in the loop on tree entries"<<endl;
    return 1;
  }
  float medianSM=v_SM.at(int(ntoysSM/2));
  float medianPS=v_PS.at(int(ntoysPS/2));
  if (!quite) {
    cout<<"Toys generated "<<ntoysSM<<"\t"<<ntoysPS<<endl;
    cout<<"Mean of SM/PS hypothesis: "<<hSM->GetMean()<<"   /   "<<hPS->GetMean()<<endl;
    cout<<"RMS  of SM/PS hypothesis: "<<hSM->GetRMS()<<"   /   "<<hPS->GetRMS()<<endl;
    cout<<"Median of SM/PS hypothesis: "<<medianSM<<"   /   "<<medianPS<<endl;
  }

  const float step=0.05;
  float coverage=0.0;
  float diff=10.0;
  float cut=v_PS.at(0)-step;
  float crosspoint=-99.0;
  int startSM=ntoysSM-1, startPS=0;

  // Possible scale location 1
//  hSM->Scale( 1/hSM->Integral() );
//  hPS->Scale( 1/hPS->Integral() );

  float integralSM=hSM->Integral();
  float integralPS=hPS->Integral();

//  cout<<integralSM<<endl; 
//  cout<<integralPS<<endl; 

  float tailSM=hSM->Integral(1,hSM->FindBin(medianPS))/integralSM;
  float tailPS=hPS->Integral(hPS->FindBin(medianSM),hPS->GetNbinsX())/integralPS;
//  cout<<tailSM<<endl;
//  cout<<tailPS<<endl;
  if (!quite) {
    if(tailSM>0)cout<<"Prob( q < median(P) | S ) = "<<tailSM<<"  ("<<ROOT::Math::normal_quantile_c(tailSM,1.0) <<" sigma)"<<endl;
    if(tailPS>0)cout<<"Prob( q > median(S) | P ) = "<<tailPS<<"  ("<<ROOT::Math::normal_quantile_c(tailPS,1.0) <<" sigma)"<<endl;
  }
  // Data point probability
  double Probability_Data_SM=hSM->Integral(0,hSM->FindBin(hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) ))-1)/integralSM;
  double Probability_Data_PS=hPS->Integral(0,hPS->FindBin(hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) ))-1)/integralPS;
  cout << "Data point at " << hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) ) << endl;

  if( Probability_Data_SM >= 1 || Probability_Data_SM <0 )
  {
    TF1 func( "f1","gaus",Hist_Range_min,Hist_Range_max );
    hSM->Fit( "f1", "QNO" );	//constant, mean, sigma
    Probability_Data_SM = ROOT::Math::normal_cdf( (hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) )-func.GetParameter(1))/func.GetParameter(2), 1, 0 );
  }
  if( Probability_Data_PS >= 1 || Probability_Data_PS <0 )
  {
    TF1 func( "f1","gaus",Hist_Range_min,Hist_Range_max );
    hPS->Fit( "f1", "QNO" );	//constant, mean, sigma
    Probability_Data_PS = ROOT::Math::normal_cdf( (hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) )-func.GetParameter(1))/func.GetParameter(2), 1, 0 );
  }

  if (!quite) {
    cout << "Prob( q >= qObs | q~SM Higgs+Bkg ) = " << static_cast<double>(1-Probability_Data_SM) << "  (";
    if( Probability_Data_SM<1 ) cout << ROOT::Math::normal_quantile_c( 1.0 - Probability_Data_SM, 1.0 ) << " sigma)\n";
    else cout << "--- sigma)\n";
    cout << "Prob( q >= qObs | q~Alt. sig. + Bkg ) = " << static_cast<double>(1-Probability_Data_PS) << "  (";
    if( Probability_Data_PS<1 ) cout << ROOT::Math::normal_quantile_c( 1.0 - Probability_Data_PS, 1.0 ) << " sigma)\n";
    else cout << "--- sigma)\n";
  }

  double  CLs = static_cast<double>((1-Probability_Data_PS)/(1-Probability_Data_SM));
  if (!quite) {
    cout << "CLs = " << CLs;
    cout << "  (";
    if( (1-Probability_Data_PS)/(1-Probability_Data_SM)>0 ) cout << ROOT::Math::normal_quantile_c( (1-Probability_Data_PS)/(1-Probability_Data_SM), 1.0 ) << " sigma)\n";
    else cout << "--- sigma)\n";
  }
  
  diff=10.0;
  coverage=0.0;
  for(int i=1;i<hSM->GetNbinsX();i++){
    
    float fracSM=hSM->Integral(1,i) / integralSM;
    float fracPS=hPS->Integral(i,hPS->GetNbinsX()) / integralPS;
    if(fabs(fracSM-fracPS)<diff){
      diff=fabs(fracSM-fracPS);
      coverage=(fracSM+fracPS)/2.0;
    }

  }

  float sepH= 2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
  if (!quite) cout<<"Separation from histograms = "<<sepH<<" with coverage "<<coverage<<endl;

  if (quite) {
      cout << "Separation = " << sepH << " sigma." << endl;
      cout << "CLs value  = " << CLs << endl;
  }

  // Possible scale location 2
  hSM->Scale( 1/hSM->Integral() );
  hPS->Scale( 1/hPS->Integral() );

  // plot
  Int_t ci;   // for color index setting
  gStyle->SetOptStat(0);
  TCanvas *c1=new TCanvas("c1","c1",800,800);

   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLeftMargin(0.15);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.15);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);

  c1->cd();
  hSM->Rebin(50);
  hPS->Rebin(50);

  if( (hSM->GetMaximum()) > (hPS->GetMaximum()) ) hSM->SetMaximum( (hSM->GetMaximum())*1.2  );
  else hSM->SetMaximum( (hPS->GetMaximum())*1.15 );

  hSM->SetXTitle("-2 #times ln(L_{"+altModelName+"}/L_{0+})");
  hSM->SetYTitle("Pseudoexperiments");
  hSM->SetLineWidth(2);

  ci = TColor::GetColor("#ff0000");
  hObs->SetLineColor( ci );
  hObs->SetFillColor( ci );
  hObs->SetLineWidth( 2 );

  // PRL style
  ci = TColor::GetColor("#ffcc33");
  hSM->SetFillColor( ci );

  ci = TColor::GetColor("#990000");
  hSM->SetLineColor(ci);
  hSM->SetLineStyle(2);
  hSM->SetLineWidth(2);
  hSM->SetMarkerStyle(20);

  ci = TColor::GetColor("#0099ff");
  hPS->SetFillColor(ci);
  hPS->SetFillStyle(3001);

  ci = TColor::GetColor("#0000ff");
  hPS->SetLineColor(ci);
  hPS->SetLineStyle(0);
  hPS->SetMarkerStyle(20);

  hSM->GetYaxis()->SetRangeUser(0,hSM->GetMaximum()*1.2);
  hSM->Draw();
  hSM->GetXaxis()->SetTitleOffset( 1.2 );
  hSM->GetXaxis()->SetTitleSize( 0.05 );
  hSM->GetYaxis()->SetTitleOffset( 1.3 );
  hSM->GetYaxis()->SetTitleSize( 0.05 );
  hPS->Draw("sames");

  
  hObs->Scale( 200 );
//  hObs->Draw("sames");

  // legend
  TLegend *leg = new TLegend(0.63,0.73,0.88,0.93,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hSM," 0^{+}","f");
  leg->AddEntry(hPS," "+altModelName+" ","f");
//  leg->AddEntry(hObs, "CMS data","l");
//  leg->AddEntry((TObject*)0, TString::Format("CLs value = %.2f", CLs),"");
//  leg->AddEntry((TObject*)0, TString::Format("Expected separation power: %.2f#sigma", fabs(sepH)),"");
  leg->SetTextFont(42);
  leg->Draw();


  TArrow Data_arrow( hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) ),0.2*hSM->GetMaximum(),hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) ),0.0, 0.05,"|->" );
  ci = TColor::GetColor("#ff0000");
//  Data_arrow.SetAngle(30);
  Data_arrow.SetLineWidth( 2 );
  Data_arrow.SetLineColor( ci );
  Data_arrow.SetFillColor( 1 );
  Data_arrow.SetFillStyle( 1001 );
//  Data_arrow.Draw();

  c1->SaveAs("plotSigSeparation_"+altHypothesis+".pdf");
  c1->SaveAs("plotSigSeparation_"+altHypothesis+".png");
  c1->SaveAs("plotSigSeparation_"+altHypothesis+".root");

}//end main

