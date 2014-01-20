void likelihoodScan1D_fa3results(){

  // open the file 
  TFile *f=new TFile("higgsCombine1D.MultiDimFit.mH126.root");
  TTree *t=(TTree*)f->Get("limit");
  t->Draw("2*deltaNLL:x", "deltaNLL > 0","PL");

  // get the graph with likelihood scan points
  TGraph *gr0 = (TGraph*) gROOT->FindObject("Graph")->Clone();
  gr0->SetName("scan1D");

  // plot the likelihood scan
  TCanvas *c1=new TCanvas("can1","CANVAS-SCAN1D",800,800);
  c1->cd();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gr0->SetTitle("");
  gr0->GetXaxis()->SetTitle("f_{a3}");
  gr0->GetYaxis()->SetTitle("2 #times #Delta(NLL)");
  gr0->GetYaxis()->SetTitleOffset(1.);
  gr0->GetXaxis()->SetLabelSize(0.04);
  gr0->GetYaxis()->SetLabelSize(0.04);
  gr0->Draw("AL");

  // save
  c1->SaveAs("plotLikelihoodScan_fa3.root");
  c1->SaveAs("plotLikelihoodScan_fa3.pdf");
  c1->SaveAs("plotLikelihoodScan_fa3.png");

  // get the inverse graph and extract exclusion limits
  TCanvas *c2=new TCanvas("can2","CANVAS-SCAN1D",800,800); c2->cd();
  t->Draw("x:2*deltaNLL", "deltaNLL > 0","PL");
  TGraph *gr0_inverse = (TGraph*) gROOT->FindObject("Graph")->Clone();
  gr0_inverse->SetName("scan1D_inv");

  cout << "fa3 upper limit @ 1 sigma = " << gr0_inverse->Eval(1.0) << endl;
  cout << "fa3 upper limit @ 95%C.L. = " << gr0_inverse->Eval(3.84) << endl;
}
