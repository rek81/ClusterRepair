#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>
#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"
void fill_pulls(std::string iFile, TH1F *h_x, TH1F *h_y, TH1F *h_pullx = nullptr, TH1F *h_pully = nullptr, bool on_edge = false, bool bad_pix = false){
  TFile *f_data = TFile::Open(iFile.c_str());
  TTree *t1 = (TTree *)f_data->Get("pixelTree");
  //read event data
  const int SIZE = 20000;
  const int TKSIZE = 10000;
  Long64_t size  =  t1->GetEntries();
  Int_t TkN, TkNHits[TKSIZE];
  Float_t TkEta[TKSIZE], TkBeta[TKSIZE][20];
  Float_t ClRhLx[SIZE], ClRhLy[SIZE],  ClRhLxE[SIZE], ClRhLyE[SIZE],  ClSimTrEta[SIZE][10], ClSimHitLx[SIZE][10], ClSimHitLy[SIZE][10];
  Int_t ClRhIsOnEdge[SIZE], ClN, ClSimHitN[SIZE], ClType[SIZE], ClRhHasBadPix[SIZE], TkClI[SIZE][20], TkClN[SIZE], ClLayer[SIZE];
  t1->SetBranchAddress("ClN", &ClN);
  t1->SetBranchAddress("ClLayer", &ClLayer);
  t1->SetBranchAddress("TkN", &TkN);
  t1->SetBranchAddress("TkNHits", &TkNHits);
  t1->SetBranchAddress("TkClN", &TkClN);
  t1->SetBranchAddress("TkClI", &TkClI);
  t1->SetBranchAddress("TkEta", &TkEta);
  t1->SetBranchAddress("TkBeta", &TkBeta);
  t1->SetBranchAddress("ClSimHitN", &ClSimHitN);
  t1->SetBranchAddress("ClType", &ClType);
  t1->SetBranchAddress("ClRhHasBadPixels", &ClRhHasBadPix);
  t1->SetBranchAddress("ClRhIsOnEdge", &ClRhIsOnEdge);
  t1->SetBranchAddress("ClRhLx", &ClRhLx);
  t1->SetBranchAddress("ClRhLy", &ClRhLy);
  t1->SetBranchAddress("ClRhLxE", &ClRhLxE);
  t1->SetBranchAddress("ClRhLyE", &ClRhLyE);
  t1->SetBranchAddress("ClSimTrEta", &ClSimTrEta);
  t1->SetBranchAddress("ClSimHitLx", &ClSimHitLx);
  t1->SetBranchAddress("ClSimHitLy", &ClSimHitLy);
  for(int i =0; i< size; i++){
    //    std::cout << "in i loop" << std::endl;
    t1->GetEntry(i);
    //    std::cout << "got event" << std::endl;
    if(ClN > SIZE) printf("WARNING: Number of clusters (%i) greater than array size (%i), memory issues likely!!! \n \n \n", ClN, SIZE);
    //    std::cout << "after if ClN statement" << std::endl;
    for(int tkI = 0; tkI < TkN; tkI++){
      //if(TKNHits[tkI] <5) continue;
      for(int iClus=0; iClus < TkClN[tkI]; iClus++){
	int j = TkClI[tkI][iClus];
	float beta = TkBeta[tkI][iClus];
	//	std::cout << "into ClN loop" << std::endl;
	if(ClType[j]==1){ // on a track
	  //std::cout << "tk eta " << fabs(TkEta[0]) << " sim eta " << fabs(ClSimTrEta[j][0]) << std::endl;
	  //	  std::cout << "got event on track" << std::endl;
	  float dx(9999), dy(9999), pullx(9999), pully(9999);
	  int iKx(0), iKy(0), q(0);
	  for(int k=0; k<ClSimHitN[j]; k++){
	    if(TkEta[tkI] > 0.5 && TkEta[tkI] < 1.5){
	    if(fabs(ClSimHitLx[j][k] - ClRhLx[j]) < fabs(dx)){ 
	      dx = ClRhLx[j] - ClSimHitLx[j][k];
	      pullx =  dx/ClRhLxE[j];
	      iKx = k;
	    }
	    if(fabs(ClSimHitLy[j][k] - ClRhLy[j]) < fabs(dy)){ 
	      dy =  ClRhLy[j]  - ClSimHitLy[j][k];
	      pully =  dy/ClRhLyE[j];
	      iKy = k;
	    }
	  }
	  bool fill = (dx<9999) && (dy<9999);
	  if(on_edge) fill = fill && ClRhIsOnEdge[j];
	  if(bad_pix) fill = fill && ClRhHasBadPix[j];
	  if(fill){
	    //mult by 10000 to convert to microns
	    float to_microns = 1e4;
	    h_x->Fill(to_microns * dx);
	    h_y->Fill(to_microns * dy);
	    if(h_pullx != nullptr) h_pullx->Fill(pullx);
	    if(h_pully != nullptr) h_pully->Fill(pully);
	    //	    std::cout << "filled hists" << std::endl;
	  }
	}
      }
    }
  }
  }

  printf("printing means and std devs for %s (N = %.0f) \n", iFile.c_str(), h_x->Integral());
  if(on_edge) printf("Edge Clusters: \n");
  if(bad_pix) printf("BadPix Clusters: \n");
  printf("X: mean %.3e std_dev %.3e \n", h_x->GetMean(1), h_x->GetStdDev(1));
  printf("Y: mean %.3e std_dev %.3e \n", h_y->GetMean(1), h_y->GetStdDev(1));
  if(h_pullx != nullptr) printf("Pull X: mean %.3e std_dev %.3e \n", h_pullx->GetMean(1), h_pullx->GetStdDev(1));
  if(h_pully != nullptr) printf("Pull Y: mean %.3e std_dev %.3e \n", h_pully->GetMean(1), h_pully->GetStdDev(1));
  //h_x->Scale(1/h_x->Integral());
  //h_y->Scale(1/h_y->Integral());
  //if(h_pullx != nullptr) h_pullx->Scale(1/h_pullx->Integral());
  //if(h_pully != nullptr) h_pully->Scale(1/h_pully->Integral());
  return;
}
void draw_pulls(){
  float range = 300.;
  int nBins = 100;
  gROOT->SetBatch(1);
  gStyle->SetOptFit(1);
  TH1F *h_x = new TH1F("h_x", "X position residual; #Deltax (#mum)" , nBins, -range, range);
  TH1F *h_y = new TH1F("h_y", "Y position residual; #Deltay (#mum)" , nBins, -range, range);
  TH1F *h_pullx = new TH1F("h_pullx", "X Pull; X Pull" , nBins, -3, 3);
  TH1F *h_pully = new TH1F("h_pully", "Y Pull; Y Pull" , nBins, -3, 3);
  h_x->SetMarkerColor(kBlack);
  h_x->SetMarkerStyle(20);
  h_x->SetLineColor(kBlack);
  h_y->SetMarkerColor(kBlack);
  h_y->SetMarkerStyle(20);
  h_y->SetLineColor(kBlack);
  h_pullx->SetMarkerColor(kBlack);
  h_pullx->SetMarkerStyle(20);
  h_pullx->SetLineColor(kBlack);
  h_pully->SetMarkerColor(kBlack);
  h_pully->SetMarkerStyle(20);
  h_pully->SetLineColor(kBlack);
  //  string f_name("root://cmsxrootd.fnal.gov///store/user/rkowalsk/ClusterRepairCrabJobs/TTbar_14TeV_TuneCP5_cfi_GEN_SIM_SiPixQualANDDynEff/crab_14TeV_TuneCP5_cfi_GEN_SIM_2024_SiPixQualANDDyn_DamClTagFalse_pixStep/200604_203732/0000/PixelTree_QualANDDyn_DCFalse_all.root");
  string f_name("PixelTree_GenericReco_August2_2020_100events_noSplitCluster_trial3.root");
  //All clusters
  bool use_edge = false;
  bool use_badpix = false;
  fill_pulls(f_name, h_x, h_y, h_pullx, h_pully, use_edge, use_badpix);
  TCanvas *c1 = new TCanvas("c1", "", 0, 0, 800, 800);
  h_x->Draw("pe");
  h_x->Fit("gaus");
  h_x->GetFunction("gaus")->SetLineColor(kBlue);
  c1->SaveAs("all_residx_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  TCanvas *c2 = new TCanvas("c2", "", 0,0 , 800, 800);
  h_y->Draw("pe");
  h_y->Fit("gaus");
  h_y->GetFunction("gaus")->SetLineColor(kBlue);
  c2->SaveAs("all_residy_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  TCanvas *c3 = new TCanvas("c3", "", 0, 0, 800, 800);
  h_pullx->Draw("pe");
  h_pullx->Fit("gaus");
  h_pullx->GetFunction("gaus")->SetLineColor(kBlue);
  c3->SaveAs("all_pullx_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  TCanvas *c4 = new TCanvas("c4", "", 0,0 , 800, 800);
  h_pully->Draw("pe");
  h_pully->Fit("gaus");
  h_pully->GetFunction("gaus")->SetLineColor(kBlue);
  c4->SaveAs("all_pully_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  h_x->Reset(); h_y->Reset(); h_pullx->Reset(); h_pully->Reset();
  //Edge clusters
  use_edge = true;
  fill_pulls(f_name, h_x, h_y, h_pullx, h_pully, use_edge, use_badpix);
  TCanvas *c1e = new TCanvas("c1", "", 0, 0, 800, 800);
  h_x->Draw("pe");
  h_x->Fit("gaus");
  h_x->GetFunction("gaus")->SetLineColor(kBlue);
  c1e->SaveAs("edge_residx_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  TCanvas *c2e = new TCanvas("c2", "", 0,0 , 800, 800);
  h_y->Draw("pe");
  h_y->Fit("gaus");
  h_y->GetFunction("gaus")->SetLineColor(kBlue);
  c2e->SaveAs("edge_residy_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  TCanvas *c3e = new TCanvas("c3", "", 0, 0, 800, 800);
  h_pullx->Draw("pe");
  h_pullx->Fit("gaus");
  h_pullx->GetFunction("gaus")->SetLineColor(kBlue);
  c3e->SaveAs("edge_pullx_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  TCanvas *c4e = new TCanvas("c4", "", 0,0 , 800, 800);
  h_pully->Draw("pe");
  h_pully->Fit("gaus");
  h_pully->GetFunction("gaus")->SetLineColor(kBlue);
  c4e->SaveAs("edge_pully_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  h_x->Reset(); h_y->Reset(); h_pullx->Reset(); h_pully->Reset();
  //BadPix clusters
  use_edge = false;
  use_badpix = true;
  fill_pulls(f_name, h_x, h_y, h_pullx, h_pully, use_edge, use_badpix);
  TCanvas *c1b = new TCanvas("c1", "", 0, 0, 800, 800);
  h_x->Draw("pe");
  h_x->Fit("gaus");
  h_x->GetFunction("gaus")->SetLineColor(kBlue);
  c1b->SaveAs("badpix_residx_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  TCanvas *c2b = new TCanvas("c2", "", 0,0 , 800, 800);
  h_y->Draw("pe");
  h_y->Fit("gaus");
  h_y->GetFunction("gaus")->SetLineColor(kBlue);
  c2b->SaveAs("badpix_residy_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  TCanvas *c3b = new TCanvas("c3", "", 0, 0, 800, 800);
  h_pullx->Draw("pe");
  h_pullx->Fit("gaus");
  h_pullx->GetFunction("gaus")->SetLineColor(kBlue);
  c3b->SaveAs("badpix_pullx_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  TCanvas *c4b = new TCanvas("c4", "", 0,0 , 800, 800);
  h_pully->Draw("pe");
  h_pully->Fit("gaus");
  h_pully->GetFunction("gaus")->SetLineColor(kBlue);
  c4b->SaveAs("badpix_pully_L1_etaBTp5AND1p5_splitGR_August2_100events_noSplitCluster.png");
  return;
}

