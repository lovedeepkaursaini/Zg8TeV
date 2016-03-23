#define anaZgTree_cxx
#include "anaZgTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include <TMath.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include "cutBphoIDRun1.cc"
#include <fstream>
#include <string>
using namespace std;
int n_Vtx;
int numPU;
float puWt;
float puWtp5;
float puWtm5;
float pho_Pt;
float pho_Phi;
float pho_Eta;
float pho_SCEta;
int pho_HasPixelSeed ;
int pho_EleVeto ;
int pho_Loose;
int pho_Tight;
int pho_Medium;
float pho_MVA;
float fatJ_Pt;
float fatJ_rawPt;
float fatJ_Eta;
float fatJ_Phi;
float fatJ_En;
float fatJ_unc;
float fatJ_Mass;
float fatJ_prdMass;
float fatJ_Tau21;
float fatJ_Tau31;
int fatJ_Id;
int fatJ_partonId;
int fatJ_Ztag;
float fatJ_CSVBtag;
float fatJ_MUF;
float SDSJ0_Pt;
float SDSJ0_Eta;
float SDSJ0_Phi;
float SDSJ0_Mass;
float SDSJ0_En;
float SDSJ0_CSV;
float SDSJ1_Pt;
float SDSJ1_Eta;
float SDSJ1_Phi;
float SDSJ1_Mass;
float SDSJ1_En;
float SDSJ1_CSV;

float J1_Pt;
float J1_Eta;
float J1_Phi;
float J1_CSVBtag;
float J1_En;
float J2_Pt;
float J2_Eta;
float J2_Phi;
float J2_CSVBtag;
float J2_En;

void anaZgTree::clearVariables(){
  n_Vtx = -9;
  numPU = -9;
  puWt = -9;
  puWtp5 = -9;
  puWtm5 = -9;
  pho_Pt = -9;
  pho_Phi = -9;
  pho_Eta = -9;
  pho_SCEta = -9;
  pho_HasPixelSeed = -9;
  pho_EleVeto = -9;
  pho_Loose = 0;
  pho_Tight = 0;
  pho_Medium = 0;
  pho_MVA = -9;
  fatJ_Pt = -9;
  fatJ_rawPt = -9; 
  fatJ_Eta = -9;
  fatJ_Phi = -9;
  fatJ_unc = -9;
  fatJ_En = -9;
  fatJ_Mass = -9;
  fatJ_prdMass = -9;
  fatJ_Tau21 = -9;
  fatJ_Tau31 = -9;
  fatJ_Id = 0;
  fatJ_partonId = -99;
  fatJ_Ztag = 0;
  fatJ_CSVBtag = -9;
  fatJ_MUF = -9;
  SDSJ0_Pt = -9;
  SDSJ0_Eta = -9;
  SDSJ0_Phi = -9;
  SDSJ0_Mass = -9;
  SDSJ0_En = -9;
  SDSJ0_CSV = -9;
  SDSJ1_Pt = -9;
  SDSJ1_Eta = -9;
  SDSJ1_Phi = -9;
  SDSJ1_Mass = -9;
  SDSJ1_En = -9;
  SDSJ1_CSV = -9;

  J1_Pt = -9;
  J1_Eta = -9;
  J1_Phi = -9;
  J1_CSVBtag = -9;
  J1_En = -9;
  J2_Pt = -9;
  J2_Eta = -9;
  J2_Phi = -9;
  J2_CSVBtag = -9;
  J2_En = -9;
}


void anaZgTree::Loop(TString name)
{
  //
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  TFile* tmp = TFile::Open(name, "RECREATE");
  TH1::SetDefaultSumw2();
  TTree* miniTree = new TTree("miniTree", "miniTree");
  miniTree->Branch("event", &event);
  miniTree->Branch("run", &run);
  miniTree->Branch("puWt" ,&puWt ,"puWt/F");
  miniTree->Branch("puWtm5" ,&puWtm5 ,"puWtm5/F");
  miniTree->Branch("puWtp5" ,&puWtp5 ,"puWtp5/F");
  miniTree->Branch("numPU", &numPU, "numPU/I");
  miniTree->Branch("n_Vtx", &n_Vtx, "n_Vtx/I");
  miniTree->Branch("pho_Pt", &pho_Pt, "pho_Pt/F");
  miniTree->Branch("pho_Eta", &pho_Eta, "pho_Eta/F");
  miniTree->Branch("pho_Phi", &pho_Phi, "pho_Phi/F");
  miniTree->Branch("pho_SCEta", &pho_SCEta, "pho_SCEta/F");
  miniTree->Branch("pho_HasPixelSeed",&pho_HasPixelSeed,"pho_HasPixelSeed/I");
  miniTree->Branch("pho_EleVeto",&pho_EleVeto,"pho_EleVeto/I");
  miniTree->Branch("pho_Loose",&pho_Loose,"pho_Loose/I");
  miniTree->Branch("pho_Tight",&pho_Tight,"pho_Tight/I");
  miniTree->Branch("pho_Medium",&pho_Medium,"pho_Medium/I");
  miniTree->Branch("pho_MVA", &pho_MVA, "pho_MVA/F");
  miniTree->Branch("fatJ_Pt", &fatJ_Pt, "fatJ_Pt/F");
  miniTree->Branch("fatJ_rawPt", &fatJ_rawPt, "fatJ_rawPt/F");
  miniTree->Branch("fatJ_Eta", &fatJ_Eta, "fatJ_Eta/F");
  miniTree->Branch("fatJ_Phi", &fatJ_Phi, "fatJ_Phi/F");
  miniTree->Branch("fatJ_En", &fatJ_En, "fatJ_En/F");
  miniTree->Branch("fatJ_unc", &fatJ_unc, "fatJ_unc/F");
  miniTree->Branch("fatJ_Mass", &fatJ_Mass, "fatJ_Mass/F");
  miniTree->Branch("fatJ_prdMass", &fatJ_prdMass, "fatJ_prdMass/F");
  miniTree->Branch("fatJ_Tau21", &fatJ_Tau21, "fatJ_Tau21/F");
  miniTree->Branch("fatJ_Tau31", &fatJ_Tau31, "fatJ_Tau31/F");
  miniTree->Branch("fatJ_Id", &fatJ_Id, "fatJ_Id/I");
  miniTree->Branch("fatJ_partonId", &fatJ_partonId, "fatJ_partonId/I");
  miniTree->Branch("fatJ_Ztag", &fatJ_Ztag, "fatJ_Ztag/I");
  miniTree->Branch("fatJ_CSVBtag", &fatJ_CSVBtag, "fatJ_CSVBtag/F");
  miniTree->Branch("fatJ_MUF",&fatJ_MUF, "fatJ_MUF/F");
  miniTree->Branch("SDSJ0_Pt",&SDSJ0_Pt,"SDSJ0_Pt/F");
  miniTree->Branch("SDSJ0_Eta",&SDSJ0_Eta,"SDSJ0_Eta/F");
  miniTree->Branch("SDSJ0_Phi",&SDSJ0_Phi,"SDSJ0_Phi/F");
  miniTree->Branch("SDSJ0_Mass",&SDSJ0_Mass,"SDSJ0_Mass/F");
  miniTree->Branch("SDSJ0_En",&SDSJ0_En,"SDSJ0_En/F");

  miniTree->Branch("SDSJ0_CSV",&SDSJ0_CSV,"SDSJ0_CSV/F");
  miniTree->Branch("SDSJ1_Pt",&SDSJ1_Pt,"SDSJ1_Pt/F");
  miniTree->Branch("SDSJ1_Eta",&SDSJ1_Eta,"SDSJ1_Eta/F");
  miniTree->Branch("SDSJ1_Phi",&SDSJ1_Phi,"SDSJ1_Phi/F");
  miniTree->Branch("SDSJ1_Mass",&SDSJ1_Mass,"SDSJ1_Mass/F");
  miniTree->Branch("SDSJ1_En",&SDSJ1_En,"SDSJ1_En/F");

  miniTree->Branch("SDSJ1_CSV",&SDSJ1_CSV,"SDSJ1_CSV/F");
/*
  miniTree->Branch("J1_Pt", &J1_Pt, "J1_Pt/F");
  miniTree->Branch("J1_Eta", &J1_Eta, "J1_Eta/F");
  miniTree->Branch("J1_Phi", &J1_Phi, "J1_Phi/F");
  miniTree->Branch("J1_CSVBtag", &J1_CSVBtag, "J1_CSVBtag/F");
  miniTree->Branch("J1_En", &J1_En, "J1_En/F");
  miniTree->Branch("J2_Pt", &J2_Pt, "J2_Pt/F");
  miniTree->Branch("J2_Eta", &J2_Eta, "J2_Eta/F");
  miniTree->Branch("J2_Phi", &J2_Phi, "J2_Phi/F");
  miniTree->Branch("J2_CSVBtag", &J2_CSVBtag, "J2_CSVBtag/F");
  miniTree->Branch("J2_En", &J2_En, "J2_En/F");
*/

  TH1F* hCounter=new TH1F("hCounter","hCounter",20,0.,20);
  TH1D* hEvents_ = new TH1D("hEvents", "total processed", 2, 0, 2);
  TH1D* hist_nVtx=new TH1D("hist_nVtx","hist_nVtx",20,0,40);
  TH1D*  hist_nVtxPUwt=new TH1D("hist_nVtxPUwt","hist_nVtxPUwt",20,0,40);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    hEvents_->Fill(1);
    if(jentry % 1000 == 0) cout << "Processed " << jentry
                                << " events out of " <<nentries<< endl;
    clearVariables();
    event = event;
    run = run;
    n_Vtx = nVtx;
    if ( processCode != 0 ) {
      puWt = puWeight;
      puWtp5 = puWeight_p5;
      puWtm5 = puWeight_m5;
      numPU = 1;//nPU->at(1);
    }
    else{
      puWt = 1;
      puWtp5 = 1;
      puWtm5 = 1;
      numPU = 1;
    }
    hist_nVtx->Fill(nVtx);
    hist_nVtxPUwt->Fill(nVtx,puWt);
    
    hCounter->Fill(0);
    if(HLT[HLTIndex[45]]==0 || HLT[HLTIndex[45]]==-1) continue;//return;
    hCounter->Fill(1);
    if(nPho<1)continue;
    hCounter->Fill(2);
    vector <int> iphotons;
    for (int ipho = 0; ipho < nPho; ++ipho){
      // PRE-PHOTON SELECTION
      if((*phoEt)[ipho] < 150) continue;
      if( fabs((*phoSCEta)[ipho])>2.5) continue;
      if( fabs((*phoSCEta)[ipho])<1.566 && fabs((*phoSCEta)[ipho])>1.4442) continue;
      //if (!passPhotonID(ipho, 1)) continue;
      //if(phoEleVeto->at(ipho)==0)continue;
      //      if(phoIDMVA->at(ipho) < 0.374)continue;
      iphotons.push_back(ipho);
    }
    if(iphotons.size() < 1 ) continue;
    hCounter->Fill(3);
    
    TLorentzVector pho;
    int j = iphotons[0];
    pho.SetPtEtaPhiM((*phoEt)[j], (*phoEta)[j], (*phoPhi)[j], 0);
    
    pho_Pt = pho.Pt();
    pho_Eta = pho.Eta();
    pho_Phi = pho.Phi();
    pho_SCEta = phoSCEta->at(j);
    pho_HasPixelSeed = phohasPixelSeed->at(j);
    pho_EleVeto = phoEleVeto->at(j);
    pho_Loose = passPhotonID(j, 0);
    pho_Medium = passPhotonID(j, 1);
    pho_Tight = passPhotonID(j, 2);
    pho_MVA = 0;//phoIDMVA->at(j);
    //check cut flow
    if(pho.Pt()>160. && fabs(pho.Eta())<1.4442) hCounter->Fill(5);
    
    //Fatjet
    vector <int> indFJ;
    for (int ijet = 0; ijet < nfatJ; ++ijet){
      if(fatJPt->at(ijet)<150. )continue;
      if(fabs(fatJEta->at(ijet))>2.5) continue;
      if(fatJPFLooseId->at(ijet)==false)continue;
      double drjetgamma = dR((*fatJEta)[ijet],(*fatJPhi)[ijet],pho.Eta(),pho.Phi());
      if(drjetgamma<0.8)continue;
      indFJ.push_back(ijet);
    }
    
    if (indFJ.size()>0)  {
      
      if(pho.Pt()>160. && fabs(pho.Eta())<1.4442
	 && fatJPt->at(indFJ[0])>160 && 
	 fabs(fatJEta->at(indFJ[0]))<2.2 && fatJPFLooseId->at(indFJ[0])==true) hCounter->Fill(6);
      
      if(pho.Pt()>160. && fabs(pho.Eta())<1.4442
	 && fatJprunedM->at(indFJ[0])>30 && fatJPt->at(indFJ[0])>160 &&
	 fabs(fatJEta->at(indFJ[0]))<2.2 && fatJPFLooseId->at(indFJ[0])==true) hCounter->Fill(7);
      
      
      fatJ_unc = fatJjecUnc->at(indFJ[0]);
      fatJ_Mass = fatJMass->at(indFJ[0]);
      fatJ_prdMass = fatJprunedM->at(indFJ[0]);
      fatJ_Tau21 = fatJ_tau2->at(indFJ[0])/fatJ_tau1->at(indFJ[0]);
      fatJ_Tau31 = fatJ_tau3->at(indFJ[0])/fatJ_tau1->at(indFJ[0]);
      fatJ_Id = fatJPFLooseId->at(indFJ[0]);
      if ( processCode != 0 )  fatJ_partonId = fatJPartonID->at(indFJ[0]);
      if(fatJ_prdMass>75 && fatJ_prdMass<105 && fatJ_Tau21<0.5 )fatJ_Ztag = 1;
      fatJ_Pt = fatJPt->at(indFJ[0]);
      fatJ_rawPt = fatJRawPt->at(indFJ[0]);
      fatJ_Eta = fatJEta->at(indFJ[0]);
      fatJ_Phi = fatJPhi->at(indFJ[0]);
      fatJ_En =  fatJEn->at(indFJ[0]);
      fatJ_CSVBtag = fatJcsv->at(indFJ[0]);
      fatJ_MUF = fatJMUF->at(indFJ[0]);
      //cout<<"subjets: "<<nAK8softdropSubjet->at(indFJ[0])<<endl;
      if(nfatJSJ->at(indFJ[0])>1){
    //   cout<<"fJ mass: "<<fatJ_Mass<<", sujet pt: "<<(*AK8softdropSubjetPt)[indFJ[0]][0]<<'\t'<<(*AK8softdropSubjetPt)[indFJ[0]][1]<<endl;  
       //lead SJ
       SDSJ0_Pt = (*fatJSJPt)[indFJ[0]][0];
       SDSJ0_Eta = (*fatJSJEta)[indFJ[0]][0];
       SDSJ0_Phi = (*fatJSJPhi)[indFJ[0]][0];
       SDSJ0_Mass = (*fatJSJMass)[indFJ[0]][0];
       SDSJ0_En = (*fatJSJE)[indFJ[0]][0];
       SDSJ0_CSV = (*fatJSJcsv)[indFJ[0]][0];
       //sublead SJ
       SDSJ1_Pt = (*fatJSJPt)[indFJ[0]][1];
       SDSJ1_Eta = (*fatJSJEta)[indFJ[0]][1];
       SDSJ1_Phi = (*fatJSJPhi)[indFJ[0]][1];
       SDSJ1_Mass = (*fatJSJMass)[indFJ[0]][1];
       SDSJ1_En = (*fatJSJE)[indFJ[0]][1];
       SDSJ1_CSV = (*fatJSJcsv)[indFJ[0]][1];
      } 
    }    
    
    //resolved jets
   /* vector <int> ijets;
    for (int ijet = 0; ijet < nJet; ++ijet){
      if(jetPt->at(ijet)<30. )continue;
      if(fabs(jetEta->at(ijet))>2.5) continue;
      if(jetPFLooseId==false)continue;
      double drjetgamma = dR((*jetEta)[ijet],(*jetPhi)[ijet],pho.Eta(),pho.Phi());
      if(drjetgamma<0.5)continue;
      ijets.push_back(ijet);
    }
    if(ijets.size()>1){
    hCounter->Fill(7);
    
    TLorentzVector jet1;
    int q1 = ijets[0];
    jet1.SetPtEtaPhiE((*jetPt)[q1], (*jetEta)[q1], (*jetPhi)[q1], (*jetEn)[q1]);
    
    TLorentzVector jet2;
    int q2 = ijets[1];
    jet2.SetPtEtaPhiE((*jetPt)[q2], (*jetEta)[q2], (*jetPhi)[q2], (*jetEn)[q2]);
    
    J1_Pt = (*jetPt)[q1];
    J1_Eta = (*jetEta)[q1];
    J1_Phi = (*jetPhi)[q1];
    J1_CSVBtag = (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[q1];
    J1_En = (*jetEn)[q1];
    J2_Pt = (*jetPt)[q2];
    J2_Eta = (*jetEta)[q2];
    J2_Phi = (*jetPhi)[q2];
    J2_CSVBtag = (*jetpfCombinedInclusiveSecondaryVertexV2BJetTags)[q2];
    J2_En = (*jetEn)[q2];}*/
    miniTree->Fill();
    
  }	
  miniTree->Write();
  hEvents_->Write();
 hCounter->Write(); 
  tmp->Write();
  tmp->Close();
}


int main(int argc, char* argv[]) {
  std::string fName = argv[1];
  
  anaZgTree t(fName);
  fName = fName.erase(0,5);
  cout<<fName<<endl;
  t.Loop("MiniTree/minitree_"+fName+".root");
  return 0;
}


