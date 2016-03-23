const int    photonID_IsConv[2][3]                = { {0, 0, 0} ,             {0, 0, 0}             };
const double photonID_HoverE[2][3]                = { {0.05, 0.05, 0.05} ,    {0.05, 0.05, 0.05}    };
const double photonID_SigmaIEtaIEta[2][3]         = { {0.012, 0.011, 0.011} , {0.034, 0.033, 0.031} };
const double photonID_RhoCorrR03ChHadIso[2][3]    = { {2.6, 1.5, 0.7} ,       {2.3, 1.2, 0.5}       };
const double photonID_RhoCorrR03NeuHadIso_0[2][3] = { {3.5, 1.0, 0.4} ,       {2.9, 1.5, 1.5}       };
const double photonID_RhoCorrR03NeuHadIso_1[2][3] = { {0.04, 0.04, 0.04} ,    {0.04, 0.04, 0.04}    };
const double photonID_RhoCorrR03PhoIso_0[2][3]    = { {1.3, 0.7, 0.5} ,       {999, 1.0, 1.0}       };
const double photonID_RhoCorrR03PhoIso_1[2][3]    = { {0.005, 0.005, 0.005} , {0.005, 0.005, 0.005} };


double anaZgTree::dR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = acos(cos(phi1 - phi2));
  double deta = eta1 - eta2;
  return sqrt(dphi*dphi + deta*deta);
}

bool anaZgTree::fidEtaPass(double Eta){

  double fabsEta = TMath::Abs(Eta);
  if( fabsEta > 2.5) return false;
  if( 1.4442 < fabsEta && fabsEta < 1.566) return false;
  return true;
}

int anaZgTree::phoRegion(double absEta){
  int region = 0;
  if( absEta >= 1.0  ) region++;
  if( absEta >= 1.479) region++;
  if( absEta >= 2.0  ) region++;
  if( absEta >= 2.2  ) region++;
  if( absEta >= 2.3  ) region++;
  if( absEta >= 2.4  ) region++;
  return region;
}
//https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012
//// Selection implementation details for SPRING15 
double anaZgTree::phoEffArea03ChHad(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.012, 0.010, 0.014, 0.012, 0.016, 0.020, 0.012};
  return area[phoRegion(eta)];
}

double anaZgTree::phoEffArea03NeuHad(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.030, 0.057, 0.039, 0.015, 0.024, 0.039, 0.072};
  return area[phoRegion(eta)];
}

double anaZgTree::phoEffArea03Pho(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.148, 0.130, 0.112, 0.216, 0.262, 0.260, 0.266};
  return area[phoRegion(eta)];
}

bool anaZgTree::passPhotonID(int phoInd, int pho_ID_ind = 0) {
//cout<<phoInd<<endl;
  // phoInd - index of the photon in the tree
  // pho_ID_ind: 0 -- loose, 1 -- medium, 2 -- tight
  double eta = (*phoSCEta)[phoInd];
  double et = (*phoEt)[phoInd];
  
  // rho corrected isolations
  double Pho03ChHadIso =     (*phoPFChIso)[phoInd]   - rho2012 * phoEffArea03ChHad(eta);
  double Pho03NeuHadIso =    (*phoPFNeuIso)[phoInd]  - rho2012 * phoEffArea03NeuHad(eta);
  double Pho03PhoIso =       (*phoPFPhoIso)[phoInd]  - rho2012 * phoEffArea03Pho(eta);
  
  
  int region = 1; //barrel
  if(TMath::Abs( eta )< 1.4442) region = 0; //endcap
  bool phoPresel = fidEtaPass( eta ) &&
    et > 15 && //  Et cut
    ((*phoHoverE12)[phoInd] < photonID_HoverE[region][pho_ID_ind]) &&
    ((*phoSigmaIEtaIEta)[phoInd]<photonID_SigmaIEtaIEta[region][pho_ID_ind]) &&
    (Pho03NeuHadIso < (photonID_RhoCorrR03NeuHadIso_0[region][pho_ID_ind] + et * photonID_RhoCorrR03NeuHadIso_1[region][pho_ID_ind])) &&
    (Pho03PhoIso < (photonID_RhoCorrR03PhoIso_0[region][pho_ID_ind] + (*phoEt)[phoInd] * photonID_RhoCorrR03PhoIso_1[region][pho_ID_ind])) &&
    (Pho03ChHadIso < photonID_RhoCorrR03ChHadIso[region][pho_ID_ind])   ;
    //    (Pho03ChHadIso < photonID_RhoCorrR03ChHadIso[region][pho_ID_ind])   ;
  /*cout<<" phoSCEta: "<< eta
<<", phoEt: "<<et
<<"phoEleVeto: "<<(*phoEleVeto)[phoInd]  
<<", phoHoverE12: "<<(*phoHoverE12)[phoInd]
<<", Pho03NeuHadIso: "<<Pho03NeuHadIso
<<", Pho03PhoIso: "<<Pho03PhoIso
<<endl;
*/
  return phoPresel;
}
