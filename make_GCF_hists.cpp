#include <iostream>
#include <math.h>
#include "TFile.h"
#include "TH1D.h"
#include "TVector3.h"
#include "TTree.h"

double m_mu = 0.10566; //mass muon
double m_proton = 0.938; //mass proton
int n_mag = 0;
int n_muon = 0;
int n_recoil = 0;
int n_leading = 0;

void make_GCF_hists()
{

  //Define input TFile and TTree
  /////////////////////////////
  TFile *inFile = new TFile("GCF_CCQE.root");
  TTree *inTree = (TTree *)inFile->Get("genT");

  //Open output file
  ////////////////////////
  TFile *outputfile = new TFile("hists_GCF_CCQE_with_weight.root", "RECREATE");

  // Define tree variables
  ////////////////////////
  double pk[3], pLead[3], pRec[3], pAm2[3];
  double Eneutrino, weight;
  int rec_type;
  inTree->SetBranchAddress("Eneutrino",&Eneutrino); // Generated neutrino energy
  inTree->SetBranchAddress("rec_type",&rec_type); // PDG particle code for recoil nucleon (proton=2212 or neutron=2112) 
  inTree->SetBranchAddress("pk",pk); // Scattered muon momentum
  inTree->SetBranchAddress("pLead",pLead); // Scattered leading proton momentum
  inTree->SetBranchAddress("pRec",pRec); // Recoil nucleon momentum
  inTree->SetBranchAddress("pAm2",pAm2); // Momentum of the residual nuclear A-2 system (you probably don't need this)
  inTree->SetBranchAddress("weight",&weight); // Event weight

  //Make histograms
  /////////////////
  vector<TH1 *> h_list;
  TH1D *h_pMiss = new TH1D("h_pMiss", "h_pMiss; pMiss [GeV/c];Counts", 30, 0, 1);
  h_list.push_back(h_pMiss);
  TH1D *h_QSq = new TH1D("h_QSq", "h_QSq; QSq [GeV^2];Counts", 30, 0, 3);
  h_list.push_back(h_QSq);

  const int num_cuts = 5;
  const char* cuts[num_cuts] = {"_b4_cuts","_pmis_cut","_muon_cut","_rec_cut","_lead_cut"};
  const int num_var = 4;
  const char* var[num_var] = {"_mom","_E","_theta","_phi"};
  int num_bins[num_var] = {50,50,50,100};
  double xlim_low[num_var] = {0.0,0.0,-1.5,-3.15};
  double xlim_high[num_var] = {2.0,3.0,1.5,3.15};
  double xlim_high_muon[num_var]={2.0,3.0,1.5,3.15};
  const char* xlabel[num_var] ={"P [GeV/c]","E [GeV]","cos(#theta)","#phi [Rad]"};
  
  TH1D* h_muon[num_cuts][num_var];
  TH1D* h_recoil[num_cuts][num_var];
  TH1D* h_leading[num_cuts][num_var];
  TH1D* h_opening_angle_protons[num_cuts];
  TH1D* h_opening_angle_mu_leading[num_cuts];
  
  for(int i = 0; i < num_cuts; i++){

    h_opening_angle_protons[i] = new TH1D(Form("h_opening_angle_protons%s",cuts[i]),Form("h_opening_angle_protons%s; Opening Angle btwn Two Protons; Counts",cuts[i]),35,-1.5,1.5); //50, 0, 1.5
    h_opening_angle_mu_leading[i] = new TH1D(Form("h_opening_angle_mu_leading%s",cuts[i]),Form("h_opening_angle_mu_leading%s;Opening Angle btwn Muon and Leading Proton; Counts",cuts[i]),35,-1.5,1.5);
    h_list.push_back(h_opening_angle_protons[i]);
    h_list.push_back(h_opening_angle_mu_leading[i]);

    for(int j = 0; j < num_var; j++){

      h_muon[i][j] = new TH1D(Form("h_muon%s%s",cuts[i],var[j]),Form(" h_muon%s%s ;%s; Counts",cuts[i],var[j],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_muon[j]);
      h_recoil[i][j] = new TH1D(Form("h_recoil%s%s",cuts[i],var[j]),Form("h_recoil%s%s ;%s; Counts",cuts[i],var[j],xlabel[j]),num_bins[j],xlim_low[j],xlim_high[j]);
      h_leading[i][j] = new TH1D(Form("h_leading%s%s",cuts[i],var[j]),Form("h_leading%s%s ;%s; Counts",cuts[i],var[j],xlabel[j]),num_bins[j],xlim_low[j],xlim_high[j]);
      h_list.push_back(h_muon[i][j]);
      h_list.push_back(h_recoil[i][j]);
      h_list.push_back(h_leading[i][j]); 

    }
  }

  // Set each histogram to "Sumw2" so that weights are handled correctly
  for (int i = 0; i < h_list.size(); i++){
    h_list[i]->Sumw2();
  }

  int total = inTree->GetEntries();
  
  // Loop over tree
  for (int i = 0; i < inTree->GetEntries(); i++) {

      inTree->GetEvent(i);

      TVector3 vBeam(0.,0.,Eneutrino); // z-direction is defined along the neutrino direction
      TVector3 vMuon(pk[0],pk[1],pk[2]); //muon momentum
      TVector3 vLead(pLead[0],pLead[1],pLead[2]); //leading proton momentum
      TVector3 vRec(pRec[0],pRec[1],pRec[2]); //recoil proton momentum
      TVector3 vq = vBeam - vMuon; // Momentum transfer
      TVector3 vmiss = vLead - vq; // Missing momentum
      double EMuon = sqrt(m_mu*m_mu + vMuon.Mag2()); // Muon Energy
      double ELead = sqrt(m_proton*m_proton + vLead.Mag2()); // Leading Energy
      double ERec = sqrt(m_proton*m_proton + vRec.Mag2()); // Recoil Energy
      double omega = Eneutrino - EMuon; // Energy transfer
      double QSq = vq.Mag2() - omega*omega; // Q-Squared
      double muon_theta = vMuon[2]/vMuon.Mag(); //note this is cos(theta)
      double lead_theta = vLead[2]/vLead.Mag(); //note this is cos(theta)
      double rec_theta = vRec[2]/vRec.Mag(); //note this is cos(theta)
      double muon_phi = atan2(vMuon[0],vMuon[1]); //this is actual angle
      double lead_phi =	atan2(vLead[0],vLead[1]); //this is actual angle
      double rec_phi =	atan2(vRec[0],vRec[1]); //this i actual angle
      double open_angle = ((vLead[0]*vRec[0])+(vLead[1]*vRec[1])+(vLead[2]*vRec[2]))/(vLead.Mag()*vRec.Mag()); //note this is the cos(opening angle)   
      double open_angle_mu = ((vLead[0]*vMuon[0])+(vLead[1]*vMuon[1])+(vLead[2]*vMuon[2]))/(vLead.Mag()*vMuon.Mag()); //note this is the cos(opening angle)  
      
      h_pMiss->Fill(vmiss.Mag(), weight); // Filling the histogram of missing mass magnitude
      h_QSq->Fill(QSq, weight); //the q2
      
      //Fill the histograms before any cuts are made
      /////////////////////////////////////////////
      h_muon[0][0]->Fill(vMuon.Mag(),weight);
      h_muon[0][1]->Fill(EMuon,weight);
      h_muon[0][2]->Fill(muon_theta,weight);
      h_muon[0][3]->Fill(muon_phi,weight);
      h_recoil[0][0]->Fill(vRec.Mag(),weight);
      h_recoil[0][1]->Fill(ERec,weight);
      h_recoil[0][2]->Fill(rec_theta,weight);
      h_recoil[0][3]->Fill(rec_phi,weight);
      h_leading[0][0]->Fill(vLead.Mag(),weight);
      h_leading[0][1]->Fill(ELead,weight);
      h_leading[0][2]->Fill(lead_theta,weight);
      h_leading[0][3]->Fill(lead_phi,weight);
      h_opening_angle_protons[0]->Fill(open_angle,weight);
      h_opening_angle_mu_leading[0]->Fill(open_angle_mu,weight);

      //Now  we apply a  bunch of cuts and fill histograms accordingly
      //////////////////////////////////////////////////////////////////
      //if(vmiss.Mag() < 0.3) continue; // Placing a cut requiring |pMiss| > 300 MeV/c
      n_mag++;
      h_muon[1][0]->Fill(vMuon.Mag(),weight);
      h_muon[1][1]->Fill(EMuon,weight);
      h_muon[1][2]->Fill(muon_theta,weight);
      h_muon[1][3]->Fill(muon_phi,weight);
      h_recoil[1][0]->Fill(vRec.Mag(),weight);
      h_recoil[1][1]->Fill(ERec,weight);
      h_recoil[1][2]->Fill(rec_theta,weight);
      h_recoil[1][3]->Fill(rec_phi,weight);
      h_leading[1][0]->Fill(vLead.Mag(),weight);
      h_leading[1][1]->Fill(ELead,weight);
      h_leading[1][2]->Fill(lead_theta,weight);
      h_leading[1][3]->Fill(lead_phi,weight);
      h_opening_angle_protons[1]->Fill(open_angle,weight);
      h_opening_angle_mu_leading[1]->Fill(open_angle_mu,weight);

      if(vMuon.Mag() < 0.1) continue; //Placing a cut  on the muon momentum 
      n_muon++;
      h_muon[2][0]->Fill(vMuon.Mag(),weight);
      h_muon[2][1]->Fill(EMuon,weight);
      h_muon[2][2]->Fill(muon_theta,weight);
      h_muon[2][3]->Fill(muon_phi,weight);
      h_recoil[2][0]->Fill(vRec.Mag(),weight);
      h_recoil[2][1]->Fill(ERec,weight);
      h_recoil[2][2]->Fill(rec_theta,weight);
      h_recoil[2][3]->Fill(rec_phi,weight);
      h_leading[2][0]->Fill(vLead.Mag(),weight);
      h_leading[2][1]->Fill(ELead,weight);
      h_leading[2][2]->Fill(lead_theta,weight);
      h_leading[2][3]->Fill(lead_phi,weight);
      h_opening_angle_protons[2]->Fill(open_angle,weight);
      h_opening_angle_mu_leading[2]->Fill(open_angle_mu,weight);
      
      if(vRec.Mag() < 0.26) continue; //Placing a cut on the Recoil Momentum
      n_recoil++;
      h_muon[3][0]->Fill(vMuon.Mag(),weight);
      h_muon[3][1]->Fill(EMuon,weight);
      h_muon[3][2]->Fill(muon_theta,weight);
      h_muon[3][3]->Fill(muon_phi,weight);
      h_recoil[3][0]->Fill(vRec.Mag(),weight);
      h_recoil[3][1]->Fill(ERec,weight);
      h_recoil[3][2]->Fill(rec_theta,weight);
      h_recoil[3][3]->Fill(rec_phi,weight);
      h_leading[3][0]->Fill(vLead.Mag(),weight);
      h_leading[3][1]->Fill(ELead,weight);
      h_leading[3][2]->Fill(lead_theta,weight);
      h_leading[3][3]->Fill(lead_phi,weight);
      h_opening_angle_protons[3]->Fill(open_angle,weight);
      h_opening_angle_mu_leading[3]->Fill(open_angle_mu,weight);
      
      if(vLead.Mag() < 0.26) continue; //Placing a cut on the Leading Momentum
      n_leading++;
      h_muon[4][0]->Fill(vMuon.Mag(),weight);
      h_muon[4][1]->Fill(EMuon,weight);
      h_muon[4][2]->Fill(muon_theta,weight);
      h_muon[4][3]->Fill(muon_phi,weight);
      h_recoil[4][0]->Fill(vRec.Mag(),weight);
      h_recoil[4][1]->Fill(ERec,weight);
      h_recoil[4][2]->Fill(rec_theta,weight);
      h_recoil[4][3]->Fill(rec_phi,weight);
      h_leading[4][0]->Fill(vLead.Mag(),weight);
      h_leading[4][1]->Fill(ELead,weight);
      h_leading[4][2]->Fill(lead_theta,weight);
      h_leading[4][3]->Fill(lead_phi,weight);
      h_opening_angle_protons[4]->Fill(open_angle,weight);
      h_opening_angle_mu_leading[4]->Fill(open_angle_mu,weight);

  } //end of loop over the entries


  std::cout<<"Total number of Events: "<<total<<std::endl;
  std::cout<<"Number of Events Remaining after pMiss Cut: "<<n_mag<<" ("<<float(100.*(float(n_mag)/float(total)))<<"%)"<<" This is a difference of "<<total-n_mag<<" Events."<<std::endl;
  std::cout<<"Number of Events Remaining after Muon Cut: "<<n_muon<<" ("<<float(100.*(float(n_muon)/float(total)))<<"%)"<<" This is a difference of "<<n_mag-n_muon<<" Events."<<std::endl;
  std::cout<<"Number of Events Remaining after Recoil Cut: "<<n_recoil<<" ("<<float(100.*(float(n_recoil)/float(total)))<<"%)"<<" This is a difference of "<<n_muon-n_recoil<<" Events."<<std::endl;
  std::cout<<"Number of Events Remaining after Leading Cut: "<<n_leading<<" ("<<float(100.*(float(n_leading)/float(total)))<<"%)"<<" This is a difference of "<<n_recoil-n_leading<<" Events."<<std::endl;
  
  // Write histograms
  outputfile->cd();

for(int i = 0; i < h_list.size(); i++){
    h_list[i]->Write();
 }

  outputfile->Close();

}


  

