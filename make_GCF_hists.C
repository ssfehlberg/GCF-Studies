#define make_GCF_hists_cxx
#include "make_GCF_hists.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TVector3.h>
#include <TFile.h>
#include <iostream>
#include <math.h>

void make_GCF_hists::Loop()
{
  
  //Open output file
  ////////////////////////
  TFile *outputfile = new TFile("hists_GCF_CCQE.root", "RECREATE");

  //Define all the Histograms
  ///////////////////////////
  Define_Histograms();
  
  //Begin Loop over all the Events:
  /////////////////////////////////////
  if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      //Remember these are also defined in the .h. Any changes you make here must also go there
      ////////////////////////////////////////////////////////////////////////////////////////
      TVector3 vBeam(0.,0.,Eneutrino); // z-direction is defined along the neutrino direction
      TVector3 vMuon(pk[0],pk[1],pk[2]); //muon momentum
      TVector3 vLead(pLead[0],pLead[1],pLead[2]); //leading proton momentum
      TVector3 vRec(pRec[0],pRec[1],pRec[2]); //recoil proton momentum
      TVector3 vq = vBeam - vMuon; // Momentum transfer
      TVector3 vmiss = vLead - vq; // Missing momentum

      Define_Parameters(vBeam, vMuon, vLead, vRec, vq, vmiss);
            
      //Fill da random histograms
      h_pMiss->Fill(vmiss.Mag(), weight); // Filling the histogram of missing mass magnitude
      h_QSq->Fill(QSq, weight); //the q2
      h_Enu->Fill(Eneutrino,weight);
      
      //Fill Histograms before any cuts:
      Fill_Histograms(0,vBeam, vMuon, vLead, vRec, vq, vmiss);
      n_events[0]++;

      //Cut Requiring |pMiss| > 300 MeV/c
      //if(vmiss.Mag() < 0.3) continue; // Placing a cut requiring |pMiss| > 300 MeV/c
      Fill_Histograms(1,vBeam, vMuon, vLead, vRec, vq, vmiss);
      n_events[1]++;
      
      //Cut on the Muon Momentum
      if(vMuon.Mag() < 0.1) continue;
      Fill_Histograms(2,vBeam, vMuon, vLead, vRec, vq, vmiss);
      n_events[2]++;
      
      //Cut on the Recoil Proton Momentum
      if(vRec.Mag() < 0.26) continue; 
      Fill_Histograms(3,vBeam, vMuon, vLead, vRec, vq, vmiss);
      n_events[3]++;
      
      //Cut on the Leading Proton Momentum
      if(vLead.Mag() < 0.26) continue;
      Fill_Histograms(4,vBeam, vMuon, vLead, vRec, vq, vmiss);
      n_events[4]++;

      //Do all the stuf for Raquel's Analysis
      compute_raquel(vq, vmiss, vBeam, vMuon, vLead, vRec);
      //Remove events with P_missT greater than 300 MeV
      if(p_missT >= 0.3) continue;
      n_pmissT++;
      h_cos_gamma_cm_with_cut->Fill(cos_gamma_cm,weight); //Figure 4 in Raquel's Note
      
   } //end loop over the entries


   std::cout<<"Total Number of Events: "<<nentries<<std::endl;
   std::cout<<"Number of Events Remaining after pMiss Cut: "<<n_events[1]<<" ("<<float(100.*(float(n_events[1])/float(nentries)))<<"%)"<<" This is a difference of "<<nentries-n_events[1]<<" Events."<<std::endl;
   std::cout<<"Number of Events Remaining after Muon Cut: "<<n_events[2]<<" ("<<float(100.*(float(n_events[2])/float(nentries)))<<"%)"<<" This is a difference of "<<n_events[1]-n_events[2]<<" Events."<<std::endl;
   std::cout<<"Number of Events Remaining after Recoil Cut: "<<n_events[3]<<" ("<<float(100.*(float(n_events[3])/float(nentries)))<<"%)"<<" This is a difference of "<<n_events[2]-n_events[3]<<" Events."<<std::endl;
   std::cout<<"Number of Events Remaining after Leading Cut: "<<n_events[4]<<" ("<<float(100.*(float(n_events[4])/float(nentries)))<<"%)"<<" This is a difference of "<<n_events[3]-n_events[4]<<" Events."<<std::endl;
   std::cout<<"Number of Events from Raquel's bit of Code: "<<n_pmissT<<std::endl;
   
   outputfile->cd();
   for(int i = 0; i < h_list.size(); i++){
     h_list[i]->Write();
   }
   for(int i = 0; i < h_list_2D.size(); i++){
     h_list_2D[i]->Write();
   }
   outputfile->Close();

   std::cout<<"----PROGRAM COMPLETED----"<<std::endl;
   
}
