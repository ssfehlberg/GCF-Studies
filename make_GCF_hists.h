#ifndef make_GCF_hists_h
#define make_GCF_hists_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TVector3.h>
#include <math.h>

// Header file for the classes stored in the TTree if any.
class make_GCF_hists {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.
   // Declaration of leaf types
  Double_t        Eneutrino;
  Int_t           rec_type; //PDG of struck nucleon 
  Double_t        pk[3]; //momentum of muon
  Double_t        pLead[3]; //mmentum of leading proton
  Double_t        pRec[3]; //momentum of the recoil proton
  Double_t        pAm2[3]; //momentum of A-2 system
  Double_t        weight; //weight
  
   // List of branches
   TBranch        *b_Eneutrino;   //!
   TBranch        *b_rec_type;   //!
   TBranch        *b_pk;   //!
   TBranch        *b_pLead;   //!
   TBranch        *b_pRec;   //!
   TBranch        *b_pAm2;   //!
   TBranch        *b_weight;   //!

  make_GCF_hists(TTree *tree=0);
  virtual ~make_GCF_hists();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  virtual void     Define_Histograms();
  virtual void     Define_Parameters(TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vRec, TVector3 vq, TVector3 vmiss);
  virtual void     compute_stvs(TVector3 p3mu, TVector3 p3p);
  virtual void     compute_raquel(TVector3 vq, TVector3 vmiss, TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vRec);
  virtual void     Fill_Histograms(int i,TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vRec, TVector3 vq, TVector3 vmiss);
  
private:

  // Mass values from GENIE v3.0.6
  /////////////////////////////////
  const double TARGET_MASS = 37.215526; // 40Ar, GeV
  const double NEUTRON_MASS = 0.93956541; // GeV
  const double PROTON_MASS = 0.93827208; // GeV
  const double BINDING_ENERGY = 0.0295; // 40Ar, GeV
  const double MUON_MASS = 0.10565837; // GeV

  //Histograms 
  /////////////////
  static const int num_cuts = 5;
  const char* cuts[num_cuts] = {"_b4_cuts","_pmis_cut","_muon_cut","_rec_cut","_lead_cut"};
  static const int num_var = 4;
  const char* var[num_var] = {"_mom","_E","_theta","_phi"};
  int num_bins[num_var] = {50,50,50,100};
  double xlim_low[num_var] = {0.0,0.0,-1.5,-3.15};
  double xlim_high[num_var] = {2.0,3.0,1.5,3.15};
  double xlim_high_muon[num_var]={2.0,3.0,1.5,3.15};
  const char* xlabel[num_var] ={"P [GeV/c]","E [GeV]","cos(#theta)","#phi [Rad]"};

  vector<TH1*> h_list;
  vector<TH2D*> h_list_2D;
  TH1D* h_pMiss;
  TH1D* h_QSq;
  TH1D* h_Enu; 
  TH1D* h_muon[num_cuts][num_var];
  TH1D* h_recoil[num_cuts][num_var];
  TH1D* h_leading[num_cuts][num_var];
  TH1D* h_opening_angle_protons[num_cuts];
  TH1D* h_opening_angle_mu_leading[num_cuts];
  TH1D* h_delta_PT[num_cuts];
  TH1D* h_delta_alphaT[num_cuts];
  TH1D* h_delta_phiT[num_cuts];

  TH1D* h_pmissT;
  TH1D* h_cos_gamma_lab;
  TH1D* h_cos_gamma_cm_raw;
  TH1D* h_cos_gamma_cm_with_cut;  
  TH2D* h_2D_cos_gamma_cm_v_p_missT;
  TH2D* h_2D_p1vp2;

  //Counters
  //////////
  int n_events[num_cuts] = {0};
  int n_pmissT = 0;

  //Parameters:
  /////////////////
  double EMuon; // Muon Energy
  double ELead; // Leading Energy
  double ERec; // Recoil Energy
  double omega; // Energy transfer
  double QSq; // Q-Squared
  double muon_theta; //note this is cos(theta)
  double lead_theta; //note this is cos(theta)
  double rec_theta; //note this is cos(theta)
  double muon_phi; //this is actual angle
  double lead_phi; //this is actual angle
  double rec_phi; //this is actual angle
  double open_angle; //note this is the cos(opening angle)
  double open_angle_mu; //note this is the cos(opening angle)  
  double delta_pT; //stv delta_pT
  double delta_alphaT; //stv delta_alphaT
  double delta_phiT; //stv delta_phiT
  double delta_phiT_sam; //to make the plot look better
  double delta_pL; //stv delta pL
  double pn; //reminant momentum ....?
  double p_missT; //transverse missing momentum
  double cos_gamma_lab; //cos(opening angle) in lab
  double cos_gamma_cm; //cos(opening angle) in cm
  double En; //energy of struck nucleon
  double betacm; //boost parameter
  //TLorentzVector lead; //leading proton in TLorentzvector to make boost easier
  //TLorentzVector rec;//recoil proton in TLorentzvector to make boost easier  
 
};

#endif

#ifdef make_GCF_hists_cxx

void make_GCF_hists::Define_Histograms(){

  h_pMiss = new TH1D("h_pMiss", "h_pMiss; pMiss [GeV/c];Counts", 30, 0, 1.6);
  h_QSq = new TH1D("h_QSq", "h_QSq; QSq [GeV^2];Counts", 30, 0, 3);
  h_Enu = new TH1D("h_Enu", "h_Enu; E_{#nu} [GeV]; Counts",10,0,2.5);
  h_pmissT = new TH1D("h_pmissT","h_pmissT; P_{miss}^{T} (GeV/C)",50,0,2);
  h_cos_gamma_lab = new TH1D("h_cos_gamma_lab","h_cos_gamma_lab; cos(#gamma_{Lab}); Counts",20,-1,1);
  h_cos_gamma_cm_raw = new TH1D("h_cos_gamma_cm_raw","h_cos_gamma_cm_raw; cos(#gamma_{COM}); Counts",20,-1.0,1.0);
  h_cos_gamma_cm_with_cut = new TH1D("h_cos_gamma_cm_with_cut","h_cos_gamma_cm_with_cut; cos(#gamma_{COM}); Counts",20,-1.0,1.0);
  h_2D_p1vp2 = new TH2D("h_2D_p1vp2","h_2D_p1vp2; Momentum of P_{Recoil} (Gev/c); Momentum of P_{Leading} (Gev/c)",50,0,3,50,0,3);
  h_2D_cos_gamma_cm_v_p_missT = new TH2D("h_2D_cos_gamma_cm_v_p_missT","h_2D_cos_gamma_cm_v_p_missT; P^{T}_{miss} (GeV/c); cos(#gamma_{COM})",20,0,1.6,20,-1.0,1.0);

  h_list.push_back(h_pMiss);
  h_list.push_back(h_QSq);
  h_list.push_back(h_Enu);
  h_list.push_back(h_pmissT);
  h_list.push_back(h_cos_gamma_lab);
  h_list.push_back(h_cos_gamma_cm_raw);
  h_list.push_back(h_cos_gamma_cm_with_cut);
  h_list_2D.push_back(h_2D_p1vp2);
  h_list_2D.push_back(h_2D_cos_gamma_cm_v_p_missT);
  
  for(int i = 0; i < num_cuts; i++){
    h_opening_angle_protons[i] = new TH1D(Form("h_opening_angle_protons%s",cuts[i]),Form("h_opening_angle_protons%s; Opening Angle btwn Two Protons; Counts",cuts[i]),35,-1.5,1.5); //50, 0, 1.5
    h_opening_angle_mu_leading[i] = new TH1D(Form("h_opening_angle_mu_leading%s",cuts[i]),Form("h_opening_angle_mu_leading%s;Opening Angle btwn Muon and Leading Proton; Counts",cuts[i]),35,-1.5,1.5);
    h_delta_PT[i] = new TH1D(Form("h_delta_PT%s",cuts[i]),Form("h_deltaPT%s;#delta P_{T} [GeV/c];Counts",cuts[i]),10,0,1);
    h_delta_alphaT[i] = new TH1D(Form("h_delta_alphaT%s",cuts[i]),Form("h_delta_alphaT%s; #delta #alpha_{T} [Deg.];Counts",cuts[i]),10,-1,1); //0,180
    h_delta_phiT[i] = new TH1D(Form("h_delta_phiT%s",cuts[i]),Form("h_delta_phiT%s; #delta #phi_{T} [Deg.];Counts",cuts[i]),10,-3.15,3.15); //0,180
    h_list.push_back(h_opening_angle_protons[i]);
    h_list.push_back(h_opening_angle_mu_leading[i]);
    h_list.push_back(h_delta_PT[i]);
    h_list.push_back(h_delta_alphaT[i]);
    h_list.push_back(h_delta_phiT[i]);

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
  for(int i = 0; i < h_list_2D.size(); i++){
    h_list_2D[i]->Sumw2();
  }
  
}
// Helper function for computing STVs (either reco or true)
//Not working at the moment
void make_GCF_hists::compute_stvs(TVector3 p3mu, TVector3 p3p){
  delta_pT = (p3mu + p3p).Perp();
  delta_phiT = acos( (-p3mu.X()*p3p.X() - p3mu.Y()*p3p.Y())
    / (p3mu.XYvector().Mod() * p3p.XYvector().Mod()) );
  TVector2 delta_pT_vec = (p3mu + p3p).XYvector();
  delta_alphaT = std::acos( (-p3mu.X()*delta_pT_vec.X()
    - p3mu.Y()*delta_pT_vec.Y())
    / (p3mu.XYvector().Mod() * delta_pT_vec.Mod()) );
  double Emu = std::sqrt(std::pow(MUON_MASS, 2) + p3mu.Mag2());
  double Ep = std::sqrt(std::pow(PROTON_MASS, 2) + p3p.Mag2());
  double R = TARGET_MASS + p3mu.Z() + p3p.Z() - Emu - Ep;
  // Estimated mass of the final remnant nucleus (CCQE assumption)
  double mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY;
  delta_pL = 0.5*R - (std::pow(mf, 2) + std::pow(delta_pT, 2)) / (2.*R);
  pn = std::sqrt( std::pow(delta_pL, 2) + std::pow(delta_pT, 2) );
}

void make_GCF_hists::Define_Parameters(TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vRec, TVector3 vq, TVector3 vmiss){ 
  //General stuff
  EMuon = sqrt(MUON_MASS*MUON_MASS + vMuon.Mag2()); // Muon Energy
  ELead = sqrt(PROTON_MASS*PROTON_MASS + vLead.Mag2()); // Leading Energy
  ERec = sqrt(PROTON_MASS*PROTON_MASS + vRec.Mag2()); // Recoil Energy
  omega = Eneutrino - EMuon; // Energy transfer
  QSq = vq.Mag2() - omega*omega; // Q-Squared
  muon_theta = vMuon[2]/vMuon.Mag(); //note this is cos(theta)
  lead_theta = vLead[2]/vLead.Mag(); //note this is cos(theta)
  rec_theta = vRec[2]/vRec.Mag(); //note this is cos(theta)
  muon_phi = atan2(vMuon[0],vMuon[1]); //this is actual angle
  lead_phi =	atan2(vLead[0],vLead[1]); //this is actual angle
  rec_phi =	atan2(vRec[0],vRec[1]); //this i actual angle
  open_angle = ((vLead[0]*vRec[0])+(vLead[1]*vRec[1])+(vLead[2]*vRec[2]))/(vLead.Mag()*vRec.Mag()); //note this is the cos(opening angle)   
  open_angle_mu = ((vLead[0]*vMuon[0])+(vLead[1]*vMuon[1])+(vLead[2]*vMuon[2]))/(vLead.Mag()*vMuon.Mag()); //note this is the cos(opening angle)

  //Stuff for STVs
  delta_pT = (vMuon + vLead).Perp();
  delta_phiT = std::acos( (-vMuon.X()*vLead.X() - vMuon.Y()*vLead.Y()) / (vMuon.XYvector().Mod() * vLead.XYvector().Mod()) );
  double x = (-vMuon.X()*vLead.X() - vMuon.Y()*vLead.Y()) / (vMuon.XYvector().Mod() * vLead.XYvector().Mod());
  delta_phiT_sam = std::atan(std::sqrt(1-std::pow(x,2))/x);
  TVector2 delta_pT_vec = (vMuon + vLead).XYvector();
  delta_alphaT = std::acos( (-vMuon.X()*delta_pT_vec.X()- vMuon.Y()*delta_pT_vec.Y()) / (vMuon.XYvector().Mod() * delta_pT_vec.Mod()) );
  double Emu = std::sqrt(std::pow(MUON_MASS, 2) + vMuon.Mag2());
  double Ep = std::sqrt(std::pow(PROTON_MASS, 2) + vLead.Mag2());
  double R = TARGET_MASS + vMuon.Z() + vLead.Z() - Emu - Ep;
  double mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY; //Estimated mass of the final remnant nucleus (CCQE assumption))
  delta_pL = 0.5*R - (std::pow(mf, 2) + std::pow(delta_pT, 2)) / (2.*R);
  pn = std::sqrt( std::pow(delta_pL, 2) + std::pow(delta_pT, 2) );
  
  //Stuff for Raquel
  p_missT = std::sqrt(std::pow((vMuon[0]+vLead[0]+vRec[0]),2)+std::pow(((vMuon[1]+vLead[1]+vRec[1])),2)); //transverse missing momentum 
  cos_gamma_lab = cos(vLead.Angle(vRec)); //cos(opening angle in lab frame) but fancy!
  En = std::sqrt(std::pow(NEUTRON_MASS,2) + vmiss.Mag2()); //energy of struck nucleon
  TLorentzVector lead(vLead[0],vLead[1],vLead[2],ELead);//leading proton TLorentzVector
  TLorentzVector rec(vRec[0],vRec[1],vRec[2],ERec); //Recoil proton TLorentzVector
  TLorentzVector betacm(vmiss[0] + vRec[0] + vBeam[0],vmiss[1] + vRec[1] + vBeam[1], vmiss[2] + vRec[2]+ vBeam[2], En + ERec + Eneutrino); //beta for CM
  TVector3 boost = betacm.BoostVector(); //the boost vector
  lead.Boost(-boost); //boost leading proton
  rec.Boost(-boost); //boost recoil proton
  cos_gamma_cm = cos(lead.Angle(rec.Vect())); //uses Lorentz Vectors

  /*not needed right now
  //see which proton has the largest momentum
  if(vLead.Mag() > vRec.Mag()){
    //std::cout<<"Leading Proton Had the Largest Momentum"<<std::endl;
      compute_stvs(vMuon,vLead);
    } else {
    //std::cout<<"Leading Proton Had the Largest Momentum"<<std::endl;
      compute_stvs(vMuon,vRec);
    }
  */
}

void make_GCF_hists::compute_raquel(TVector3 vq, TVector3 vmiss, TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vRec){

  h_2D_p1vp2->Fill(vRec.Mag(),vLead.Mag(),weight); //figure 1 in Raquels note
  h_pmissT->Fill(p_missT,weight); //plot of the transverse missing momentum for my sake
  h_cos_gamma_lab->Fill(cos_gamma_lab,weight); //cos(gamma of lab frame) for sanity check
  h_cos_gamma_cm_raw->Fill(cos_gamma_cm,weight); //plots of the cos of the opening angle in the CM frame for my sakeo
  h_2D_cos_gamma_cm_v_p_missT->Fill(p_missT,cos_gamma_cm,weight); //Figure 3 in Raquel's Note  
  //Figure 5 in Raquel's Note: Dont' know how to get P perp pn-p2
  //Figure 7 in Raquel's Note //need cut on previous plot
}

void make_GCF_hists::Fill_Histograms(int i,TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vRec, TVector3 vq, TVector3 vmiss){
  h_muon[i][0]->Fill(vMuon.Mag(),weight);
  h_muon[i][1]->Fill(EMuon,weight);
  h_muon[i][2]->Fill(muon_theta,weight);
  h_muon[i][3]->Fill(muon_phi,weight);
  h_recoil[i][0]->Fill(vRec.Mag(),weight);
  h_recoil[i][1]->Fill(ERec,weight);
  h_recoil[i][2]->Fill(rec_theta,weight);
  h_recoil[i][3]->Fill(rec_phi,weight);
  h_leading[i][0]->Fill(vLead.Mag(),weight);
  h_leading[i][1]->Fill(ELead,weight);
  h_leading[i][2]->Fill(lead_theta,weight);
  h_leading[i][3]->Fill(lead_phi,weight);
  h_opening_angle_protons[i]->Fill(open_angle,weight);
  h_opening_angle_mu_leading[i]->Fill(open_angle_mu,weight);
  h_delta_PT[i]->Fill(delta_pT,weight);
  h_delta_alphaT[i]->Fill(cos(delta_alphaT),weight);
  h_delta_phiT[i]->Fill(delta_phiT_sam,weight);
}
make_GCF_hists::make_GCF_hists(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("GCF_CCQE.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("GCF_CCQE.root");
      }
      f->GetObject("genT",tree);

   }
   Init(tree);
}

make_GCF_hists::~make_GCF_hists()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t make_GCF_hists::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t make_GCF_hists::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void make_GCF_hists::Init(TTree *tree)
{
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Eneutrino", &Eneutrino, &b_Eneutrino);
   fChain->SetBranchAddress("rec_type", &rec_type, &b_rec_type);
   fChain->SetBranchAddress("pk", pk, &b_pk);
   fChain->SetBranchAddress("pLead", pLead, &b_pLead);
   fChain->SetBranchAddress("pRec", pRec, &b_pRec);
   fChain->SetBranchAddress("pAm2", pAm2, &b_pAm2);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   Notify();
}

Bool_t make_GCF_hists::Notify()
{

   return kTRUE;
}

void make_GCF_hists::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t make_GCF_hists::Cut(Long64_t entry)
{
   return 1;
}
#endif // #ifdef make_GCF_hists_cxx
