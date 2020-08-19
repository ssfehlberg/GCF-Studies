void analysis(){

  gStyle->SetPaintTextFormat("4.2f");gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetOptStat(0);

  //Load in the histogram file
  TFile *f=new TFile("hists_GCF_CCQE.root");

  const int num_cuts = 5;
  const char* cuts[num_cuts] = {"_b4_cuts","_pmis_cut","_muon_cut","_rec_cut","_lead_cut"};
  const char* titles_cuts[num_cuts] = {"Before Cuts","After P_{Miss} Cut","After Muon Mom. Cut","After Recoil Mom. Cut","After Leading Mom. Cut"};
  const int num_var = 4;
  const char* var[num_var] = {"_mom","_E","_theta","_phi"};
  const char* titles_var[num_var] = {"Momentum (GeV/c)","Energy (GeV)","cos(#theta)","#phi (Rad)"};
  const char* in_prog= "#scale[0.6]{MicroBooNE In-Progress}";
  //int proton_ylim[] = {24000,50000,18000,12000};
  //int muon_ylim[] = {17000,18000,14000,12000};
  double proton_ylim[] = {0.03,0.08,0.03,0.004};
  double muon_ylim[] = {0.02,0.03,0.03,0.004};
  TLatex *t = new TLatex(); //T Latex stuff                                                                                                                                                        
  t->SetNDC();
  t->SetTextAlign(22);
  
  TH1D* h_muon[num_cuts][num_var];
  TH1D* h_recoil[num_cuts][num_var];
  TH1D* h_leading[num_cuts][num_var];
  TCanvas* canv[num_cuts][num_var];
  TLegend* leg_proton[num_cuts][num_var];
  TLegend* leg_muon[num_cuts][num_var];
  TH1D* h_opening_angle_proton[num_cuts];
  TH1D* h_opening_angle_mu_leading[num_cuts];
  TH1D* h_delta_PT[num_cuts];
  TH1D* h_delta_alphaT[num_cuts];
  TH1D* h_delta_phiT[num_cuts];
  TCanvas* canv0[num_cuts];
  TCanvas* canv1[num_cuts];
  TCanvas* canv2[num_cuts];
  
  //Grab the histograms
  //////////////////////
  TH1D* h_Enu = (TH1D*)f->Get("h_Enu");
  TH1D* h_pmissT = (TH1D*)f->Get("h_pmissT");
  TH1D* h_pMiss = (TH1D*)f->Get("h_pMiss");
  TH1D* h_cos_gamma_lab = (TH1D*)f->Get("h_cos_gamma_lab");
  TH1D* h_cos_gamma_cm = (TH1D*)f->Get("h_cos_gamma_cm_raw");
  TH1D* h_cos_gamma_cm_with_cut = (TH1D*)f->Get("h_cos_gamma_cm_with_cut");
  TH2D* h_2D_p1vp2 = (TH2D*)f->Get("h_2D_p1vp2");
  TH2D* h_2D_cos_gamma_cm_vs_p_missT = (TH2D*)f->Get("h_2D_cos_gamma_cm_v_p_missT"); 
  
  for(int i = 0; i < num_cuts; i++){
    h_opening_angle_proton[i] = (TH1D*)f->Get(Form("h_opening_angle_protons%s",cuts[i]));
    h_opening_angle_mu_leading[i] = (TH1D*)f->Get(Form("h_opening_angle_mu_leading%s",cuts[i]));
    h_delta_PT[i] = (TH1D*)f->Get(Form("h_delta_PT%s",cuts[i]));
    h_delta_alphaT[i] = (TH1D*)f->Get(Form("h_delta_alphaT%s",cuts[i]));
    h_delta_phiT[i] = (TH1D*)f->Get(Form("h_delta_phiT%s",cuts[i]));
    
    for(int j = 0; j < num_var; j++){
      h_muon[i][j] = (TH1D*)f->Get(Form("h_muon%s%s",cuts[i],var[j]));
      h_recoil[i][j] = (TH1D*)f->Get(Form("h_recoil%s%s",cuts[i],var[j]));
      h_leading[i][j] = (TH1D*)f->Get(Form("h_leading%s%s",cuts[i],var[j]));    
    }
  }

  //Plotting Time!
  ////////////////
  TCanvas* fig0 = new TCanvas("fig0","fig0",3400,800);
  fig0->Divide(3,1);
  fig0->cd(1);
  h_pmissT->Draw("e1");
  h_pmissT->SetLineColor(kRed-7);
  h_pmissT->SetLineWidth(4);
  h_pmissT->SetTitle("Missing Transverse Momentum");
  h_pmissT->GetXaxis()->SetTitle("P_{miss}^{T} (GeV/c)");
  h_pmissT->GetYaxis()->SetTitle("Counts");
  //h_pmissT->SetMaximum();
  t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");
  fig0->cd(2);
  h_Enu->Draw("E1");
  h_Enu->SetLineColor(kRed-7);
  h_Enu->SetLineWidth(4);
  h_Enu->SetTitle("Energy of the Neutrino");
  h_Enu->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  h_Enu->GetYaxis()->SetTitle("Counts");
  //h_Enu->SetMaximum();
  t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");
  fig0->cd(3);
  h_pMiss->Draw("E1");
  h_pMiss->SetLineColor(kRed-7);
  h_pMiss->SetLineWidth(4);
  h_pMiss->SetTitle("Struck Neutron Momentum");
  h_pMiss->GetXaxis()->SetTitle("Momentum of Struck Neutron (GeV/c)");
  h_pMiss->GetYaxis()->SetTitle("Counts");
  //h_pMiss->SetMaximum();
  t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");
  fig0->cd();
  fig0->Print("images/raquel_24673.png");
  fig0->Print("images/raquel_24673.pdf");

  TCanvas* fig1 = new TCanvas("fig1","fig1",2000,1500);
  fig1->cd();
  h_2D_p1vp2->Draw("colz");
  h_2D_p1vp2->GetXaxis()->SetTitle("P_{Recoil} (GeV/c)");
  h_2D_p1vp2->GetYaxis()->SetTitle("P_{Leading} (GeV/c)");
  h_2D_p1vp2->SetTitle("");
  fig1->Print("images/raquel_Fig1.png");
  fig1->Print("images/raquel_Fig1.pdf");

  TCanvas* fig2 = new TCanvas("fig2","fig2",2000,1500);
  fig2->cd();
  h_cos_gamma_lab->Draw("E1");
  h_cos_gamma_lab->GetXaxis()->SetTitle("cos(#gamma_{Lab})");
  h_cos_gamma_lab->GetYaxis()->SetTitle("Counts");
  h_cos_gamma_lab->SetTitle("cos(#gamma_{Lab})");
  //h_cos_gamma_lab->SetMaximum();
  fig2->Print("images/raquel_Fig2.png");
  fig2->Print("images/raquel_Fig2.pdf");

  TCanvas* fig3 = new TCanvas("fig3","fig3",2000,1500);
  fig3->cd();
  h_2D_cos_gamma_cm_vs_p_missT->Draw("colz");
  h_2D_cos_gamma_cm_vs_p_missT->GetXaxis()->SetTitle("P_{miss}^{T} (GeV/c)");
  h_2D_cos_gamma_cm_vs_p_missT->GetYaxis()->SetTitle("cos(#gamma_{COM})");
  h_2D_cos_gamma_cm_vs_p_missT->SetTitle("");
  TLine* a = new TLine(0.3,-1,0.3,1);
  a->Draw("same");
  a->SetLineColor(kRed);
  a->SetLineWidth(4);
  fig3->Print("images/raquel_Fig3.png");
  fig3->Print("images/raquel_Fig3.pdf");

  TCanvas* sanity = new TCanvas("sanity","sanity",2000,1500);
  sanity->cd();
  h_cos_gamma_cm->Draw("e1");
  h_cos_gamma_cm->GetXaxis()->SetTitle("cos(#gamma_{COM})");
  h_cos_gamma_cm->GetYaxis()->SetTitle("Counts");
  h_cos_gamma_cm->SetTitle("cos(#gamma_{COM}): No Cut");
  h_cos_gamma_cm->SetMaximum(0.023);
  sanity->Print("images/raquel_cosgamma_cm.png");
  sanity->Print("images/raquel_cosgamma_cm.pdf");

  TCanvas* fig4 = new TCanvas("fig4","fig4",2000,1500);
  fig4->cd();
  h_cos_gamma_cm_with_cut->Draw("e1");
  h_cos_gamma_cm_with_cut->GetXaxis()->SetTitle("cos(#gamma_{COM})");
  h_cos_gamma_cm_with_cut->GetYaxis()->SetTitle("Counts");
  h_cos_gamma_cm_with_cut->SetMaximum(0.023);                                                                                                                                                                                                                                                    
  h_cos_gamma_cm_with_cut->SetTitle("cos(#gamma_{COM}) with Cut on P_{miss}^{T}");
  fig4->Print("images/raquel_Fig4.png");
  fig4->Print("images/raquel_Fig4.pdf");
  

  for(int i = 0; i < num_cuts; i++){

    canv0[i] = new TCanvas(Form("c0%s",cuts[i]),Form("c0%s",cuts[i]),2000,1500);
    canv0[i]->cd();
    h_opening_angle_proton[i]->Draw("hist");
    h_opening_angle_proton[i]->SetLineColor(kBlue-7);
    h_opening_angle_proton[i]->GetXaxis()->SetTitle("cos(Opening Angle btwn the Protons)");
    h_opening_angle_proton[i]->GetYaxis()->SetTitle("Counts"); //No. Events
    h_opening_angle_proton[i]->GetXaxis()->SetTitleSize(0.035);
    h_opening_angle_proton[i]->GetYaxis()->SetTitleSize(0.035);
    h_opening_angle_proton[i]->SetMaximum(0.015); //10000 //1700
    h_opening_angle_proton[i]->SetTitle(Form("%s: Opening Angle btwn the Protons",titles_cuts[i]));
    t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");
    canv0[i]->Print(Form("images/%s_opening_angle_protons.png",cuts[i]));
    canv0[i]->Print(Form("images/%s_opening_angle_protons.pdf",cuts[i]));

    canv1[i] = new TCanvas(Form("c1%s",cuts[i]),Form("c1%s",cuts[i]),2000,1500);
    canv1[i]->cd();
    h_opening_angle_mu_leading[i]->Draw("hist");
    h_opening_angle_mu_leading[i]->SetLineColor(kBlue-7);
    h_opening_angle_mu_leading[i]->GetXaxis()->SetTitle("cos(Opening Angle btwn the Muon and Leading Proton)");
    h_opening_angle_mu_leading[i]->GetYaxis()->SetTitle("Counts"); //No. Events
    h_opening_angle_mu_leading[i]->GetXaxis()->SetTitleSize(0.035);
    h_opening_angle_mu_leading[i]->GetYaxis()->SetTitleSize(0.035);
    h_opening_angle_mu_leading[i]->SetMaximum(0.015); //0.015//10000,//1900
    h_opening_angle_mu_leading[i]->SetTitle(Form("%s: cos(Opening Angle btwn the Muon and Leading Proton)",titles_cuts[i]));
    t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");
    canv1[i]->Print(Form("images/%s_opening_angle_mu_leading.png",cuts[i]));
    canv1[i]->Print(Form("images/%s_opening_angle_mu_leading.pdf",cuts[i]));

    canv2[i] = new TCanvas(Form("c2%s",cuts[i]),Form("c2%s",cuts[i]),3400,800);
    canv2[i]->Divide(3,1);
    canv2[i]->cd(1);
    h_delta_PT[i]->Draw("hist");
    h_delta_PT[i]->SetLineColor(kBlue-7);
    h_delta_PT[i]->GetXaxis()->SetTitle("#delta P_{T} (Gev/c)");
    h_delta_PT[i]->GetYaxis()->SetTitle("Counts"); //No. Events
    h_delta_PT[i]->GetXaxis()->SetTitleSize(0.04);
    h_delta_PT[i]->GetYaxis()->SetTitleSize(0.04);
    h_delta_PT[i]->SetMaximum(0.05);
    h_delta_PT[i]->SetTitle("");
    TLine* a1 = new TLine(0.259,0,0.259,0.05);
    a1->Draw("same");
    a1->SetLineColor(kBlack);
    a1->SetLineWidth(3);
    t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");

    canv2[i]->cd(2);
    h_delta_alphaT[i]->Draw("hist");
    h_delta_alphaT[i]->SetLineColor(kRed-7);
    h_delta_alphaT[i]->GetXaxis()->SetTitle("cos(#delta #alpha_{T})");
    h_delta_alphaT[i]->GetYaxis()->SetTitle("Counts"); //No. Events
    h_delta_alphaT[i]->GetXaxis()->SetTitleSize(0.04);
    h_delta_alphaT[i]->GetYaxis()->SetTitleSize(0.04);
    h_delta_alphaT[i]->SetMaximum(0.04);
    h_delta_alphaT[i]->SetTitle("");
    t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");

    canv2[i]->cd(3);
    h_delta_phiT[i]->Draw("hist");
    h_delta_phiT[i]->SetLineColor(kGreen+3);
    h_delta_phiT[i]->GetXaxis()->SetTitle("#delta #phi_{T} (Radians)");
    h_delta_phiT[i]->GetYaxis()->SetTitle("Counts"); //No. Events
    h_delta_phiT[i]->GetXaxis()->SetTitleSize(0.04); //0.04
    h_delta_phiT[i]->GetYaxis()->SetTitleSize(0.04);
    h_delta_phiT[i]->SetMaximum(0.04);
    h_delta_phiT[i]->SetTitle("");
    t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");

    canv2[i]->cd();
    t->DrawLatex(0.5,0.96,Form("Single Transverse Variables: %s",titles_cuts[i]));
    canv2[i]->Print(Form("images/stv%s.png",cuts[i]));
    canv2[i]->Print(Form("images/stv%s.pdf",cuts[i]));
    
    for(int j = 0; j < num_var; j++){
      
      canv[i][j] = new TCanvas(Form("c%s%s",cuts[i],var[j]),Form("c%s%s",cuts[i],var[j]),2800,800);
      canv[i][j]->Divide(2,1);
      canv[i][j]->cd(1);
      h_recoil[i][j]->Draw("hist");
      h_recoil[i][j]->SetLineColor(kBlue-7);
      h_recoil[i][j]->GetXaxis()->SetTitle(Form("%s",titles_var[j]));
      h_recoil[i][j]->GetYaxis()->SetTitle("Counts"); //No. Events
      h_recoil[i][j]->GetXaxis()->SetTitleSize(0.04);
      h_recoil[i][j]->GetYaxis()->SetTitleSize(0.04);
      h_recoil[i][j]->SetMaximum(proton_ylim[j]);
      h_recoil[i][j]->SetTitle("");
      h_leading[i][j]->Draw("histsame");
      h_leading[i][j]->SetLineColor(kRed-7);
      leg_proton[i][j] = new TLegend(0.6,0.7,0.8,0.83);
      leg_proton[i][j]->AddEntry(h_recoil[i][j],"Proton 1 (Recoil)","L");
      leg_proton[i][j]->AddEntry(h_leading[i][j],"Proton 2 (Leading)","L");
      leg_proton[i][j]->SetBorderSize(0);
      leg_proton[i][j]->SetTextSize(0.03);
      leg_proton[i][j]->SetFillColor(0);
      leg_proton[i][j]->Draw("same");
      t->DrawLatex(0.55,0.95,Form("#scale[0.8]{%s: Leading and Recoil Protons}",titles_var[j]));
      t->DrawLatex(0.8,0.88,"#scale[0.6]{MicroBooNE In-Progress}");

      canv[i][j]->cd(2);
      h_muon[i][j]->Draw("hist");
      h_muon[i][j]->SetLineColor(kGreen+3);
      h_muon[i][j]->GetXaxis()->SetTitle(Form("%s",titles_var[j]));
      h_muon[i][j]->GetYaxis()->SetTitle("Counts"); //No. Events
      h_muon[i][j]->GetXaxis()->SetTitleSize(0.04);
      h_muon[i][j]->GetYaxis()->SetTitleSize(0.04);
      h_muon[i][j]->SetTitle("");
      h_muon[i][j]->SetMaximum(muon_ylim[j]);
      leg_muon[i][j] = new TLegend(0.6,0.7,0.8,0.83,"");
      leg_muon[i][j]->AddEntry(h_muon[i][j],"Muon","L");
      leg_muon[i][j]->SetBorderSize(0);
      leg_muon[i][j]->SetTextSize(0.03);
      leg_muon[i][j]->SetFillColor(0);
      leg_muon[i][j]->Draw("same");
      t->DrawLatex(0.55,0.95,Form("#scale[0.8]{%s: Muon}",titles_var[j]));
      t->DrawLatex(0.8,0.88,"#scale[0.6]{MicroBooNE In-Progress}");

      canv[i][j]->cd();
      t->DrawLatex(0.5,0.96,Form("#scale[0.95]{%s}",titles_cuts[i]));      
      canv[i][j]->Print(Form("images/%s%s.png",cuts[i],var[j]));
      canv[i][j]->Print(Form("images/%s%s.pdf",cuts[i],var[j]));

    }
  }

  std::cout<<"----PROGRAM HAS FINISHED-----"<<std::endl;
  
  ///h_opening_angle_proton = (TH1D*)f->Get("h_opening_angle_proton");
  //h_opening_angle_mu_leading = (TH1D*)f->Get("h_opening_angle_mu_leading");




}
