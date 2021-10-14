#include <iostream>
#include <vector>
#include <string>
#include <TH1F.h>
#include <fstream>
using namespace std;

void merge()
{
    vector<std::string> names_files {"ttbarh.root", "ttbarZ.root", "Zprime_tata_350.root", "Zprime_tata_1500.root"};

    vector<std::string> names_files_wo_ext {"$t\\bar{t}h$", "$t\\bar{t}Z$", "$m(Z')=350\\,GeV$", "$m(Z')=1500\\,GeV$"};

    vector<std::string> names {"t#bar{t}h", "t#bar{t}Z", "m(Z')= 350 GeV", "m(Z')= 1500 GeV"};

    vector<std::string> plots {"PT_b1", "ETA_b1", "PHI_b1", "PT_b2", "ETA_b2", "PHI_b2", "PT_tau1", "ETA_tau1", "PHI_tau1", "PT_tau2", "ETA_tau2", "PHI_tau2", "PT_j1", "ETA_j1", "PHI_j1", "PT_j2", "ETA_j2", "PHI_j2", "PT_j3", "ETA_j3", "PHI_j3", "PT_j4", "ETA_j4", "PHI_j4", "DeltaR_tau1tau2", "DeltaR_b1b2", "DeltaR_j1j2", "DeltaR_j1j3", "DeltaR_j1j4", "DeltaR_j2j3", "DeltaR_j2j4", "DeltaR_j3j4", "DeltaPhi_tau1tau2", "DeltaPhi_b1b2", "DeltaPhi_j1j2", "DeltaPhi_j1j3", "DeltaPhi_j1j4", "DeltaPhi_j2j3", "DeltaPhi_j2j4", "DeltaPhi_j3j4", "DeltaPt_tau1tau2", "DeltaPt_b1b2", "DeltaPt_j1j2", "DeltaPt_j1j3", "DeltaPt_j1j4", "DeltaPt_j2j3", "DeltaPt_j2j4", "DeltaPt_j3j4", "M_Zp", "M_j1j2", "M_j1j3", "M_j1j4", "M_j2j3", "M_j2j4", "M_j3j4", "M_j1j2b1", "M_j1j3b1", "M_j1j4b1", "M_j2j3b1", "M_j2j4b1", "M_j3j4b1", "M_j1j2b2", "M_j1j3b2", "M_j1j4b2", "M_j2j3b2", "M_j2j4b2", "M_j3j4b2", "StMet"};

    vector<std::string> plots_names {"p_{T}(b_{1})", "#eta(b_{1})", "#phi(b_{1})", "p_{T}(b_{2})", "#eta(b_{2})", "#phi(b_{2})", "p_{T}(#tau_{1})", "#eta(#tau_{1})", "#phi(#tau_{1})", "p_{T}(#tau_{2})", "#eta(#tau_{2})", "#phi(#tau_{2})", "p_{T}(j_{1})", "#eta(j_{1})", "#phi(j_{1})", "p_{T}(j_{2})", "#eta(j_{2})", "#phi(j_{2})", "p_{T}(j_{3})", "#eta(j_{3})", "#phi(j_{3})", "p_{T}(j_{4})", "#eta(j_{4})", "#phi(j_{4})", "#Delta R(#tau_{1},#tau_{2})", "#Delta R(b_{1},b_{2})", "#Delta R(j_{1},j_{2})", "#Delta R(j_{1},j_{3})", "#Delta R(j_{1},j_{4})", "#Delta R(j_{2},j_{3})", "#Delta R(j_{2},j_{4})", "#Delta R(j_{3},j_{4})", "#Delta #phi(#tau_{1},#tau_{2})", "#Delta #phi(b_{1},b_{2})", "#Delta #phi(j_{1},j_{2})", "#Delta #phi(j_{1},j_{3})", "#Delta #phi(j_{1},j_{4})", "#Delta #phi(j_{2},j_{3})", "#Delta #phi(j_{2},j_{4})", "#Delta #phi(j_{3},j_{4})", "#Delta p_{T}(#tau_{1},#tau_{2})", "#Delta p_{T}(b_{1},b_{2})", "#Delta p_{T}(j_{1},j_{2})", "#Delta p_{T}(j_{1},j_{3})", "#Delta p_{T}(j_{1},j_{4})", "#Delta p_{T}(j_{2},j_{3})", "#Delta p_{T}(j_{2},j_{4})", "#Delta p_{T}(j_{3},j_{4})", "M(Zp)", "M(j_{1},j_{2})", "M(j_{1},j_{3})", "M(j_{1},j_{4})", "M(j_{2},j_{3})", "M(j_{2},j_{4})", "M(j_{3},j_{4})", "M(j_{1},j_{2},b_{1})", "M(j_{1},j_{3},b_{1})", "M(j_{1},j_{4},b_{1})", "M(j_{2},j_{3},b_{1})", "M(j_{2},j_{4},b_{1})", "M(j_{3},j_{4},b_{1})", "M(j_{1},j_{2},b_{2})", "M(j_{1},j_{3},b_{2})", "M(j_{1},j_{4},b_{2})", "M(j_{2},j_{3},b_{2})", "M(j_{2},j_{4},b_{2})", "M(j_{3},j_{4},b_{2})", "StMet"}; 
    
    vector<std::string> x_labels {"p_{T}(b_{1}) [GeV]", "#eta(b_{1}) [a.u]", "#phi(b_{1}) [a.u]", "p_{T}(b_{2}) [GeV]", "#eta(b_{2}) [a.u]", "#phi(b_{2}) [a.u]", "p_{T}(#tau_{1}) [GeV]", "#eta(#tau_{1}) [a.u]", "#phi(#tau_{1}) [a.u]", "p_{T}(#tau_{2}) [GeV]", "#eta(#tau_{2}) [a.u]", "#phi(#tau_{2}) [a.u]", "p_{T}(j_{1}) [GeV]", "#eta(j_{1}) [a.u]", "#phi(j_{1}) [a.u]", "p_{T}(j_{2}) [GeV]", "#eta(j_{2}) [a.u]", "#phi(j_{2}) [a.u]", "p_{T}(j_{3}) [GeV]", "#eta(j_{3}) [a.u]", "#phi(j_{3}) [a.u]", "p_{T}(j_{4}) [GeV]", "#eta(j_{4}) [a.u]", "#phi(j_{4}) [a.u]", "#Delta R(#tau_{1},#tau_{2}) [a.u]", "#Delta R(b_{1},b_{2}) [a.u]", "#Delta R(j_{1},j_{2}) [a.u]", "#Delta R(j_{1},j_{3}) [a.u]", "#Delta R(j_{1},j_{4}) [a.u]", "#Delta R(j_{2},j_{3}) [a.u]", "#Delta R(j_{2},j_{4}) [a.u]", "#Delta R(j_{3},j_{4}) [a.u]", "#Delta #phi(#tau_{1},#tau_{2}) [a.u]", "#Delta #phi(b_{1},b_{2}) [a.u]", "#Delta #phi(j_{1},j_{2}) [a.u]", "#Delta #phi(j_{1},j_{3}) [a.u]", "#Delta #phi(j_{1},j_{4}) [a.u]", "#Delta #phi(j_{2},j_{3}) [a.u]", "#Delta #phi(j_{2},j_{4}) [a.u]", "#Delta #phi(j_{3},j_{4}) [a.u]", "#Delta p_{T}(#tau_{1},#tau_{2}) [GeV]", "#Delta p_{T}(b_{1},b_{2}) [GeV]", "#Delta p_{T}(j_{1},j_{2}) [GeV]", "#Delta p_{T}(j_{1},j_{3}) [GeV]", "#Delta p_{T}(j_{1},j_{4}) [GeV]", "#Delta p_{T}(j_{2},j_{3}) [GeV]", "#Delta p_{T}(j_{2},j_{4}) [GeV]", "#Delta p_{T}(j_{3},j_{4}) [GeV]", "M(Zp) [GeV]", "M(j_{1},j_{2}) [GeV]", "M(j_{1},j_{3}) [GeV]", "M(j_{1},j_{4}) [GeV]", "M(j_{2},j_{3}) [GeV]", "M(j_{2},j_{4}) [GeV]", "M(j_{3},j_{4}) [GeV]", "M(j_{1},j_{2},b_{1}) [GeV]", "M(j_{1},j_{3},b_{1}) [GeV]", "M(j_{1},j_{4},b_{1}) [GeV]", "M(j_{2},j_{3},b_{1}) [GeV]", "M(j_{2},j_{4},b_{1}) [GeV]", "M(j_{3},j_{4},b_{1}) [GeV]", "M(j_{1},j_{2},b_{2}) [GeV]", "M(j_{1},j_{3},b_{2}) [GeV]", "M(j_{1},j_{4},b_{2}) [GeV]", "M(j_{2},j_{3},b_{2}) [GeV]", "M(j_{2},j_{4},b_{2}) [GeV]", "M(j_{3},j_{4},b_{2}) [GeV]", "StMet [GeV]"}; 
    
    
    vector<int> colors {3, 6, 5, 2};

    vector<int> linestyles {1, 1, 10, 10};

    TList *l = new TList();

    for(int i=0; i<plots.size(); i++)
    {
        THStack *hs = new THStack("hs", plots_names[i].c_str());
        TCanvas *c2 = new TCanvas(plots[i].c_str(),"Histos",1280,1024);  
        Double_t x_1,x_2;
        //if (plots[i]=="M_dijet_partially" || plots[i]=="M_b_dijet_partially" || plots[i]=="M_b_dijet_fully" || plots[i]=="M_dijet_no" || plots[i]=="M_b_not_used_diff" || plots[i]=="M_b1b2"){
            //x_1 = 0.15;
            //x_2 = 0.35;
        //}

        //else if (plots[i]=="M_b_not_used_diff_before" || plots[i]=="M_b1b2_before"){
            //x_1 = 0.45;
            //x_2 = 0.65;
        //}
        //else{
            x_1 = 0.65;
            x_2 = 0.85;
        //}

        auto legend = new TLegend(x_1,0.65,x_2,0.85);
    for (int j=0; j<names.size(); j++)
    {    
        TFile f(names_files[j].c_str());
        TH1F *h = (TH1F*)f.Get(plots[i].c_str());
        h->SetDirectory(0);

        if ((plots[i] == "N_Merged") || (plots[i] == "Eff_mu") || (plots[i] == "Eff_e"))
        {
          h->SetLineColor(colors[j]);
          h->SetLineStyle(linestyles[j]);
          h->SetLineWidth(2);   
        }

        else
        {
          if((names_files[j] == "Zprime_tata_350.root") || (names_files[j] == "Zprime_tata_1500.root"))
            {  
              h->SetLineColor(colors[j]);
              h->SetLineStyle(linestyles[j]);
              h->SetLineWidth(2);
            }

            else
            {
              h->SetFillColor(colors[j]);
              h->SetFillStyle(1001);
              h->SetLineColor(0);
            }
        }

        h->Scale(1.0/h->Integral());
    
        if((plots[i] == "N_Merged") || (plots[i] == "Eff_mu") || (plots[i] == "Eff_e"))
        {
          legend->AddEntry(h,names[j].c_str(),"l");   
        }

        else
        {
          if((names_files[j] == "Zprime_tata_350.root") || (names_files[j] == "Zprime_tata_1500.root"))
            {
              legend->AddEntry(h,names[j].c_str(),"l");
            }

            else
            {
              legend->AddEntry(h,names[j].c_str(),"f");
            }
        }
     
        legend->SetBorderSize(0);
        hs->Add(h);
         
    }   

        hs->Draw("NOSTACK HIST");
        hs->GetXaxis()->SetTitle(x_labels[i].c_str());
        hs->GetYaxis()->SetTitle("a.u.");
        legend->Draw();
        l->Add(c2);
        std::string filename = plots[i] + ".png";
        c2->SaveAs(filename.c_str());
    }

    TFile* Output = new TFile("joined.root", "RECREATE"); 
    l->Write();
    Output->Close();


    /*ofstream myfile_mu;
    ofstream myfile_e;

    myfile_mu.open ("Eff_mu.txt");
    myfile_e.open ("Eff_e.txt");

    for (int j=0; j<names.size(); j++)
    {    
        TFile f(names_files[j].c_str());
        TH1F *h_mu = (TH1F*)f.Get("Eff_mu");
        TH1F *h_e = (TH1F*)f.Get("Eff_e");
        myfile_mu << names_files_wo_ext[j] << " " << h_mu->GetBinContent(1) << " " << h_mu->GetBinContent(2) << " " << h_mu->GetBinContent(3) << " " << h_mu->GetBinContent(4) << " " << h_mu->GetBinContent(5) << " " << h_mu->GetBinContent(6) << "\n";
        myfile_e << names_files_wo_ext[j] << " " <<  h_e->GetBinContent(1) << " " << h_e->GetBinContent(2) << " " << h_e->GetBinContent(3) << " " << h_e->GetBinContent(4) << " " << h_e->GetBinContent(5) << " " << h_e->GetBinContent(6) << "\n";
    } 
    
    myfile_mu.close();
    myfile_e.close();
    */

       

    
}
