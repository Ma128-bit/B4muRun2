#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <filesystem>

using namespace RooFit;
namespace fs = std::filesystem;

void FitBsJPsiPhi_Mass(TString year="2022", TString label="") {
    // Aprire il file root contenente l'albero
    TFile *file = new TFile("ROOTFiles_"+label+"/AllControl"+year+".root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening the file" << std::endl;
        return;
    }

    // Ottenere l'albero dal file
    TTree *tree = (TTree*)file->Get("FinalTree");
    if (!tree) {
        std::cerr << "Error opening Tree" << std::endl;
        file->Close();
        return;
    }

    std::string label_str = label.Data();
    fs::path dir_path = "BsJPsiPhi_MassFit_"+label_str;

    // Creare la directory
    try {
        if (fs::create_directory(dir_path)) {
            std::cout << "Directory created successfully: " << dir_path << std::endl;
        } else {
            std::cout << "The directory already exists or could not be created." << std::endl;
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    
    //tree->Draw("Quadruplet_Mass_eq>>h1(52,5.0, 5.9)","isMC==0 && (isMedium[0]+isMedium[1]+isMedium[2]+isMedium[3]==4)");
    tree->Draw("Quadruplet_Mass_eq>>h1(30, 5.05, 5.7)","isMC==0");
    //tree->Draw("Quadruplet_Mass_eq>>h1(80, 4.8, 6.1)","isMC==0");
    TH1F *h1 = (TH1F*)gDirectory->Get("h1");
    
    tree->Draw("Quadruplet_Mass_eq>>h2(80, 5.05, 5.7)","isMC!=0");
    //tree->Draw("Quadruplet_Mass_eq>>h1(80, 4.8, 6.1)","isMC==0");
    TH1F *h2 = (TH1F*)gDirectory->Get("h2");

    RooRealVar x("NewMassEqation", "NewMassEqation", 5.05, 5.7);
    //RooRealVar isMC("isMC", "isMC", -10, 10); 
    //RooDataSet dataset("dataset", "dataset", RooArgSet(x, isMC), RooFit::Import(*tree));
    //RooDataSet *filteredDataset = (RooDataSet*)dataset.reduce("isMC==0");
    //RooDataSet data = *static_cast<RooDataSet*>((RooDataSet*)filteredDataset->reduce(RooArgSet(x)));

    x.setRange("R1", 5.05, 5.25);
    x.setRange("R2", 5.55, 5.7);
    x.setRange("RT", 5.05, 5.7);
    x.setRange("RS", 5.275, 5.45);
    //x.setBins(100);
    
    RooDataHist data("data", h1->GetTitle(), RooArgSet(x), Import(*h1, kFALSE));
    
    RooDataHist MC("MC","MC", RooArgSet(x), Import(*h2, kFALSE));
    
    // Creare il fondo
    RooRealVar gamma("#Gamma", "Gamma", -0.9, -10, 10);
    RooExponential exp_bkg("exp_bkg", "exp_bkg", x, gamma);
    exp_bkg.fitTo(data,Range("R1,R2"));

    RooRealVar mu("mu", "mu", 5.36, 4.50, 6.0);
    RooRealVar lambd("lambd", "lambd", 0.02, 0.001, 1.5);
    RooRealVar gamm("gamm", "gamm", 0.14, 0.01, 1.5);
    RooRealVar delta("delta", "delta", 1.45, 0.1, 10);

    //RooRealVar mu("mu", "mu", 5.46, 5.0, 5.90);
    //RooRealVar lambd("lambd", "lambd", 0.01, 100);
    //RooRealVar gamm("gamm", "gamm", 0.01, 100);
    //RooRealVar delta("delta", "delta", 0.01, 400);
    RooJohnson gauss_pdf("signal_Bs", "signal_Bs", x, mu, lambd, gamm, delta);
         
    // Creare la gaussiana
    //RooRealVar mean("mean", "Media gaussiana", 5.367, 5.33, 5.40);
    //RooRealVar sigma("sigma", "Deviazione standard gaussiana", 0.02, 0.01, 0.06);
    //RooGaussian gauss_pdf("gauss_pdf", "Signal Gaussian PDF", x, mean, sigma);
        
    // Creare il modello di fit combinando fondo e gaussiana
    RooRealVar nsig("nsig", "Numero di segnali", 60, 30, 1000);
    RooRealVar nbkg("nbkg", "Numero di background", h1->GetEntries(), 1, 2*h1->GetEntries());

    gauss_pdf.fitTo(MC, Save(true), Range("RS"));
    RooPlot *frame2 = x.frame();
    MC.plotOn(frame2);
    gauss_pdf.plotOn(frame2, LineStyle(kDashed), LineColor(kRed));  
    gauss_pdf.paramOn(frame2, Parameters(RooArgSet(mu, gamm, lambd, delta)), Layout(0.6,0.9,0.9));
    TCanvas *canvas2 = new TCanvas("canvas2", "Fit Result2", 900, 600);
    frame2->Draw();
    canvas2->SaveAs("BsJPsiPhi_MassFit_"+label+"/Fit_MC_"+year+".png");
    canvas2->Clear();
    canvas2->Delete();

    //mu.setConstant(kTRUE); 
    lambd.setConstant(kTRUE); 
    delta.setConstant(kTRUE); 
    gamm.setConstant(kTRUE); 

    RooAddPdf model("model", "Signal + Background", RooArgList(gauss_pdf,  exp_bkg), RooArgList(nsig, nbkg));
    
    RooFitResult *result = model.fitTo(data, Save(true), Range("RT"));
    RooPlot *frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame, Components(gauss_pdf), LineStyle(kDashed), LineColor(kRed));
    //model.paramOn(frame, Parameters(RooArgSet(nsig, nbkg, mean, sigma, gamma)), Layout(0.6,0.9,0.9));
    model.paramOn(frame, Parameters(RooArgSet(nsig, nbkg, mu, lambd, gamm ,delta , gamma)), Layout(0.6,0.9,0.9));
    model.plotOn(frame, Components(exp_bkg), LineStyle(kDashed), LineColor(kGreen));
    model.plotOn(frame);
    
    TCanvas *canvas = new TCanvas("canvas", "Fit Result", 900, 600);
    frame->Draw();
    canvas->SaveAs("BsJPsiPhi_MassFit_"+label+"/Fit_BsJPsiPhi_"+year+".png");

    canvas->Clear();
    canvas->Delete();
    // Chiudere il file
    file->Close();
}
