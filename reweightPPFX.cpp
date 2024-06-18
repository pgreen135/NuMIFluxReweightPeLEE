// macro to re-scale PPFX CV + systematic universes
// original ntuples --> g4_10_4 + geometry bug fixes
// pelee format

#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"

float getScalingFactor(float nu_e, float nu_l, TH2F &ppfx_ratio);

int main(int argc, char *argv[]) {

    // parse arguments
    if ( argc != 4 ) {
        std::cout << "Usage: ./Reweight INPUT_PELEE_NTUPLE_FILE OUPUT_FILE_NAME HORN_CURR" << std::endl;
        return 1;
    }

    std::string inputFileName( argv[1] );
    std::string outputFileName( argv[2] );
    std::string hornCurr( argv[3] );

    if (hornCurr != "FHC" && hornCurr != "RHC") {
        std::cout << "Invalid horn current: " << hornCurr << ". Accepted arguments: FHC, RHC." << std::endl;
        return 1;
    }

    std::cout << "Processing File: " << inputFileName << ", Horn Current Mode: " << hornCurr << std::endl;

    // read PPFX ratio histograms
    std::vector<std::string> ratioFileNames;
    ratioFileNames.push_back("weights/ratios_g4_10_4_ppfx_le_g3Chase_standardbin_set1.root");
    ratioFileNames.push_back("weights/ratios_g4_10_4_ppfx_le_g3Chase_standardbin_set2.root");
    ratioFileNames.push_back("weights/ratios_g4_10_4_ppfx_le_g3Chase_standardbin_set3.root");
    ratioFileNames.push_back("weights/ratios_g4_10_4_ppfx_le_g3Chase_standardbin_set4.root");
    ratioFileNames.push_back("weights/ratios_g4_10_4_ppfx_le_g3Chase_standardbin_set5.root");
    ratioFileNames.push_back("weights/ratios_g4_10_4_ppfx_le_g3Chase_standardbin_set6.root");

    std::vector<TFile*> f_ppfx_ratio;
    for (int i = 0; i < ratioFileNames.size(); i++) {
        f_ppfx_ratio.push_back(new TFile(ratioFileNames[i].c_str()));
    }

    // CV
    TH2F ppfx_cv_nue;
    TH2F ppfx_cv_nuebar;
    TH2F ppfx_cv_numu;
    TH2F ppfx_cv_numubar;
    
    if (hornCurr == "FHC") {
        ppfx_cv_nue = *((TH2F*)f_ppfx_ratio[0]->Get("ratio_cv_fhc_enu_length_nue"));
        ppfx_cv_nuebar = *((TH2F*)f_ppfx_ratio[0]->Get("ratio_cv_fhc_enu_length_nuebar"));
        ppfx_cv_numu = *((TH2F*)f_ppfx_ratio[0]->Get("ratio_cv_fhc_enu_length_numu"));
        ppfx_cv_numubar = *((TH2F*)f_ppfx_ratio[0]->Get("ratio_cv_fhc_enu_length_numubar"));
    }
    else {
        ppfx_cv_nue = *((TH2F*)f_ppfx_ratio[0]->Get("ratio_cv_rhc_enu_length_nue"));
        ppfx_cv_nuebar = *((TH2F*)f_ppfx_ratio[0]->Get("ratio_cv_rhc_enu_length_nuebar"));
        ppfx_cv_numu = *((TH2F*)f_ppfx_ratio[0]->Get("ratio_cv_rhc_enu_length_numu"));
        ppfx_cv_numubar = *((TH2F*)f_ppfx_ratio[0]->Get("ratio_cv_rhc_enu_length_numubar"));
    }

    // Universes
    std::vector<TH2F> ppfx_universes_nue;
    std::vector<TH2F> ppfx_universes_nuebar;
    std::vector<TH2F> ppfx_universes_numu;
    std::vector<TH2F> ppfx_universes_numubar;
    
    for (size_t f_idx = 0; f_idx < f_ppfx_ratio.size(); f_idx++) {
        for (size_t i = 0; i < 100; i++) {
            
            if (hornCurr == "FHC") {
                ppfx_universes_nue.push_back( *((TH2F*)f_ppfx_ratio[f_idx]->Get( ("univs/nue/ratio_univ" + std::to_string(i) + "_fhc_enu_length_nue").c_str())) );
                ppfx_universes_nuebar.push_back( *((TH2F*)f_ppfx_ratio[f_idx]->Get( ("univs/nuebar/ratio_univ" + std::to_string(i) + "_fhc_enu_length_nuebar").c_str())) );
                ppfx_universes_numu.push_back( *((TH2F*)f_ppfx_ratio[f_idx]->Get( ("univs/numu/ratio_univ" + std::to_string(i) + "_fhc_enu_length_numu").c_str())) );
                ppfx_universes_numubar.push_back( *((TH2F*)f_ppfx_ratio[f_idx]->Get( ("univs/numubar/ratio_univ" + std::to_string(i) + "_fhc_enu_length_numubar").c_str())) );
            }
            else {
                ppfx_universes_nue.push_back( *((TH2F*)f_ppfx_ratio[f_idx]->Get( ("univs/nue/ratio_univ" + std::to_string(i) + "_rhc_enu_length_nue").c_str())) );
                ppfx_universes_nuebar.push_back( *((TH2F*)f_ppfx_ratio[f_idx]->Get( ("univs/nuebar/ratio_univ" + std::to_string(i) + "_rhc_enu_length_nuebar").c_str())) );
                ppfx_universes_numu.push_back( *((TH2F*)f_ppfx_ratio[f_idx]->Get( ("univs/numu/ratio_univ" + std::to_string(i) + "_rhc_enu_length_numu").c_str())) );
                ppfx_universes_numubar.push_back( *((TH2F*)f_ppfx_ratio[f_idx]->Get( ("univs/numubar/ratio_univ" + std::to_string(i) + "_rhc_enu_length_numubar").c_str())) );
            }
        }
    }

    if (ppfx_universes_nue.size() != 600 || ppfx_universes_nuebar.size() != 600) {
        std::cout << "Warning unexpected number of PPFX universes: " << ppfx_universes_nue.size() << std::endl;
    }
    
    // input file
    TFile *f = new TFile(inputFileName.c_str());
    TTree *t_nu = (TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
    TTree *t_pot = (TTree*)f->Get("nuselection/SubRun");
    
    // neutrino truth information
    int nu_pdg;
    float nu_e;
    float nu_l;

    // neutrino parent information
    double mcflux_vx;
    double mcflux_vy;
    double mcflux_vz;

    // PPFX weights
    // CV
    float ppfx_cv;
    // Universes
    std::map<std::string, std::vector<double>> *mc_weights_map = nullptr;

    // branch addresses
    t_nu->SetBranchAddress("nu_pdg", &nu_pdg);
    t_nu->SetBranchAddress("nu_e", &nu_e);
    t_nu->SetBranchAddress("nu_l", &nu_l);

    t_nu->SetBranchAddress("par_decay_vx", &mcflux_vx);
    t_nu->SetBranchAddress("par_decay_vy", &mcflux_vy);
    t_nu->SetBranchAddress("par_decay_vz", &mcflux_vz);

    t_nu->SetBranchAddress("ppfx_cv", &ppfx_cv);
    t_nu->SetBranchAddress("weights", &mc_weights_map);

    // output file
    //std::string::size_type pos = inputFileName.find('.');
    //std::string filename_new =  inputFileName.substr(0, pos) + "_reweightedPPFX.root";
    std::string filename_new = outputFileName;    

    TFile *f_out = new TFile(filename_new.c_str(), "recreate");
    TTree *t_out_nu = t_nu->CloneTree(0);
    TTree *t_out_pot = t_pot->CloneTree();  // clone existing tree unchanged

    // loop over nu tree adjusting weights and filling
    int nEntries = t_nu->GetEntries();

    std::cout << "Number entries: " << nEntries << std::endl;
    
    for (int iEntry = 0; iEntry < nEntries; iEntry++) {

        t_nu->GetEntry(iEntry);

        if ( (iEntry != 0) && (nEntries >= 10) &&  (iEntry % (nEntries/10) == 0) ) {
            std::cout << Form("%i0%% Completed...\n", iEntry / (nEntries/10));
        }

        // remove events in the problem geometry region
        if ( (std::abs(mcflux_vx) > 58.42 || ((mcflux_vy > 35.56 && mcflux_vz > 350) || (mcflux_vy > 65 && mcflux_vz < 350))) && mcflux_vz < 4569.9 ) continue;

        // adjust weights here
        // CV weight
        float cv_scaling_factor = 1;

        if (nu_pdg == 12) cv_scaling_factor = getScalingFactor(nu_e, nu_l, ppfx_cv_nue);
        else if (nu_pdg == -12) cv_scaling_factor = getScalingFactor(nu_e, nu_l, ppfx_cv_nuebar);
        else if (nu_pdg == 14) cv_scaling_factor = getScalingFactor(nu_e, nu_l, ppfx_cv_numu);
        else if (nu_pdg == -14) cv_scaling_factor = getScalingFactor(nu_e, nu_l, ppfx_cv_numubar);
        else {
            std::cout << "Error. Invalid Neutrino Flavour: " << nu_pdg << std::endl;
            exit(1);
        }

        ppfx_cv = ppfx_cv*cv_scaling_factor;

        // Also adjust CV weight in universes map
        auto ppfx_cv_index = mc_weights_map->find("ppfx_cv_UBPPFXCV");
        if (ppfx_cv_index != mc_weights_map->end()) ppfx_cv_index->second[0] = ppfx_cv_index->second[0] * cv_scaling_factor;
        else {
            std::cout << "Error. Could not find ppfx cv universe in map." << std::endl;
            exit(1);
        }

        // Universes
        auto ppfx_universes_index = mc_weights_map->find("ppfx_all");

        std::vector<double> ppfx_universes;
        if (ppfx_universes_index != mc_weights_map->end()) ppfx_universes = ppfx_universes_index->second;
        else {
            std::cout << "Error. Could not find ppfx systematic universes." << std::endl;
            exit(1);
        }

        // reweight universes
        for (int idx = 0; idx < ppfx_universes_nue.size(); idx++) {
            float universe_scaling_factor = 1;

            if (nu_pdg == 12) universe_scaling_factor = getScalingFactor(nu_e, nu_l, ppfx_universes_nue[idx]);
            else if (nu_pdg == -12) universe_scaling_factor = getScalingFactor(nu_e, nu_l, ppfx_universes_nuebar[idx]);
            else if (nu_pdg == 14) universe_scaling_factor = getScalingFactor(nu_e, nu_l, ppfx_universes_numu[idx]);
            else if (nu_pdg == -14) universe_scaling_factor = getScalingFactor(nu_e, nu_l, ppfx_universes_numubar[idx]);
            else {
                std::cout << "Error. Invalid Neutrino Flavour: " << nu_pdg << std::endl;
                exit(1);
            }

            ppfx_universes[idx] = ppfx_universes[idx] * universe_scaling_factor;
        }

        // add updated universes into map
        ppfx_universes_index->second = ppfx_universes;

        // fill
        t_out_nu->Fill();
    }

    // close input file
    f->Close();

    // write outputs to new file
    f_out->cd();
    f_out->mkdir("nuselection");
    f_out->cd("/nuselection/");
    t_out_nu->Write("NeutrinoSelectionFilter");
    t_out_pot->Write("SubRun");
    f_out->Close();

    // clean up
    delete f;
    delete f_out;
    for (int i = 0; i < f_ppfx_ratio.size(); i++) delete f_ppfx_ratio[i];

    return 0;
}

float getScalingFactor(float nu_e, float nu_l, TH2F &ppfx_ratio) {

    int binx = ppfx_ratio.GetXaxis()->FindBin(nu_e);
    int biny = ppfx_ratio.GetYaxis()->FindBin(nu_l);
    float scaling_factor = ppfx_ratio.GetBinContent(binx,biny);
    return scaling_factor;
}
