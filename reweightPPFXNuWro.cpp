// macro to re-scale NuWro fake data using alternative reweighting scheme 
// original ntuples, without PPFX CV --> g4_10_4 + geometry bug fixes
// only PPFX CV, PPFX universes not present for NuWro
// pelee format

#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TRotation.h"

float getScalingFactor(float nu_e, float nu_angle, TH2F &ratio);
float GetNuMIAngle(float px, float py, float pz, std::string direction);

int main(int argc, char *argv[]) {

    // parse arguments
    if ( argc != 4 ) {
        std::cout << "Usage: ./Reweight INPUT_NUWRO_PELEE_NTUPLE_FILE OUPUT_FILE_NAME HORN_CURR" << std::endl;
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

    TFile *PPFXCVFile = new TFile("weights/ratios_g4_10_4_cvonly_fornuwrofakedata.root"); // old flux without ppfx --> new flux with ppfx, E/Theta binning

    TH2F ppfx_cv_nue;
    TH2F ppfx_cv_nuebar;
    TH2F ppfx_cv_numu;
    TH2F ppfx_cv_numubar;
    
    if (hornCurr == "FHC") {
        ppfx_cv_nue = *((TH2F*)PPFXCVFile->Get("ratio_cv_fhc_enu_angle_nue"));
        ppfx_cv_nuebar = *((TH2F*)PPFXCVFile->Get("ratio_cv_fhc_enu_angle_nuebar"));
        ppfx_cv_numu = *((TH2F*)PPFXCVFile->Get("ratio_cv_fhc_enu_angle_numu"));
        ppfx_cv_numubar = *((TH2F*)PPFXCVFile->Get("ratio_cv_fhc_enu_angle_numubar"));
    }
    else {
        ppfx_cv_nue = *((TH2F*)PPFXCVFile->Get("ratio_cv_rhc_enu_angle_nue"));
        ppfx_cv_nuebar = *((TH2F*)PPFXCVFile->Get("ratio_cv_rhc_enu_angle_nuebar"));
        ppfx_cv_numu = *((TH2F*)PPFXCVFile->Get("ratio_cv_rhc_enu_angle_numu"));
        ppfx_cv_numubar = *((TH2F*)PPFXCVFile->Get("ratio_cv_rhc_enu_angle_numubar"));
    }
    
    // input file
    TFile *f = new TFile(inputFileName.c_str());
    TTree *t_nu = (TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
    TTree *t_pot = (TTree*)f->Get("nuselection/SubRun");
    
    // neutrino truth information
    int nu_pdg;
    float nu_e;

    float true_nu_px;
    float true_nu_py;
    float true_nu_pz;

    // PPFX weights
    // CV
    float ppfx_cv;

    // branch addresses
    t_nu->SetBranchAddress("nu_pdg", &nu_pdg);
    t_nu->SetBranchAddress("nu_e", &nu_e);

    t_nu->SetBranchAddress("true_nu_px", &true_nu_px);
    t_nu->SetBranchAddress("true_nu_py", &true_nu_py);
    t_nu->SetBranchAddress("true_nu_pz", &true_nu_pz);

    t_nu->SetBranchAddress("ppfx_cv", &ppfx_cv);

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

        // calculate true nu angle from numi beamline 
        float nu_angle = GetNuMIAngle(true_nu_px, true_nu_py, true_nu_pz, "beam");

        // adjust weights here
        // CV weight
        float cv_scaling_factor = 1;

        if (nu_pdg == 12) cv_scaling_factor = getScalingFactor(nu_e, nu_angle, ppfx_cv_nue);
        else if (nu_pdg == -12) cv_scaling_factor = getScalingFactor(nu_e, nu_angle, ppfx_cv_nuebar);
        else if (nu_pdg == 14) cv_scaling_factor = getScalingFactor(nu_e, nu_angle, ppfx_cv_numu);
        else if (nu_pdg == -14) cv_scaling_factor = getScalingFactor(nu_e, nu_angle, ppfx_cv_numubar);
        else {
            std::cout << "Error. Invalid Neutrino Flavour: " << nu_pdg << std::endl;
            exit(1);
        }

        ppfx_cv = cv_scaling_factor; // scaling factor from unweighted old CV == 1 --> ppfx weighted new CV

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

    return 0;
}

float getScalingFactor(float nu_e, float nu_angle, TH2F &ratio) {

    int binx = ratio.GetXaxis()->FindBin(nu_e);
    int biny = ratio.GetYaxis()->FindBin(nu_angle);
    float scaling_factor = ratio.GetBinContent(binx,biny);
    return scaling_factor;
}

// NuMI angle calculation [Krishan]
float GetNuMIAngle(float px, float py, float pz, std::string direction) {

    // Variables
    TRotation RotDet2Beam;             // Rotations
    TVector3  detxyz, BeamCoords;      // Translations
    std::vector<double> rotmatrix;     // Inputs

    // input detector coordinates to translate
    detxyz = {px, py, pz};     

    // From beam to detector rotation matrix
    rotmatrix = {
        0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021,
        4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359,
        -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291 };

    // Return the TRotation
    TVector3 newX, newY, newZ;
    newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
    newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
    newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

    RotDet2Beam.RotateAxes(newX, newY, newZ); // Return the TRotation now det to beam
    // RotDet2Beam.Invert(); // Invert back to the beam to det

    // Rotate to beam coords
    BeamCoords = RotDet2Beam * detxyz;

    TVector3 beamdir = {0 , 0 , 1};;
    
    // Get the angle wrt to the beam
    if (direction == "beam") beamdir = {0 , 0 , 1};
    
    // Get the angle wrt to the target to detector direction
    else if (direction == "target") {
        beamdir = {5502, 7259, 67270};
        beamdir = beamdir.Unit(); // Get the direction
    }
    else {
        std::cout << "Warning unknown angle type specified, you should check this" << std::endl;
    }
    
    double angle = BeamCoords.Angle(beamdir) * 180 / 3.1415926;

    // Create vectors to get the angle in the yz and xz planes
    TVector3 BeamCoords_yz = { 0, BeamCoords.Y(), BeamCoords.Z() }; // Angle upwards
    TVector3 BeamCoords_xz = { BeamCoords.X(), 0, BeamCoords.Z() }; // Angle across

    return angle;
}
