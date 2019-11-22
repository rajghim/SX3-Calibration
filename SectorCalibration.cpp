/*************************************************************************
 * SectorCalibration.cpp script can be used to alpha calibrate energies  *
 * for the back sectors of SX3 detectors.             			 *
 *                                                                       *
 * Event by event pedestal substraction is done and secotr Energy is     *
 * projected onto 1D Energy histograms (12 dets*4 strips*4 sectors)	 *
 * Energy calibration (for 12 dets*4 strips*4sectors) is done by gaussian*
 * fit of each histograms. The channel numbers are then saved in         *
 * SX3SectorEnCalChannels.dat file.                			 *
 * 									 *
 * To see if the gaussian fit is correct, SectorEnCal.root file is       *
 * created in the Output directory. 					 *
 *************************************************************************/

#define Analysis_cxx

#include "Analysis.h"


#include <string>
#include <fstream>

TChain* MakeChain();

int main() {
    TChain *chain = MakeChain();

    Analysis t(chain);
    t.Loop();
}

TChain* MakeChain() {
    auto *chain = new TChain("data");
    TString PathToFiles = "/mnt/e/goddessSort-master/Output/Run";

    chain->Add(PathToFiles + "0446.root"); //SX3 upstream 0-3
    //chain->Add(PathToFiles + "0448.root"); //SX3 upstream 2-4 (detector 5 is empty. I am calibrating #4 with this file)
    //chain->Add(PathToFiles + "0449.root"); // SX3 upstream 6-8
    //chain->Add(PathToFiles + "0447.root"); //SX3 upstream 8-10 (I am calibrating #9 and #10 with this file. Also, 11 is empty)



    return chain;
}

void Analysis::Loop() {
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    int prevRunNumber = -1;



	//Open the pedestals File 
	std::ifstream file;
	file.open("Sectorpedestals.dat");
	Double_t SectorPedestals[192] = {0};
	for (Int_t i = 0; i<192;i++){
		file >> SectorPedestals[i];	
	}



	//Define Histograms here
	//Histograms without the gains applied
	TH1F *SX3_SectorEnCal[12][4][4]; // Raw Energy histograms for Sector Energy Calibration
	for (Int_t i=0; i<12; i++){ // Loop over detectors
		for (Int_t j=0; j<4; j++){ //Loop over the strips
			for (Int_t k=0; k<4; k++){ // Loop over the backs
				std::string nameSX3_SectorEnCal = Form("SX3_SectorEnCal_Det%i_Strip_%i_Back%i",i,j,k);
				SX3_SectorEnCal[i][j][k] = new TH1F(nameSX3_SectorEnCal.c_str(), "SX3 Energy Spectrum for Sector Energy Calibration",1000,1000,4000);
			}//End of loop over backs
		}//End of Loop over Fronts
	}//End of Loop over Detectors


	//Create Output File   
	TFile* outputFile = new TFile("/mnt/e/Analysis/SX3 Calibration/Output/SectorEnCal.root", "recreate");


    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

		
	   	
		//Loop over the Multiplicity
		for(Int_t j=0; j<SX3Mul; j++){
			
			//Pedestals Substracted
			Float_t SX3RawSector = SX3SectorADC[j] - SectorPedestals[(SX3Det[j]*16) + (SX3Strip[j]*4) + (SX3Sector[j])];

			//Filling the histograms here
			if (SX3Upstream[j])SX3_SectorEnCal[SX3Det[j]][SX3Strip[j]][SX3Sector[j]]->Fill(SX3RawSector);

		
	
		}// End of the multiplicity Loop	
    }// End of event by event analysis

    //Cleaning the dat file to save the gains
   std::ofstream pfile;
    pfile.open("/mnt/e/Analysis/SX3 Calibration/SX3SectorEnCalChannels.dat", std::ofstream::out | std::ofstream::trunc);
    pfile.close();
    
    //Define your fit function here
    TF1* fit_gaus = new TF1 ("fit_gaus","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))",1000,4000);
    
    Int_t npeaks = 5;
    TSpectrum *s = new TSpectrum(npeaks);
    Float_t xp[5] = {0};
	
    outputFile->cd();

    for (Int_t i=0; i<12; i++){
   	for (Int_t j=0; j<4; j++){
		for (Int_t k=0; k<4;k++){

			Int_t nfound = s->Search(SX3_SectorEnCal[i][j][k],2.5,"",0.4);
			Double_t *xpeaks = s->GetPositionX();
			for (Int_t p=0; p<nfound; p++){
 				xp[p] = xpeaks[p];
			}
			std::sort(xp, xp+5);
			Float_t xmin = std::min(std::min(std::min(std::min(xp[0],xp[1]),xp[2]),xp[3]),xp[4]);
			Float_t xmax = std::max(std::max(std::max(std::max(xp[0],xp[1]),xp[2]),xp[3]),xp[4]);
			fit_gaus->SetParLimits (1,xp[0]-5,xp[0]+5);
			fit_gaus->SetParLimits(4,xp[1]-5,xp[1]+5);
			fit_gaus->SetParLimits(6,xp[2]-5,xp[2]+5);
			fit_gaus->SetParLimits(8,xp[3]-5,xp[3]+5);
			fit_gaus->SetParLimits(10,xp[4]-5,xp[4]+5);
			fit_gaus->SetParLimits(2,1,50);	
			
			SX3_SectorEnCal[i][j][k]->Fit("fit_gaus","Q","",xmin-50,xmax+50);
			std::ofstream pfile;
			pfile.open("/mnt/e/Analysis/SX3 Calibration/SX3SectorEnCalChannels.dat", std::ofstream::app);
			pfile << i << '\t' << j << '\t' << k << '\t' << fit_gaus->GetParameter (1) << '\t' << fit_gaus->GetParameter (4) << '\t' << fit_gaus->GetParameter (6) << '\t' << fit_gaus->GetParameter (8) << '\t' << fit_gaus->GetParameter (10) << std::endl;
			pfile.close();
	
			SX3_SectorEnCal[i][j][k]->Write();   
  		}
	}
    }
      

    outputFile->Close();
   
}















