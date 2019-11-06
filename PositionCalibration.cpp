

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

    chain->Add(PathToFiles + "0000.root");



    return chain;
}

void Analysis::Loop() {
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    int prevRunNumber = -1;



	//Open the pedestals File 
	std::ifstream file;
	file.open("SX3pedestals.dat");
	Double_t Pedestals[12][8] = {0};
	for (Int_t i = 0; i<12;i++){
		for (Int_t k=0; k<8; k++){
			file >> Pedestals[i][k];	
			//std::cout << Pedestals[10][6];
		}
	}


	//Open the gains file
	std::ifstream gainfile;
	gainfile.open("/mnt/e/Analysis/SX3 Calibration/SX3gains.dat");
	Double_t Gains[192] = {0}; 
	for (Int_t i=0; i<192; i++){
			gainfile >>Gains[i];
			//std::cout << Gains[2] << std::endl;
	}

	//Open the Energy Calibration file
	std::ifstream EnCalfile;
	EnCalfile.open("/mnt/e/Analysis/SX3 Calibration/SX3EnCal.dat");
	Double_t EnCalSlope[192] = {0};
	Double_t EnCalIntercept[192] = {0}; 
	for (Int_t i=0; i<192; i++){
			EnCalfile >>EnCalSlope[i] >> EnCalIntercept[i];
			//std::cout << EnCalSlope[2] << " " << EnCalIntercept[3] << std::endl;
	}
	
	

	//Define Histograms here
	//Histograms without the gains applied
	TH1F *SX3_PosCal[12][4][4]; // E = L+R Histogram for the Calibration of the Energy
	for (Int_t i=0; i<12; i++){ // Loop over detectors
		for (Int_t j=0; j<4; j++){ //Loop over the strips
			for (Int_t k=0; k<4; k++){ // Loop over the backs
				std::string nameSX3_PosCal = Form("SX3_PosCal_Det%i_Strip_%i_Back%i",i,j,k);
				SX3_PosCal[i][j][k] = new TH1F(nameSX3_PosCal.c_str(), "SX3 Position Spectrum for Front Side Position Calibration",1000,-1,1);
			}//End of loop over backs
		}//End of Loop over Fronts
	}//End of Loop over Detectors*/


	//Create Output File   
	TFile* outputFile = new TFile("/mnt/e/Analysis/SX3 Calibration/Output/PosCaloutput.root", "recreate");


    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

		
	   	
		//Loop over the Multiplicity
		for(Int_t j=0; j<SX3Mul; j++){
			
			//Without gains (Pedestals Substracted)
			Float_t SX3RawStripRight = SX3StripRightADC[j] - Pedestals[SX3Det[j]][SX3Strip[j]*2];
			Float_t SX3RawStripLeft = SX3StripLeftADC[j] - Pedestals[SX3Det[j]][SX3Strip[j]*2+1];
		
			//Gain Adjustement and Calibration Applied	
			Float_t RawStripLeft = SX3RawStripLeft; 
			Float_t RawStripRight = -1. * SX3RawStripRight * (Gains[(SX3Det[j]*12)+(SX3Strip[j]*4)+(SX3Sector[j])]); //Gains applied 

			Float_t RawEnergy = RawStripRight + RawStripLeft; // Gain matched Energy (No calibration though)
			Float_t Energy = RawEnergy * EnCalSlope[(SX3Det[j]*12)+(SX3Strip[j]*4)+(SX3Sector[j])] + EnCalIntercept[(SX3Det[j]*12)+(SX3Strip[j]*4)+(SX3Sector[j])]; //Energy Calibrated
			Float_t RawPosition = ((RawStripRight - RawStripLeft) / Energy ); // Gain matched Position( Energy has been Calibrated but No Calibration on Left and Right though)


			//Filling the histograms here
			if (SX3Upstream[j])SX3_PosCal[SX3Det[j]][SX3Strip[j]][SX3Sector[j]]->Fill(RawPosition);
		
	
		}// End of the multiplicity Loop	
    }// End of event by event analysis
   
    
    outputFile->cd();

    for (Int_t i=0; i<12; i++){
   	for (Int_t j=0; j<4; j++){
		for (Int_t k=0; k<4;k++){
			SX3_PosCal[i][j][k]->Write();   
  		}
	}
    }
      

    outputFile->Close();
}













