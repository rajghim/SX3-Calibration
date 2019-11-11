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
    chain->Add(PathToFiles + "0446.root");
    return chain;
}

void Analysis::Loop() {
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    TString CutPath = "/mnt/e/Goddess/Cut/";
    TString CutPrefix = "cut_p30dp";

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



	//Open the Energy Calibration file for Upstream SX3 
	std::ifstream SX3EnCalfile;
	SX3EnCalfile.open("/mnt/e/Analysis/SX3 Calibration/SX3EnCal.dat");
	Double_t SX3EnCalSlope[192] = {0};
	Double_t SX3EnCalIntercept[192] = {0}; 
	for (Int_t i=0; i<192; i++){
			SX3EnCalfile >> SX3EnCalSlope[i] >> SX3EnCalIntercept[i];
	}

	//Open the Position Calibration file for Upstream SX3
	std::ifstream SX3PosCalfile;
	SX3PosCalfile.open("/mnt/e/Analysis/SX3 Calibration/SX3PosCal.dat");
	Double_t SX3PosCalSlope[192] = {0};
	Double_t SX3PosCalIntercept[192] = {0}; 
	for (Int_t i=0; i<192; i++){
			SX3PosCalfile >> SX3PosCalSlope[i] >> SX3PosCalIntercept[i];
	}


	//Define Histograms here
	//Defining histograms detector by detector
	TH2F *SX3_EvA_TDC[12]; //Energy vs Angle Histograms for 12 upstream SX3 detectors
	TH1F *SX3_Ex_TDC[12]; //Excitation Energy Histograms for 12 upstream SX3 detectors
	TH2F *SX3_EvP[12][4]; //Energy vs Position Histograms for 12 upstream SX3 Detectors
 	for (int i=0; i<12; i++){
		for (int k=0; k<4; k++){
			//std::string nameSX3_EvA_TDC = Form("SX3u%i_EvA_TDC",i);
			//std::string nameSX3_Ex_TDC = Form("SX3u%i_Ex_TDC",i);
			std::string nameSX3_EvP = Form("SX3u%i_Strip%i_EvP",i,k);
			//SX3_EvA_TDC[i] = new TH2F(nameSX3_EvA_TDC.c_str(),"SX3 upstream Energy vs Angle with TDC",200,0,180,1000,0,10);
			//SX3_Ex_TDC[i] = new TH1F(nameSX3_Ex_TDC.c_str(),"SX3 upstream Excitation Energy with TDC",250,-2,8);
			SX3_EvP[i][k] = new TH2F(nameSX3_EvP.c_str(), "SX3 upstream Energy vs Position",1000,-10,100,1000,0,12000);
		}
	}
	


	TH2F *SX3EvA = new TH2F("SX3EvA", "SX3 Energy vs Angle",200,0,180,1000,0,10);
	TH2F *SX3EvA_TDC = new TH2F("SX3EvA_TDC","SX3 Energy vs Angle with TDC",1000,0,180,1000,0,10000);
	TH1F *SX3Ex = new TH1F("SX3Ex","SX3 Excitation Energy",1000,-2,10);
	TH1F *SX3qval_TDC = new TH1F("SX3qval_TDC","SX3 Q-Value",250,-10,15);
	TH1F *SX3Ex_TDC = new TH1F("SX3Ex_TDC","SX3 Excitation Energy with TDC",250,-2,8);
    	TH2F* hICdEE = new TH2F("icdEE", "Ionization Chamber dE vs E; E; dE", 500, 0, 4000, 500, 0, 4000);
    	TH2F* hICdEE_cut = new TH2F("icdEE_cut", "Ionization Chamber dE vs E (with IC Cut); E; dE", 500, 0, 4000, 500, 0, 4000);
	
	//Create Output File   
	TFile* outputFile = new TFile("/mnt/e/Analysis/SX3 Calibration/Output/Fronttest.root", "recreate");


    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
	
		//Loop over the Multiplicity
		for(Int_t j=0; j<SX3Mul; j++){


			//Pedestals Substracted
			Float_t SX3RawStripRight = SX3StripRightADC[j] - Pedestals[SX3Det[j]][SX3Strip[j]*2];
			Float_t SX3RawStripLeft = SX3StripLeftADC[j] - Pedestals[SX3Det[j]][SX3Strip[j]*2+1];
			
			//Gains applied
			Float_t RawStripLeft = SX3RawStripLeft; 
			Float_t RawStripRight = -1. * SX3RawStripRight * (Gains[(SX3Det[j]*12)+(SX3Strip[j]*4)+(SX3Sector[j])]);  
		
			//Energy and Position Calculations			
			Float_t RawEnergy = RawStripRight + RawStripLeft; //Gain Matched Energy (No Calibration though)
			Float_t Energy = RawEnergy * SX3EnCalSlope[(SX3Det[j]*12)+(SX3Strip[j]*4)+(SX3Sector[j])] + SX3EnCalIntercept[(SX3Det[j]*12)+(SX3Strip[j]*4)+(SX3Sector[j])]; //Energy Calibrated
			Float_t RawPosition = ((RawStripRight - RawStripLeft) / Energy ); // Gain matched Position( Energy has been Calibrated but No Calibration in Position)
			Float_t Position = (RawPosition * SX3PosCalSlope[(SX3Det[j]*12)+(SX3Strip[j]*4)+(SX3Sector[j])]) + SX3PosCalIntercept[(SX3Det[j]*12)+(SX3Strip[j]*4)+(SX3Sector[j])]; //Position Calibration applied

			//std::cout << RawPosition << '\t' << Position << std::endl;
		
			//Fill your histograms here
			if (SX3Upstream[j])SX3_EvP[SX3Det[j]][SX3Strip[j]]->Fill(Position,Energy);
			
			


			
		}//End of Loop Over Multiplicity	
    }// End of event by event analysis
  
    
    outputFile->cd();

    
    //Writing Histograms Detector by Detector
    for (Int_t i=0; i<12; i++){
	for (Int_t k=0; k<4; k++)
		SX3_EvP[i][k]->Write();
    }
    outputFile->Close();
}
