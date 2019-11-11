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
	
	
	//Define Histograms here
	//Histograms for all Upstream 12 SX3 detectors with 4 strips and 4 backs
	TH2F *SX3_LvR[12][4][4]; //Left vs Right for 12*4*4 Upstream SX3
	TH2F *SX3_EvP[12][4][4]; //Energy vs Position for 12*4*4 Upstream SX3
	for (Int_t i=0; i<12; i++){ // Loop over detectors
		for (Int_t j=0; j<4; j++){ //Loop over the strips
			for (Int_t k=0; k<4; k++){ // Loop over the backs
				std::string nameSX3_LvR = Form("SX3_%i_Strip_%i_Back_%i_LvR",i,j,k);
				std::string nameSX3_EvP = Form("SX3_%i_Strip_%i_Back_%i_EvP",i,j,k);
				SX3_LvR[i][j][k] = new TH2F(nameSX3_LvR.c_str(), "SX3 detector Left vs Right", 1000,0,3000,1000,0,3000);
				SX3_EvP[i][j][k] = new TH2F(nameSX3_EvP.c_str(), "SX3 detector Energy vs Position", 2000,-2,2,2000,0,7000);
			}//End of loop over backs
		}//End of Loop over Fronts
	}//End of Loop over Detectors
	
	//Create Output File   
	TFile* outputFile = new TFile("/mnt/e/Analysis/SX3 Calibration/Output/Gains.root", "recreate");
	


    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

		
	   	
		//Loop over the Multiplicity
		for(Int_t j=0; j<SX3Mul; j++){
			
			//Pedestal Substraction
			Float_t SX3RawStripRight = SX3StripRightADC[j] - Pedestals[SX3Det[j]][SX3Strip[j]*2];
			Float_t SX3RawStripLeft = SX3StripLeftADC[j] - Pedestals[SX3Det[j]][SX3Strip[j]*2+1];

			Float_t SX3RawEnergy = SX3RawStripRight + SX3RawStripLeft;
		
			
			if (SX3Upstream[j] && SX3RawEnergy > 2670.)SX3_LvR[SX3Det[j]][SX3Strip[j]][SX3Sector[j]]->Fill(SX3RawStripRight,SX3RawStripLeft);
			
	
		}// End of the multiplicity Loop	
    }// End of event by event analysis
    
    //Defining linear fit to get the gains(slope)
    TF1* linefit = new TF1("linefit", "[0]*x + [1]",0,3000);
    outputFile->cd();

    //Cleaning the dat file to save the gains
    std::ofstream pfile;
    pfile.open("/mnt/e/Analysis/SX3 Calibration/SX3gains.dat", std::ofstream::out | std::ofstream::trunc);
    pfile.close();

    //Writing Histogram for 12*4*4 upstream SX3
   for (Int_t i=0; i<12; i++){
   	for (Int_t j=0; j<4; j++){
		for (Int_t k=0; k<4;k++){
			
			 //Linear fitting to get the Slopes for gain matching
			Int_t Bin = SX3_LvR[i][j][k]->FindFirstBinAbove(3,1);
			Float_t xvalue = ((TAxis*)SX3_LvR[i][j][k]->GetXaxis())->GetBinCenter(Bin);
			//std::cout << xvalue << std::endl;
			SX3_LvR[i][j][k]->Fit("linefit", "Q" , "", xvalue, xvalue+500);
			std::ofstream pfile;
			pfile.open("/mnt/e/Analysis/SX3 Calibration/SX3gains.dat", std::ofstream::app);
			pfile << linefit->GetParameter (0) << std::endl;
			pfile.close();
			SX3_LvR[i][j][k]->Write();
			//std::cout << linefit->GetParameter (0) << '\t' << linefit->GetParameter (1) << std::endl;
		}
	}
    }

    outputFile->Close();
}














