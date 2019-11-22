/*************************************************************************
 * Gainmatch.cpp script can be used to get slopes of Left vs right signal*
 * in the front strips of position sensitive SX3 detectors.              *
 *                                                                       *
 * Event by event pedestal substraction is done and LvR histograms are   *
 * drawn for different strips and sectors (12 dets *4 strips * 4 sectors)*
 * Highest energy alpha line is then separated and linear fit gives slope*
 * At the current stage, the average of 4 slopes (from 4 sectors) is     *
 * calculated for each strip and is stored in a dat file.                *
 * 									 *
 * To see if the linear fit is correct, Gains.root file is created inside*
 * the Output directory							 *
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
    //chain->Add(PathToFiles + "0448.root"); //SX3 upstream 2-4 (Currently, I am calibrating #4 with this file. Also, #5 is empty)
    //chain->Add(PathToFiles + "0449.root"); // SX3 upstream 6-8
    //chain->Add(PathToFiles + "0447.root"); //SX3 upstream 8-10 (Currently, I am calibrating #9 and #10 with this file. Also, #11 is empty)

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
	for (Int_t i=0; i<12; i++){ // Loop over detectors
		for (Int_t j=0; j<4; j++){ //Loop over the strips
			for (Int_t k=0; k<4; k++){ // Loop over the backs
				std::string nameSX3_LvR = Form("SX3_%i_Strip_%i_Back_%i_LvR",i,j,k);
				std::string nameSX3_EvP = Form("SX3_%i_Strip_%i_Back_%i_EvP",i,j,k);
				SX3_LvR[i][j][k] = new TH2F(nameSX3_LvR.c_str(), "SX3 detector Left vs Right", 1000,0,3000,1000,0,3000);
			}//End of loop over backs
		}//End of Loop over Fronts
	}//End of Loop over Detectors
	
	//Create Output File   
	TFile* outputFile = new TFile("/mnt/e/Analysis/SX3 Calibration/Output/Gains.root", "recreate");
	


    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) { //Event by event Loop
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

		
	   	
		//Loop over the Multiplicity
		for(Int_t j=0; j<SX3Mul; j++){
			
			//Pedestal Substraction
			Float_t SX3RawStripRight = SX3StripRightADC[j] - Pedestals[SX3Det[j]][SX3Strip[j]*2];
			Float_t SX3RawStripLeft = SX3StripLeftADC[j] - Pedestals[SX3Det[j]][SX3Strip[j]*2+1];

			Float_t SX3RawEnergy = SX3RawStripRight + SX3RawStripLeft;
		
			//Fill your histograms here
			if (SX3Upstream[j] && SX3RawEnergy > 2890.)SX3_LvR[SX3Det[j]][SX3Strip[j]][SX3Sector[j]]->Fill(SX3RawStripRight,SX3RawStripLeft); 
			
	
		}// End of the multiplicity Loop	
    }// End of event by event analysis
    
    //Defining linear fit to get the gains(slope)
    TF1* linefit = new TF1("linefit", "[0]*x + [1]",0,3000);
    outputFile->cd();

    //Cleaning the dat file to save the gains
    std::ofstream pfile;
    pfile.open("/mnt/e/Analysis/SX3 Calibration/SX3gains.dat", std::ofstream::out | std::ofstream::trunc);
    pfile.close();
    Float_t gain[4]={0};

    //Writing Histogram for 12*4*4 upstream SX3
   for (Int_t i=0; i<12; i++){ //Loop over detectors
   	for (Int_t j=0; j<4; j++){ //Loop over strips
		for (Int_t k=0; k<4;k++){ //Loop over Backs
			
			 //Linear fitting to get the Slopes for gain matching
			Int_t Bin = SX3_LvR[i][j][k]->FindFirstBinAbove(2,1);
			Float_t xvalue = ((TAxis*)SX3_LvR[i][j][k]->GetXaxis())->GetBinCenter(Bin);
			//std::cout << xvalue << std::endl;
			SX3_LvR[i][j][k]->Fit("linefit", "Q" , "", xvalue, xvalue+500);
			gain[k] = linefit->GetParameter (0);
			SX3_LvR[i][j][k]->Write();
		} //End of loop over backs

		Float_t average = ( gain[0] + gain[1] + gain[2] + gain[3] ) / 4.;
		//std::cout << gain[0] << '\t' << gain[1] << '\t' << gain[2] << '\t' << gain[3] << '\t' << average << std::endl;
		std::ofstream pfile;
		pfile.open("/mnt/e/Analysis/SX3 Calibration/SX3gains.dat", std::ofstream::app);
		pfile << average << std::endl;
		pfile.close();
	} //End of loop over strips
    } //End of loop over detectors

    outputFile->Close();
    
}














