#include <iostream>
#include <cmath>
#include <chrono>

#include <TRandom3.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObject.h>
#include <TFile.h>
#include <TAxis.h>
#include <TMath.h>


bool CosSquare(TRandom& Random, double& obs , double min , double max , double Omega = 1){ //accept and reject of cos^2
	while( true ){
		double RandomY =  Random.Rndm();
		double RandomX= max - (max - min) * Random.Rndm();
		if(RandomY <= TMath::Cos(Omega * RandomX) * TMath::Cos(Omega * RandomX)){
			obs = RandomX;
			return true;
		}
	}
	return false;
}

bool SinSquare(TRandom3& Random ,double& obs,  double min , double max , double Omega = 1){	//accept and reject of sin^2
	while( true ){
		double WertY =  Random.Rndm();
		double WertX = max - (max - min) * Random.Rndm();
		if(WertY <= TMath::Sin(Omega * WertX) * TMath::Sin(Omega * WertX)){
			obs = WertX;
			return true;
		}
	}
	return false;
}

bool Sin(TRandom3& Random , double& obs, double min , double max , double Omega = 1){
	if (TMath::Abs(Omega *( min - max ) ) < 2*TMath::Pi()){
		double Alpha = TMath::Cos(Omega * min);
		double Beta = TMath::Cos(Omega * max);
		obs = 1/Omega * TMath::ACos((Beta - Alpha )*Random.Rndm()-Alpha);
		return true;
	}else {
		while( true ){
			double WertY = 1 - 2 *Random.Rndm();
			double WertX = max - (max - min) * Random.Rndm();
			if(WertY <= TMath::Sin(WertX * Omega)){
				obs = WertX;
				return true;
			}
		}
	}	
	return false;
}

bool Cos(TRandom3& Random , double& obs, double min , double max , double Omega = 1){
	if (TMath::Abs(Omega *( min - max ) ) < 2*TMath::Pi()){
		double Alpha = TMath::Sin(Omega * min);
		double Beta = TMath::Sin(Omega * max);
		obs = 1/Omega * TMath::ASin((Beta - Alpha )*Random.Rndm()-Alpha);
		return true;
	}else {
		while( true ){
			double WertY = 1 - 2 * Random.Rndm();
			double WertX = max - (max - min) * Random.Rndm();
			if(WertY <= TMath::Cos(WertX * Omega)){
				obs = WertX;
				return true;
			}
		}
	}
	return false;
}

bool ConstTrafo(TRandom3& Random , double& obs , double min , double max){
	obs = (max - min) * Random.Rndm() + min;
	return true;
}

bool ExpCos(TRandom3& Random , double& obs , double Lambda  , double Omega  , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * Lambda );
	double Beta = TMath::Exp(-max * Lambda);
	while( true ){
		double WertX = -1/Lambda*TMath::Log(Alpha - (Alpha - Beta) *  Random.Rndm());
		double WertY = 1 - 2 * Random.Rndm();
		if(WertY <= TMath::Cos(Omega * WertX) ){
			obs = WertX;
			return true;
		}
	}
	return false;
}

bool ExpSin(TRandom3& Random , double& obs , double Lambda  , double Omega  , double min = 0  , double max = 18 ){ 
	double Alpha = TMath::Exp(-min * (Lambda));
	double Beta = TMath::Exp(-max * (Lambda));	
	while( true ){
		double WertX = -1/Lambda*TMath::Log( Alpha - (Alpha - Beta) * Random.Rndm());
		double WertY = 1 - 2 * Random.Rndm() ;
		if(WertY <= TMath::Sin(Omega * WertX) ){
			obs = WertX;
			return true;
		}
	}
	return false;
}

bool ExpCosH(TRandom3& Random , double& obs , double Lambda  , double Omega  , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * (Lambda + Omega));
	double Beta = TMath::Exp(-max * (Lambda + Omega ));	
	while( true ){
		double WertX = -1/(Lambda + Omega)*TMath::Log(Alpha - (Alpha - Beta) *  Random.Rndm());
		double WertY = Random.Rndm() * TMath::Exp(WertX);
		if(WertY <= TMath::Cos(Omega * WertX) * TMath::Exp(WertX) ){
			obs = WertX;
			return true;
		}
	}
	return false;
}

bool ExpSinH(TRandom3& Random , double& obs , double Lambda  , double Omega  , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * (Lambda + Omega));
	double Beta = TMath::Exp(-max * (Lambda + Omega) );	
	while( true ){
		double WertX = -1/(Lambda + Omega)*TMath::Log( Alpha - (Alpha - Beta) * Random.Rndm());
		double WertY = Random.Rndm() * TMath::Exp(WertX);
		if(WertY <= TMath::SinH(Omega * WertX) * TMath::Exp(WertX) ){
			obs = WertX;
			return true;
		}
	}
	return false;
}

bool Time_Trafo(TRandom3& Random, double& obs,
				  double a ,double b , double c , double d , 
				  double Gamma_s , double Delta_Gamma_s , double Delta_Mass_s , 
				  double min = 0 ,double max = 18)
{
	/*double AlphaExp = TMath::Exp(-Gamma_s * min);
	double BetaExp = TMath::Exp(-Gamma_s * max);
	double AlphaCosh = TMath::CosH(Delta_Gamma_s/2 * min);
	double BetaCosh = TMath::CosH(Delta_Gamma_s/2 * max);
	double AlphaSinh = TMath::SinH(Delta_Gamma_s/2 * min);
	double BetaSinh = TMath::SinH(Delta_Gamma_s/2 * max);
	double AlphaCos = TMath::Cos(Delta_Mass_s * min);
	double BetaCos = TMath::Cos(Delta_Mass_s * max);
	double AlphaSin = TMath::Sin(Delta_Mass_s * min);
	double BetaSin = TMath::Sin(Delta_Mass_s * max);

	double IntCosh = (AlphaExp *( Gamma_s *AlphaCosh + Delta_Gamma_s * AlphaSinh)-BetaExp*( Gamma_s/2 *BetaCosh + Delta_Gamma_s * BetaSinh ))/(pow(Gamma_s,2)-pow( Delta_Gamma_s , 2 ));
	double IntSinh = (AlphaExp *( Gamma_s *AlphaSinh + Delta_Gamma_s * AlphaCosh)-BetaExp*( Gamma_s/2 *BetaSinh + Delta_Gamma_s * BetaCosh ))/(pow(Gamma_s,2)-pow( Delta_Gamma_s , 2 ));
	double IntCos = (AlphaExp *( Gamma_s *AlphaCos - Delta_Mass_s * AlphaSin) + BetaExp*( Gamma_s/2 *BetaSin - Delta_Mass_s * BetaCos ))/(pow(Gamma_s,2)+pow( Delta_Gamma_s , 2 ));
	double IntSin = (AlphaExp *( Gamma_s *AlphaSin + Delta_Mass_s * AlphaCos) - BetaExp*( Gamma_s/2 *BetaSin + Delta_Mass_s * BetaCos ))/(pow(Gamma_s,2)+pow( Delta_Gamma_s , 2 ));
	*/
	double IntCosh = 1 ;
	double IntSinh = 1 ;
	double IntSin = 1 ;
	double IntCos = 1 ;
	double sum_Time = 	TMath::Abs(a * IntCosh) + TMath::Abs(b * IntSinh) + TMath::Abs(c * IntCosh) + TMath::Abs(d * IntSin);
	double TimeRndm = Random.Rndm();
	if(sum_Time  * TimeRndm < TMath::Abs(a * IntCosh)){
		ExpCosH( Random , obs , Gamma_s , Delta_Gamma_s/2);
		return true;
	}else if(sum_Time  * TimeRndm < TMath::Abs(a * IntCosh) + TMath::Abs(b * IntSinh)){
		ExpSinH( Random , obs , Gamma_s , Delta_Gamma_s/2);
		return true;
	}else if(sum_Time  * TimeRndm < TMath::Abs(a * IntCosh) + TMath::Abs(b * IntSinh) + TMath::Abs(c * IntCosh)){
		ExpCos( Random , obs , Gamma_s , Delta_Mass_s);
		return true;
	}else{	
		ExpSin( Random , obs , Gamma_s , Delta_Mass_s); 
		return true;
	}
	return false;
}

int main()
{
	auto start_time = std::chrono::high_resolution_clock::now();

	unsigned long seed =1234567; 
	double Daten = 1e6;
	TRandom3 Random(seed);

	//Const	
		double Delta_Mass_s = 17.757;
		double Gamma_s= 0.6629;
		double Delta_Gamma_s = 0.077;
		const double PI = TMath::Pi(); 

	// Configuration
		double Amplitude_0 = TMath::Sqrt(0.348);
		double Amplitude_Parallel = TMath::Sqrt(0.287);
		double Amplitude_Vertical = TMath::Sqrt(0.365);
		//double Amplitude_0 = 0.701;
		//double Amplitude_Parallel = 0.506;
		//double Amplitude_Vertical = 0.502;

		double Phase_0 = 0 ;
		double Phase_Parallel = 2.71;
		double Phase_Vertical = 2.39;

		double CP_Lambda = 1;
		double CP_Phi_sss = 0; // = 0.17

	// CP-Violation Parameter
		double CPV_C = (1- pow(CP_Lambda,2))/(1+pow(CP_Lambda,2));
		double CPV_S = -2 * CP_Lambda * TMath::Sin(CP_Phi_sss)/(1+pow(CP_Lambda,2));
		double CPV_D = -2 * CP_Lambda * TMath::Cos(CP_Phi_sss)/(1+pow(CP_Lambda,2));

	//Defination 
		double Phase_1 = Phase_Vertical - Phase_Parallel;
		double Phase_2 = Phase_Vertical - Phase_0;

		double Value_1 = 4 * Amplitude_0 * Amplitude_0;
		double Value_2 = 2 * Amplitude_Parallel * Amplitude_Parallel;
		double Value_3 = 2 * Amplitude_Vertical * Amplitude_Vertical;
		double Value_4 = 2 * Amplitude_Parallel * Amplitude_Vertical;
		double Value_5 = TMath::Sqrt(2)*Amplitude_Parallel * Amplitude_0;
		double Value_6 = TMath::Sqrt(2)*Amplitude_0 * Amplitude_Vertical;

		double Sum = Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6;

	//Creation of the Tree

   		TFile *root_file= new TFile("P2VV.root", "RECREATE");
   		TTree *tree = new TTree("P2VV" , "");
   		double Test;
   		tree->Branch("Test" , &Test , "Test/D" );
   		double Time;
   		tree->Branch("Time" , &Time , "Time/D" );
		double theta1;
   		tree->Branch("theta1",&theta1,"theta1/D");
   		double theta2;
   		tree->Branch("theta2",&theta2,"theta2/D");
   		double Phi;
   		tree->Branch("Phi" , &Phi , "Phi/D");
   		int cut;
   		tree->Branch("cut" , &cut , "cut/I");

	for (int i = 0 ; i < Daten ; ++i)
	{
		double RndmValue = Random.Rndm();

		if(( Sum * RndmValue < Value_1)){
			CosSquare(Random , theta1 , 0 , PI);
			CosSquare(Random , theta2 , 0 , PI);
			ConstTrafo(Random , Phi , -PI , PI);
			Time_Trafo(Random , Time ,
							  1 , CPV_D , CPV_C , -CPV_S ,
							  Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);
			cut = 1;
		}else if( Sum * RndmValue < Value_1 + Value_2){
			SinSquare(Random , theta1 , 0 , PI);
			SinSquare(Random , theta2 , 0 , PI);
			SinSquare(Random , Phi , -PI , PI);
			Time_Trafo(Random , Time ,
							  1 , CPV_D , CPV_C , -CPV_S ,
							  Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);
			cut = 2;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3){
			SinSquare(Random , theta1 , 0 , PI);
			SinSquare(Random , theta2 , 0 , PI);
			CosSquare(Random ,Phi, -PI , PI);
			Time_Trafo(Random , Time ,
							  1 , -CPV_D , CPV_C , CPV_S ,
							  Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);
			cut = 3;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4){
			SinSquare(Random , theta1 , 0 , PI, 1);
			SinSquare(Random , theta2 , 0 , PI , 1);
			Sin(Random , Phi , -PI , PI , -2);
			Time_Trafo(Random , Time ,
							  CPV_C * TMath::Sin(Phase_1) , CPV_S * TMath::Sin(Phase_1), TMath::Sin(Phase_1) , CPV_D*TMath::Sin(Phase_1) ,
							  Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);

			cut = 4;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5){
			Sin(Random , theta1 , 0 , PI , 2);
			Sin(Random , theta2 , 0 , PI , 2);
			Cos(Random , Phi , -PI , PI , 1);
			Time_Trafo(Random , Time ,
					   TMath::Sin(Phase_Parallel) ,CPV_D * TMath::Sin(Phase_Parallel), TMath::Sin(Phase_2) , -CPV_S *TMath::Sin(Phase_Parallel) ,
					   Gamma_s, Delta_Gamma_s , Delta_Mass_s ,  0 , 18);

			cut = 5;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6){
			Sin(Random , theta1 , 0 , PI , -2);
			Sin(Random , theta2 , 0 , PI , -2);
			Sin(Random , Phi , -PI , PI , -1);
			Time_Trafo(Random , Time ,
					   CPV_C * TMath::Sin(Phase_2) , CPV_S * TMath::Cos(Phase_2) , TMath::Sin(Phase_2) , CPV_D * TMath::Cos(Phase_2) ,
					   Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);

			cut = 6;
		}
		ExpCos( Random , Test , Gamma_s , Delta_Mass_s); 

		tree->Fill();

	}
   	root_file->Write("P2VV");

   	auto end_time  = std::chrono::high_resolution_clock::now();     
 	
   	double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();

   	std::cout << "Time: " << duration << " ms"<< std::endl;

   	root_file->Close();
	return 0;
}