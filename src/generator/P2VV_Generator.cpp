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


double CosSquare(TRandom& Random, double min , double max , double Omega = 1){ //accept and reject of cos^2
	while( true ){
		double RandomY =  Random.Rndm();
		double RandomX= max - (max - min) * Random.Rndm();
		if(RandomY <= TMath::Cos(Omega * RandomX) * TMath::Cos(Omega * RandomX)){
			return RandomX;
		}
	}
	return 4;
}

double SinSquare(TRandom3& Random , double min , double max , double Omega = 1){	//accept and reject of sin^2
	while( true ){
		double WertY =  Random.Rndm();
		double WertX = max - (max - min) * Random.Rndm();
		if(WertY <= TMath::Sin(Omega * WertX) * TMath::Sin(Omega * WertX)){
			return WertX;
		}
	}
	return 4;
}

double Sin(TRandom3& Random , double min , double max , double Omega = 1){
	if (TMath::Abs(Omega *( min - max ) ) < 2*TMath::Pi()){
		double Alpha = TMath::Cos(Omega * min);
		double Beta = TMath::Cos(Omega * max);
		return 1/Omega * TMath::ACos((Beta - Alpha )*Random.Rndm()-Alpha);
	}else {
		while( true ){
			double WertY = 1 - 2*Random.Rndm();
			double WertX = max - (max - min) * Random.Rndm();
			if(WertY <= TMath::Sin(WertX * Omega)){
				return WertX;
			}
		}
	}	
	return 4;
}

double Cos(TRandom3& Random , double min , double max , double Omega = 1){
	if (TMath::Abs(Omega *( min - max ) ) < 2*TMath::Pi()){
		double Alpha = TMath::Sin(Omega * min);
		double Beta = TMath::Sin(Omega * max);
		return 1/Omega * TMath::ASin((Beta - Alpha )*Random.Rndm()-Alpha);
	}else {
		while( true ){
			double WertY = 1 - 2*Random.Rndm();
			double WertX = max - (max - min) * Random.Rndm();
			if(WertY <= TMath::Cos(WertX * Omega)){
				return WertX;
			}
		}
	}
	return 4;
}

double ConstTrafo(TRandom3& Random , double min , double max){
	return (max - min) * Random.Rndm() + min;
}

double ExpCos(TRandom3& Random , double Lambda  , double Omega  , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * Lambda );
	double Beta = TMath::Exp(-max * Lambda);
	while( true ){
		double WertX = -1/Lambda*TMath::Log(Alpha - (Alpha - Beta) *  Random.Rndm());
		double WertY =1 -2 * Random.Rndm();
		if(WertY <= TMath::Cos(Omega * WertX) ){
			return WertX;
		}
	}
	return 4;
}

double ExpSin(TRandom3& Random , double Lambda  , double Omega  , double min = 0  , double max = 18 ){ 
	double Alpha = TMath::Exp(-min * (Lambda));
	double Beta = TMath::Exp(-max * (Lambda));	
	while( true ){
		double WertX = -1/Lambda*TMath::Log( Alpha - (Alpha - Beta) * Random.Rndm());
		double WertY = 1 - 2 * Random.Rndm() ;
		if(WertY <= TMath::Sin(Omega * WertX) ){
			return WertX;
		}
	}
	return 4;
}

double ExpCosH(TRandom3& Random , double Lambda  , double Omega  , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * (Lambda + Omega));
	double Beta = TMath::Exp(-max * (Lambda + Omega ));	
	while( true ){
		double WertX = -1/(Lambda + Omega)*TMath::Log(Alpha - (Alpha - Beta) *  Random.Rndm());
		double WertY = Random.Rndm() * TMath::Exp(WertX);
		if(WertY <= TMath::Cos(Omega * WertX) * TMath::Exp(WertX) ){
			return WertX;
		}
	}
	return 4;
}

double ExpSinH(TRandom3& Random , double Lambda  , double Omega  , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * (Lambda + Omega));
	double Beta = TMath::Exp(-max * (Lambda + Omega) );	
	while( true ){
		double WertX = -1/(Lambda + Omega)*TMath::Log( Alpha - (Alpha - Beta) * Random.Rndm());
		double WertY = Random.Rndm() * TMath::Exp(WertX);
		if(WertY <= TMath::SinH(Omega * WertX) * TMath::Exp(WertX) ){
			return WertX;
		}
	}
	return 4;
}

double Time_Trafo(TRandom3& Random, 
				  double a ,double b , double c , double d , 
				  double Gamma_s , double Delta_Gamma_s , double Delta_Mass_s , 
				  double min = 0 ,double max = 18)
{
	double sum_Time = 	TMath::Abs(a) + TMath::Abs(b) + TMath::Abs(c) + TMath::Abs(d);
	double TimeRndm = Random.Rndm();
	if(sum_Time  * TimeRndm < TMath::Abs(a)){
		return ExpCosH( Random , Gamma_s , Delta_Gamma_s/2);
	}else if(sum_Time  * TimeRndm < TMath::Abs(a) + TMath::Abs(b)){
		return ExpSinH( Random , Gamma_s , Delta_Gamma_s/2);
	}else if(sum_Time  * TimeRndm < TMath::Abs(a) + TMath::Abs(b) + TMath::Abs(c)){
		return ExpCos( Random , Gamma_s , Delta_Mass_s);
	}else{	
		return ExpSin( Random , Gamma_s , Delta_Mass_s); 
	}
	return 0;
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

   		TFile *root_file= new TFile("~/Daten/Bachelor/P2VV.root", "RECREATE");
   		TTree *tree = new TTree("P2VV");
   		//double Test;
   		//tree->Branch("Test" , &Test , "Test/D" );
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
			theta1 = CosSquare(Random , 0 , PI);
			theta2 = CosSquare(Random , 0 , PI);
			Phi = ConstTrafo(Random , -PI , PI);
			Time = Time_Trafo(Random ,
							  1 , CPV_D , CPV_C , -CPV_S ,
							  Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);
			cut = 1;
		}else if( Sum * RndmValue < Value_1 + Value_2){
			theta1 = SinSquare(Random , 0 , PI);
			theta2 = SinSquare(Random , 0 , PI);
			Phi = SinSquare(Random , -PI , PI);
			Time = Time_Trafo(Random ,
							  1 , CPV_D , CPV_C , -CPV_S ,
							  Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);

			cut = 2;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3){
			theta1 = SinSquare(Random , 0 , PI);
			theta2 = SinSquare(Random , 0 , PI);
			Phi = CosSquare(Random , -PI , PI);
			Time = Time_Trafo(Random ,
							  1 , -CPV_D , CPV_C , CPV_S ,
							  Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);
			cut = 3;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4){
			theta1 = SinSquare(Random , 0 , PI);
			theta2 = SinSquare(Random , 0 , PI);
			Phi = Sin(Random , -PI , PI , -2);
			Time = Time_Trafo(Random ,
							  CPV_C * TMath::Sin(Phase_1) , CPV_S * TMath::Sin(Phase_1), TMath::Sin(Phase_1) , CPV_D*TMath::Sin(Phase_1) ,
							  Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);

			cut = 4;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5){
			theta1 = Sin(Random , 0 , PI , 2);
			theta2 = Sin(Random , 0 , PI , 2);
			Phi = Cos(Random , -PI , PI , 1);
			Time = Time_Trafo(Random ,
							  TMath::Sin(Phase_Parallel) ,CPV_D * TMath::Sin(Phase_Parallel), TMath::Sin(Phase_2) , -CPV_S *TMath::Sin(Phase_Parallel) ,
							  Gamma_s, Delta_Gamma_s , Delta_Mass_s ,  0 , 18);

			cut = 5;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6){
			theta1 = Sin(Random , 0 , PI , -2);
			theta2 = Sin(Random , 0 , PI , -2);
			Phi = Sin(Random , -PI , PI , -1);
			Time = Time_Trafo(Random ,
							  CPV_C * TMath::Sin(Phase_2) , CPV_S * TMath::Cos(Phase_2) , TMath::Sin(Phase_2) , CPV_D * TMath::Cos(Phase_2) ,
							   Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);

			cut = 6;
		}
		//Test = ExpCos( Random , 1 , 1 , 0 ,PI);

		tree->Fill();

	}
   	root_file->Write("P2VV");

   	auto end_time  = std::chrono::high_resolution_clock::now();     
 	
   	double zeit = std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();

   	std::cout << "Zeit: " << zeit << " ms"<< std::endl;

   	root_file->Close();
	return 0;
}