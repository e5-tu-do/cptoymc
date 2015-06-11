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


bool CosSquare(TRandom& Random, double& obs , double min , double max , double sign = 1,double Omega = 1){ //accept and reject of cos^2
	while( true ){
		double RandomY =  sign * Random.Rndm();
		double RandomX= max - (max - min) * Random.Rndm();
		if(RandomY <= sign * TMath::Cos(Omega * RandomX) * TMath::Cos(Omega * RandomX)){
			obs = RandomX;
			return true;
		}
	}
	return false;
}

bool SinSquare(TRandom3& Random ,double& obs,  double min , double max , double sign = 1 , double Omega = 1){	//accept and reject of sin^2
	while( true ){
		double RandomY = sign * Random.Rndm();
		double RandomX = max - (max - min) * Random.Rndm();
		if(RandomY <= sign * TMath::Sin(Omega * RandomX) * TMath::Sin(Omega * RandomX)){
			obs = RandomX;
			return true;
		}
	}
	return false;
}

bool Sin_Trafo(TRandom3& Random , double& obs, double min , double max , double sign = 1, double Omega = 1){
	/*if (TMath::Abs(Omega *( min - max ) ) < 2*TMath::Pi()){
		double Alpha = TMath::Cos(Omega * min);
		double Beta = TMath::Cos(Omega * max);
		obs = 1/Omega * TMath::ACos((Beta - Alpha )*Random.Rndm()-Alpha);
		return true;
	}else {*/
		while( true ){
			double RandomY = sign * (1 - 2 *Random.Rndm());
			double RandomX = max - (max - min) * Random.Rndm();
			if( RandomY <= sign * TMath::Sin(  RandomX * Omega)){
				obs = RandomX;
				return true;
			}
		}
	return false;
}

bool Cos_Trafo(TRandom3& Random , double& obs, double min , double max , double sign = 1 , double Omega = 1){
	/*if (TMath::Abs(Omega *( min - max ) ) < 2*TMath::Pi()){
		double Alpha = TMath::Sin(Omega * min);
		double Beta = TMath::Sin(Omega * max);
		obs = 1/Omega * TMath::ASin((Beta - Alpha )*Random.Rndm()-Alpha);
		return true;
	}else {*/
		while( true ){
			double RandomY = sign * (1 - 2 * Random.Rndm());
			double RandomX = max - (max - min) * Random.Rndm();
			if(  RandomY <=   sign *TMath::Cos(RandomX * Omega)){
				obs = RandomX;
				return true;
			}
		}
	return false;
}

bool ConstTrafo(TRandom3& Random , double& obs , double min , double max){
	obs = (max - min) * Random.Rndm() + min;
	return true;
}

bool ExpCos(TRandom3& Random , double& obs , double Lambda  , double Omega , int sign , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * Lambda );
	double Beta = TMath::Exp(-max * Lambda);
	while( true ){
		double RandomX = -1/Lambda*TMath::Log(Alpha - (Alpha - Beta) *  Random.Rndm());
		double RandomY = 1 - 2 * Random.Rndm();
		if(RandomY <= sign * TMath::Cos(Omega * RandomX) ){
			obs = RandomX;
			return true;
		}
	}
	return false;
}

bool ExpSin(TRandom3& Random , double& obs , double Lambda  , double Omega  , int sign , double min = 0  , double max = 18 ){ 
	double Alpha = TMath::Exp(-min * (Lambda));
	double Beta = TMath::Exp(-max * (Lambda));	
	while( true ){
		double RandomX = -1/Lambda*TMath::Log( Alpha - (Alpha - Beta) * Random.Rndm());
		double RandomY = 1 - 2 * Random.Rndm() ;
		if(RandomY <= TMath::Sin(sign * Omega * RandomX) ){
			obs = RandomX;
			return true;
		}
	}
	return false;
}

bool ExpCosH(TRandom3& Random , double& obs , double Lambda  , double Omega , int sign = 1 , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * (Lambda + Omega));
	double Beta = TMath::Exp(-max * (Lambda + Omega ));	
	double RandomX;
	double RandomY;
	while( true ){
		RandomX = -1/(Lambda + Omega)*TMath::Log(Alpha - (Alpha - Beta) *  Random.Rndm());
		RandomY = sign * Random.Rndm() * TMath::Exp( ( Lambda + Omega ) * RandomX);
		if(RandomY <= TMath::Cos(Omega * RandomX) * TMath::Exp(Lambda * RandomX) ){
			obs = RandomX;
			return true;
		}
	}
	return false;
}

bool ExpSinH(TRandom3& Random , double& obs , double Lambda  , double Omega  , int sign = 1 , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * (Lambda + Omega));
	double Beta = TMath::Exp(-max * (Lambda + Omega));	
	double RandomX;
	double RandomY;
	while( true ){
		RandomX = -1/(Lambda + TMath::Abs(Omega))*TMath::Log( Alpha - (Alpha - Beta) * Random.Rndm());
		RandomY = sign *Random.Rndm() * TMath::Exp( ( Lambda + TMath::Abs(Omega) ) * RandomX);
		if(RandomY <= TMath::SinH(sign * Omega * RandomX) * TMath::Exp(Lambda * RandomX) ){
			obs = RandomX;
			return true;
		}
	}
	return false;
}

int Time_Trafo(TRandom3& Random, double& obs,
				  double a ,double b , double c , double d , 
				  double Gamma_s , double Delta_Gamma_s , double Delta_Mass_s , 
				  double min = 0 ,double max = 18)
{
	double sum_Time = 	TMath::Abs(a) + TMath::Abs(b) + TMath::Abs(c) + TMath::Abs(d);
	double TimeRndm = Random.Rndm();
	if(sum_Time  * TimeRndm < TMath::Abs(a)){
		if (a < 0 )
		{
			ExpCosH( Random , obs , Gamma_s , Delta_Gamma_s/2 , -1 , min , max );			
		} else{
			ExpCosH( Random , obs , Gamma_s , Delta_Gamma_s/2 , +1 , min , max );
		}
		return 1;
	}else if(sum_Time  * TimeRndm < TMath::Abs(a) + TMath::Abs(b)){
		if( b < 0){
			ExpSinH( Random , obs , Gamma_s , Delta_Gamma_s/2 , -1 , min , max );	
		} else{
			ExpSinH( Random , obs , Gamma_s , Delta_Gamma_s/2 , +1 , min , max );
		}
		return 2;
	}else if(sum_Time  * TimeRndm < TMath::Abs(a) + TMath::Abs(b) + TMath::Abs(c)){
		if (c < 0)
		{
			ExpCos( Random , obs , Gamma_s , Delta_Mass_s , -1 , min , max );
		}else{
			ExpCos( Random , obs , Gamma_s , Delta_Mass_s , +1 , min , max );	
		}
		return 3;
	}else{	
		if(c < 0){
			ExpSin( Random , obs , Gamma_s , Delta_Mass_s , -1 , min , max );
		}else{	
			ExpSin( Random , obs , Gamma_s , Delta_Mass_s , +1 , min , max ); 
		}
		return 4;
	}
	return 5;
}

int main()
{
	auto start_time = std::chrono::high_resolution_clock::now();

		unsigned long seed =123; 
		double Daten = 1e6;
		double DatenTime = 1;
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
		//double Phase_Parallel = 2.40;
		//double Phase_Vertical = 2.39;

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
		double Value_4 = 2*Amplitude_Parallel * Amplitude_Vertical;
		double Value_5 = TMath::Sqrt(2) * Amplitude_Parallel * Amplitude_0;
		double Value_6 = TMath::Sqrt(2) * Amplitude_0 * Amplitude_Vertical;

		double Sum = Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6;

	//Creation of the Tree

   		TFile *root_file= new TFile("P2VV.root", "RECREATE");
   		TTree *tree = new TTree("P2VV" , "P2VV");
   		//double Test;
   		//tree->Branch("Test" , &Test , "Test/D" );
   		double Time;
   		tree->Branch("Time" , &Time , "Time/D" );
   		int TimeCut;
   		tree->Branch("TimeCut" , &TimeCut , "TimeCut/I" );
		double theta1;
   		tree->Branch("theta1",&theta1,"theta1/D");
   		double theta2;
   		tree->Branch("theta2",&theta2,"theta2/D");
   		double Phi;
   		tree->Branch("Phi" , &Phi , "Phi/D");
   		int cut;
   		tree->Branch("cut" , &cut , "cut/I");
   		//std::cout << Value_1/4 +Value_2/2 +Value_3/2 << std::endl;

	for (int i = 0 ; i < Daten ; ++i)
	{
		double RndmValue = Random.Rndm();

		if(( Sum * RndmValue < Value_1)){
			for(int j = 0 ; j < DatenTime ; ++j){
				CosSquare(Random , theta1 ,  0  , PI );
				CosSquare(Random , theta2 ,  0  , PI );
				ConstTrafo(Random , Phi   , -PI , PI );
				TimeCut = Time_Trafo(Random , Time ,
								1 , CPV_D , CPV_C , -CPV_S ,
								Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);
				cut = 1;
				tree->Fill();
			}
		}else if( Sum * RndmValue < Value_1 + Value_2){
			for(int j = 0 ; j < DatenTime ; ++j){
				SinSquare(Random , theta1 ,  0  , PI );
				SinSquare(Random , theta2 ,  0  , PI );
				CosSquare(Random , Phi    , -PI , PI );
				TimeCut = Time_Trafo(Random , Time ,
									 1 , CPV_D , CPV_C , -CPV_S ,
									 Gamma_s , Delta_Gamma_s , Delta_Mass_s,  0 , 18);
				cut = 2;
				tree->Fill();
			}
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3){
			for(int j = 0 ; j < DatenTime ; ++j){
				SinSquare(Random , theta1 ,  0  , PI );
				SinSquare(Random , theta2 ,  0  , PI );
				SinSquare(Random , Phi    , -PI , PI );
				TimeCut = Time_Trafo(Random , Time ,
									 1 , -CPV_D , CPV_C , CPV_S ,
									 Gamma_s    , Delta_Gamma_s , Delta_Mass_s,  0 , 18);
				cut = 3;
				tree->Fill();
			}
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4){
			for(int j = 0 ; j < DatenTime ; ++j){
				SinSquare(Random , theta1 ,  0  , PI );
				SinSquare(Random , theta2 ,  0  , PI );
				Sin_Trafo(Random , Phi    , -PI , PI , -1 , 2);
				TimeCut = Time_Trafo(Random , Time ,
								CPV_C * TMath::Sin(Phase_1) , CPV_S * TMath::Sin(Phase_1), TMath::Sin(Phase_1) , CPV_D*TMath::Sin(Phase_1) ,
								Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);
				cut = 4;
				tree->Fill();
			}
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5){
			for(int j = 0 ; j < DatenTime ; ++j){
				Sin_Trafo(Random , theta1 ,  0  , PI , 1 , 2);
				Sin_Trafo(Random , theta2 ,  0  , PI , 1 , 2);
				Cos_Trafo(Random , Phi    , -PI , PI , 1 , 1);
				TimeCut = Time_Trafo(Random , Time ,
						   TMath::Sin(Phase_Parallel) ,CPV_D * TMath::Sin(Phase_Parallel), TMath::Sin(Phase_2) , -CPV_S *TMath::Sin(Phase_Parallel) ,
						   Gamma_s, Delta_Gamma_s , Delta_Mass_s ,  0 , 18);
				cut = 5;
				tree->Fill();
			}
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6){
			for(int j = 0 ; j < DatenTime ; ++j){			
				Sin_Trafo(Random , theta1 , 0   , PI , -1 , 2);
				Sin_Trafo(Random , theta2 , 0   , PI , -1 , 2);
				Sin_Trafo(Random , Phi    , -PI , PI , -1 , 1);
				TimeCut = Time_Trafo(Random , Time ,
						   CPV_C * TMath::Sin(Phase_2) , CPV_S * TMath::Cos(Phase_2) , TMath::Sin(Phase_2) , CPV_D * TMath::Cos(Phase_2) ,
						   Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);
				cut = 6;
				tree->Fill();
			}
		}
		//ExpCos( Random , Test , Gamma_s , Delta_Mass_s); 

		//tree->Fill();

	}
   	root_file->Write("P2VV");

   	auto end_time  = std::chrono::high_resolution_clock::now();     
 	
   	double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();

   	std::cout << "Time: " << duration << " ms"<< std::endl;

   	root_file->Close();
	return 0;
}