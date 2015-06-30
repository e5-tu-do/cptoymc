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
	for( int i = 0 ; i < 1000 ; ++i){
		if((i+1)%100 == 0 ){
			std::cout << std::endl;
			std::cout << "CosSquare.Warning:" <<'\t'<< (i + 1) <<  " times missed" << std::endl;
		}
		double RandomY =  sign * Random.Rndm();
		double RandomX= max - (max - min) * Random.Rndm();
		if(RandomY <= sign * TMath::Cos(Omega * RandomX) * TMath::Cos(Omega * RandomX)){
			obs = RandomX;
			return true;
		}
	}
	std::cout << std::endl;
	std::cerr << "CosSquare: Fatal Error" << std::endl;
	return false;
}
bool SinSquare(TRandom3& Random ,double& obs,  double min , double max , double sign = 1 , double Omega = 1){	//accept and reject of sin^2
	for( int i = 0 ; i < 1000 ; ++i){
		if((i+1)%20 == 0 ){
			std::cout << std::endl;
			std::cout <<"SinSquare.Warning:" << '\t' << (i + 1) <<  " times missed" << std::endl;
		}
		double RandomY = sign * Random.Rndm();
		double RandomX = max - (max - min) * Random.Rndm();
		if(RandomY <= sign * TMath::Sin(Omega * RandomX) * TMath::Sin(Omega * RandomX)){
			obs = RandomX;
			return true;
		}
	}
	std::cout << std::endl;
	std::cerr << "SinSquare: Fatal Error" << std::endl;
	return false;
}
bool Sin_Trafo(TRandom3& Random , double& obs, double min , double max , double sign = 1, double Omega = 1){
	double Alpha = TMath::Cos(Omega * min);
	double Beta = TMath::Cos(Omega * max);
	if (TMath::Abs(Omega *( min - max ) ) < 2*TMath::Pi() && sign > 0 ){
		double Alpha = TMath::Cos(Omega * min);
		double Beta = TMath::Cos(Omega * max);
		obs = 1/Omega * TMath::ACos((Beta - Alpha )*Random.Rndm()+Alpha);
		return true;
	}else /*if (TMath::Abs(Omega*( min - max )) == 2*TMath::Pi() && sign > 0){
		double Alpha = TMath::Cos(Omega * min);
		double Mu = TMath::Cos(Omega * (min + max)/2);
		double Beta = TMath::Cos( Omega * max);
		if(2 * Random.Rndm() < 1 ){
			obs = 1/Omega * TMath::ACos((Mu - Alpha) * Random.Rndm() + Alpha);
			return true;
		}else{
			obs = 1/Omega * TMath::ACos((Beta - Mu )* Random.Rndm() + Mu)+1/Omega*(min + max)/2;
			return true;
		}
	}else*/{
		//if(min == 0 && max == TMath::Pi() && Omega == 1){
	//		obs = 1/Omega*TMath::ACos( ( Beta - Alpha ) *Random.Rndm()+Alpha );
	//}
		for( int i = 0 ; i < 1000 ; ++i){
			if((i+1)%20 == 0 ){
				std::cout << std::endl;
				std::cout <<"Sin.Warning:" << '\t' << (i + 1) <<  " times missed" << std::endl;
			}
			//double RandomY = sign * (1 - 2*Random.Rndm());
			double RandomY = sign * Random.Rndm();
			double RandomX = max - (max - min) * Random.Rndm();
			if( RandomY <= sign * TMath::Abs(TMath::Sin(  RandomX * Omega))){
				obs = RandomX;
				return true;
			}
		}
	}
	std::cout << std::endl;
	std::cerr << "Sin: Fatal Error"<< obs << std::endl;
	return false;
}
bool Cos_Trafo(TRandom3& Random , double& obs, double min , double max , double sign = 1 , double Omega = 1){
	if (TMath::Abs(Omega *( min - max ) ) < 2*TMath::Pi() && sign == 1 ){
		double Alpha = TMath::Sin(Omega * min);
		double Beta = TMath::Sin(Omega * max);
		obs = 1/Omega * TMath::ASin( -(Beta - Alpha )*Random.Rndm() - Alpha);
		return true;
	}else
	{
		for( int i = 0 ; i < 1000 ; ++i){
			if(( i + 1 ) % 20 == 0 ){
				std::cout << std::endl;
				std::cout <<"Cos.Warning:" << '\t' << (i + 1) <<  " times missed" << std::endl;
			}
			double RandomY = sign * Random.Rndm();
			//double RandomY = sign *(1-2*Random.Rndm());
			double RandomX = max - (max - min) * Random.Rndm();
			if(  RandomY <=   sign * TMath::Abs( TMath::Cos( RandomX * Omega ) ) ){
				obs = RandomX;
				return true;
			}
		}
	}
	std::cout << std::endl;
	std::cerr << "Cos: Fatal Error" << std::endl;
	return false;
}
bool ConstTrafo(TRandom3& Random , double& obs , double min , double max){
	obs = (max - min) * Random.Rndm() + min;
	return true;
}
bool ExpCos(TRandom3& Random , double& obs , double Lambda  , double Omega , int sign , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * Lambda );
	double Beta = TMath::Exp(-max * Lambda);
	for( int i = 0 ; i < 1000 ; ++i){
		if(( i + 1 )%600 == 0){
			std::cout << std::endl;
			std::cout << "ExpCos.Warning:" << '\t' << (i + 1) <<  " times missed" << std::endl;
		}
		double RandomX = -1/Lambda*TMath::Log(Alpha - (Alpha - Beta) *  Random.Rndm());
		double RandomY = 1 - 2 * Random.Rndm();
		if(RandomY <= sign * TMath::Cos(Omega * RandomX) ){
			obs = RandomX;
			return true;
		}
	}
	std::cout << std::endl;
	std::cerr << "ExpCos: Fatal Error" << std::endl;
	return false;
}
bool ExpSin(TRandom3& Random , double& obs , double Lambda  , double Omega  , int sign , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * (Lambda));
	double Beta = TMath::Exp(-max * (Lambda));
	for( int i = 0 ; i < 1000 ; ++i){
		if((i + 1)%600 == 0 ){
			std::cout << std::endl;
			std::cout << "ExpSin.Warning:" << '\t' << (i + 1) <<  " times missed" << std::endl;
		}
		double RandomX = -1/Lambda*TMath::Log( Alpha - (Alpha - Beta) * Random.Rndm());
		double RandomY =  Random.Rndm() ;
		if(RandomY <= TMath::Sin(sign * Omega * RandomX) ){
			obs = RandomX;
			return true;
		}
	}
	std::cout << std::endl;
	std::cerr << "ExpSin: Fatal Error" << std::endl;
	return false;
}
bool ExpCosH(TRandom3& Random , double& obs , double Lambda  , double Omega , int sign = 1 , double min = 0  , double max = 18 ){
	double Alpha = TMath::Exp(-min * (Lambda + Omega));
	double Beta = TMath::Exp(-max * (Lambda + Omega ));
	double RandomX;
	double RandomY;
	for( int i = 0 ; i < 1000 ; ++i){
		if((i + 1 )%600 == 0 ){
			std::cout << std::endl;
			std::cout << "ExpCosH.Warning:" << '\t' << (i + 1) <<  " times missed" << std::endl;
		}
		RandomX = -1/(Lambda + Omega)*TMath::Log(Alpha - (Alpha - Beta) *  Random.Rndm());
		RandomY = sign * Random.Rndm() * TMath::Exp( ( Lambda + Omega ) * RandomX);
		if(RandomY <= sign * TMath::Cos(Omega * RandomX) * TMath::Exp(Lambda * RandomX) ){
			obs = RandomX;
			return true;
		}
	}
	std::cout << std::endl;
	std::cerr << "ExpCosH: Fatal Error" << std::endl;
	return false;
}
bool ExpSinH(TRandom3& Random , double& obs , double Lambda  , double Omega  , int sign = 1 , double min = 0  , double max = 18 ){
	//double Alpha = TMath::Exp(-min * (Lambda + Omega));
	//double Beta = TMath::Exp(-max * (Lambda + Omega));	
	double Alpha = TMath::Exp(-min * (Lambda ));
	double Beta = TMath::Exp(-max * (Lambda ));	
	double RandomX;
	double RandomY;
	bool Warning =false;
	for( int i = 0 ; i < 2000 ; ++i){
		if((i + 1 )%1600 == 0 && Warning == false){
			Warning = true;
		}
		RandomX = -1/(Lambda)*TMath::Log( Alpha - (Alpha - Beta) * Random.Rndm());
		RandomY = sign *Random.Rndm() * TMath::Exp( ( Lambda ) * RandomX);
		if(RandomY <=  sign*TMath::SinH( Omega * RandomX) * TMath::Exp(-Lambda * RandomX) ){
			if (Warning == true)
			{
				std::cout << std::endl;
				std::cout << "ExpSinH.Warning:" << '\t'<< (i + 1) <<  " times missed" << std::endl;
			}
			obs = RandomX;
			return true;
		}
	}
	std::cout << std::endl;
	std::cerr << "ExpSinH: Fatal Error" << std::endl;
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
		unsigned int Daten = 1e6;
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
		double Amplitude_S = 0;
		double Amplitude_SS = 0;

		double Phase_0 = 0 ;
		double Phase_Parallel = 2.71;
		double Phase_Vertical = 2.39;
		double Phase_S = 0;
		double Phase_SS = 0;

		double CP_Lambda = 1;
		double CP_Phi_sss = 0; // = 0.17

	// CP-Violation Parameter
		double CPV_C = (1- pow(CP_Lambda,2))/(1+pow(CP_Lambda,2));
		double CPV_S = -2 * CP_Lambda * TMath::Sin(CP_Phi_sss)/(1+pow(CP_Lambda,2));
		double CPV_D = -2 * CP_Lambda * TMath::Cos(CP_Phi_sss)/(1+pow(CP_Lambda,2));

	//Defination 
		double Phase_1 = Phase_Vertical - Phase_Parallel;
		double Phase_2 = Phase_Vertical - Phase_0;
		double Phase_21 = Phase_2 - Phase_1 ;
		double Time_Begin = 0;
		double Time_End = 20;
		double theta_Begin = 0;
		double theta_End = TMath::Pi();
		double Phi_Begin = -TMath::Pi();
		double Phi_End = TMath::Pi();

		double Value_1 = 1/(Amplitude_0  		* Amplitude_0 			*4)										;
		double Value_2 = 1/(Amplitude_Parallel * Amplitude_Parallel 	*2)										;
		double Value_3 = 1/(Amplitude_Vertical * Amplitude_Vertical 	*2)										;
		double Value_4 = 1/(Amplitude_Parallel * Amplitude_Vertical 	*2)										;
		double Value_5 = 1/(Amplitude_Parallel * Amplitude_0			*TMath::Sqrt(2))						;
		double Value_6 = 1/(Amplitude_0 		* Amplitude_Vertical 	*TMath::Sqrt(2))						;

		double Value_7 = Amplitude_SS		* Amplitude_SS			*4/9 														;
		double Value_8 = Amplitude_S 		* Amplitude_SS 			*4/3 									  					;
		double Value_9 = Amplitude_SS 		* Amplitude_SS			*8/(3*TMath::Sqrt(3))										;
		double Value_10= Amplitude_0 		* Amplitude_SS			*8/3 					* TMath::Cos(Phase_SS)				;
		double Value_11= Amplitude_Parallel * Amplitude_SS			*4* TMath::Sqrt(2)/3 	* TMath::Cos(Phase_21 - Phase_SS)	;
		double Value_12= Amplitude_Vertical * Amplitude_SS			*4* TMath::Sqrt(2)/3 										;
		double Value_13= Amplitude_0 		* Amplitude_S 			*8/TMath::Sqrt(3) 											;
		double Value_14= Amplitude_Parallel * Amplitude_S 			*4*TMath::Sqrt(2/3)											;
		double Value_15= Amplitude_Vertical * Amplitude_S 			*4*TMath::Sqrt(2/3)		* TMath::Sin(Phase_2 - Phase_S)		;//*/

		/*double Value_1 = Amplitude_0  		* Amplitude_0 			;
		double Value_2 = Amplitude_Parallel * Amplitude_Parallel 	;
		double Value_3 = Amplitude_Vertical * Amplitude_Vertical 	;
		double Value_4 = Amplitude_Parallel * Amplitude_Vertical 	;
		double Value_5 = Amplitude_Parallel * Amplitude_0			* TMath::Cos(Phase_21);
		double Value_6 = Amplitude_0 		* Amplitude_Vertical 	;

		double Value_7 = Amplitude_SS		* Amplitude_SS			;
		double Value_8 = Amplitude_S 		* Amplitude_SS 			;
		double Value_9 = Amplitude_SS 		* Amplitude_SS			;
		double Value_10= Amplitude_0 		* Amplitude_SS			;
		double Value_11= Amplitude_Parallel * Amplitude_SS			;
		double Value_12= Amplitude_Vertical * Amplitude_SS			;
		double Value_13= Amplitude_0 		* Amplitude_S 			;
		double Value_14= Amplitude_Parallel * Amplitude_S 			;
		double Value_15= Amplitude_Vertical * Amplitude_S 			;//*/

		double Sum = Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6 + Value_7 + Value_8 + Value_9 + Value_10 + Value_11 + Value_12 + Value_13 + Value_14 + Value_15;

	//Creation of the Tree

   		TFile *root_file= new TFile("P2VV.root", "RECREATE");
   		TTree *tree = new TTree("P2VV" , "P2VV");

   		double Test_Cos;
   		tree->Branch("Test_Cos" , &Test_Cos , "Test_Cos/D");
   		double Test_Sin;
   		tree->Branch("Test_Sin" , &Test_Sin , "Test_Sin/D");
   		double Test_CosSquare;
   		tree->Branch("Test_CosSquare" , &Test_CosSquare , "Test_CosSquare/D");
   		double Test_SinSquare;
   		tree->Branch("Test_SinSquare" , &Test_SinSquare , "Test_SinSquare/D");
   		double Test_ExpSin;
   		tree->Branch("Test_ExpSin" , &Test_ExpSin , "Test_ExpSin/D");
   		double Test_ExpCos;
   		tree->Branch("Test_ExpCos" , &Test_ExpCos , "Test_ExpCos/D");

   		double Time;
   		tree->Branch("Time" , &Time , "Time/D" );
		double theta1;
   		tree->Branch("theta1",&theta1,"theta1/D");
   		double theta2;
   		tree->Branch("theta2",&theta2,"theta2/D");
   		double Phi;
   		tree->Branch("Phi" , &Phi , "Phi/D");
   		int TimeCut;
   		tree->Branch("TimeCut" , &TimeCut , "TimeCut/I");
   		int cut;
   		tree->Branch("cut" , &cut , "cut/I");
   		double cos_theta1;
   		tree->Branch("cos_theta1" , &cos_theta1 , "cos_theta1/D");
   		double cos_theta2;
   		tree->Branch("cos_theta2" , &cos_theta2 , "cos_theta2/D");
   		double cos_Phi;
   		tree->Branch("cos_Phi" , &cos_Phi , "cos_Phi/D");
   		double sin_Phi;
   		tree->Branch("sin_Phi" , &sin_Phi , "sin_Phi/D");
   		double B_s0_f1;
   		tree->Branch("B_s0_f1" , &B_s0_f1 , "B_s0_f1/D");
   		double B_s0_f2;
   		tree->Branch("B_s0_f2" , &B_s0_f2 , "B_s0_f2/D");
   		double B_s0_f3;
   		tree->Branch("B_s0_f3" , &B_s0_f3 , "B_s0_f3/D");
   		double B_s0_f4;
   		tree->Branch("B_s0_f4" , &B_s0_f4 , "B_s0_f4/D");
   		double B_s0_f5;
   		tree->Branch("B_s0_f5" , &B_s0_f5 , "B_s0_f5/D");
   		double B_s0_f6;
   		tree->Branch("B_s0_f6" , &B_s0_f6 , "B_s0_f6/D");

   		/*double B_s0_f7;
   		tree->Branch("B_s0_f7" , &B_s0_f7 , "B_s0_f7/D");
   		double B_s0_f8;
   		tree->Branch("B_s0_f8" , &B_s0_f8 , "B_s0_f8/D");
   		double B_s0_f9;
   		tree->Branch("B_s0_f9" , &B_s0_f9 , "B_s0_f9/D");
   		double B_s0_f10;
   		tree->Branch("B_s0_f10" , &B_s0_f10 , "B_s0_f10/D");
   		double B_s0_f11;
   		tree->Branch("B_s0_f11" , &B_s0_f11 , "B_s0_f11/D");
   		double B_s0_f12;
   		tree->Branch("B_s0_f12" , &B_s0_f12 , "B_s0_f12/D");
   		double B_s0_f13;
   		tree->Branch("B_s0_f13" , &B_s0_f13 , "B_s0_f13/D");
   		double B_s0_f14;
   		tree->Branch("B_s0_f14" , &B_s0_f14 , "B_s0_f14/D");
   		double B_s0_f15;
   		tree->Branch("B_s0_f15" , &B_s0_f15 , "B_s0_f15/D");*/

	for (unsigned int i = 0 ; i < Daten ; ++i)
	{
		if ((i+1)%10000 == 0)
		{
			std::cout << (double)(i+1)*100 / (Daten)<< '%' << '\t' << std::flush;
			if ((i+1)%100000 == 0)
			{
				//std::cout << '\r' << std::flush;
				std::cout << std::endl;
			}
		}
		double RndmValue = Random.Rndm();

		if(( Sum * RndmValue < Value_1)){
			CosSquare(Random , theta1 ,  theta_Begin  , theta_End );
			CosSquare(Random , theta2 ,  theta_Begin  , theta_End );
			ConstTrafo(Random , Phi   , Phi_Begin , Phi_End );
			TimeCut = Time_Trafo(Random , Time ,
								 1 , CPV_D , CPV_C , -CPV_S ,
								 Gamma_s, Delta_Gamma_s , Delta_Mass_s,  Time_Begin , Time_End );
			cut = 1;
		}else if( Sum * RndmValue < Value_1 + Value_2){
			SinSquare(Random , theta1 , theta_Begin  , theta_End );
			SinSquare(Random , theta2 , theta_Begin  , theta_End );
			CosSquare(Random , Phi    , Phi_Begin , Phi_End );
			TimeCut = Time_Trafo(Random , Time ,
								 1 , CPV_D , CPV_C , -CPV_S ,
								 Gamma_s , Delta_Gamma_s , Delta_Mass_s,  Time_Begin , Time_End );
			cut = 2;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3){
			SinSquare(Random , theta1 , theta_Begin  , theta_End );
			SinSquare(Random , theta2 , theta_Begin  , theta_End );
			SinSquare(Random , Phi    , Phi_Begin , Phi_End );
			TimeCut = Time_Trafo(Random , Time ,
								 1 , -CPV_D , CPV_C , CPV_S ,
								 Gamma_s    , Delta_Gamma_s , Delta_Mass_s,  Time_Begin , Time_End );
			cut = 3;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4){
			SinSquare(Random , theta1 ,  theta_Begin  , theta_End );
			SinSquare(Random , theta2 ,  theta_Begin  , theta_End );
			Sin_Trafo(Random , Phi    , Phi_Begin , Phi_End , 1 , -2);
			TimeCut = Time_Trafo(Random , Time ,
								 CPV_C * TMath::Sin(Phase_1) , CPV_S * TMath::Sin(Phase_1), TMath::Sin(Phase_1) , CPV_D*TMath::Sin(Phase_1) ,
								 Gamma_s, Delta_Gamma_s , Delta_Mass_s,  Time_Begin , Time_End );
			cut = 4;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5){
			Sin_Trafo(Random , theta1 , theta_Begin  , theta_End , 1 , 1);
			Sin_Trafo(Random , theta2 , theta_Begin  , theta_End , 1 , 1);
			Cos_Trafo(Random , Phi    , Phi_Begin , Phi_End , 1 , 1);
			TimeCut = Time_Trafo(Random , Time ,
							 	  1,CPV_D , 1 , -CPV_S  ,
						   		  Gamma_s, Delta_Gamma_s , Delta_Mass_s ,  Time_Begin , Time_End );
			cut = 5;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6){
			Sin_Trafo(Random , theta1 , theta_Begin  , theta_End ,  1 , 1);
			Sin_Trafo(Random , theta2 , theta_Begin  , theta_End ,  1 , 1);
			Sin_Trafo(Random , Phi    , Phi_Begin , Phi_End ,  1 , -1);
			TimeCut = Time_Trafo(Random , Time ,
					   CPV_C * TMath::Sin(Phase_2) , CPV_S * TMath::Cos(Phase_2) , TMath::Sin(Phase_2) , CPV_D * TMath::Cos(Phase_2) ,
					   Gamma_s, Delta_Gamma_s , Delta_Mass_s,  Time_Begin , Time_End );
			cut = 6;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6 + Value_7){
			ConstTrafo(Random , theta1 , theta_Begin  , theta_End );
			ConstTrafo(Random , theta2 , theta_Begin  , theta_End );
			ConstTrafo(Random , Phi ,  Phi_Begin , Phi_End );
			TimeCut = Time_Trafo(Random , Time ,
					   1 , CPV_D , CPV_C , -CPV_S ,
					   Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 16);
			cut = 7;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6 + Value_7 + Value_8){
			double RandomValue = 4*Random.Rndm();
			if(RandomValue < 1){
				CosSquare(Random , theta1 , theta_Begin  , theta_End);
				ConstTrafo(Random , theta2 , theta_Begin  , theta_End);
			}else if(RandomValue < 2){
				CosSquare(Random , theta2 , theta_Begin  , theta_End );
				ConstTrafo(Random , theta1 , theta_Begin  , theta_End );
			} else {
				Cos_Trafo(Random , theta1 , theta_Begin  , theta_End );
				Cos_Trafo(Random , theta2 , theta_Begin  , theta_End );
			}
			ConstTrafo(Random , Phi ,  Phi_Begin , Phi_End );
			TimeCut = Time_Trafo(Random , Time ,
					   1 , -CPV_D , CPV_C , CPV_S ,
					   Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 16);
			cut = 8;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6 + Value_7 + Value_8 + Value_9){
			if (Random.Rndm() < 0.5)
			{
				Cos_Trafo(Random , theta1 , theta_Begin  , theta_End );
				ConstTrafo(Random , theta2 , theta_Begin  , theta_End );
			}else{
				Cos_Trafo(Random , theta2 , theta_Begin  , theta_End );
				ConstTrafo(Random , theta1 , theta_Begin  , theta_End );
			}
			ConstTrafo(Random , Phi ,  Phi_Begin , Phi_End );
			TimeCut = Time_Trafo(Random , Time ,
					   CPV_C*TMath::Cos(Phase_S - Phase_SS) , CPV_S*TMath::Sin(Phase_S - Phase_SS), TMath::Cos(Phase_S - Phase_SS) , TMath::Sin(Phase_S - Phase_SS) ,
					   Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 16);
			cut = 9;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6 + Value_7 + Value_8 + Value_9 + Value_10){
			Cos_Trafo(Random , theta1 , theta_Begin  , theta_End );
			Cos_Trafo(Random , theta2 , theta_Begin  , theta_End );
			ConstTrafo(Random , Phi ,  Phi_Begin , Phi_End );
			TimeCut = Time_Trafo(Random , Time ,
					   1 , CPV_D , CPV_C , -CPV_S ,
					   Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 16);
			cut = 10;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6 + Value_7 + Value_8 + Value_9 + Value_10 + Value_11){
			Sin_Trafo(Random , theta1 , theta_Begin  , theta_End );
			Sin_Trafo(Random , theta2 , theta_Begin  , theta_End );
			Cos_Trafo(Random , Phi ,  Phi_Begin , Phi_End );
			TimeCut = Time_Trafo(Random , Time ,
					   1 , CPV_D , CPV_C , -CPV_S ,
					   Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 16);
			cut = 11;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6 + Value_7 + Value_8 + Value_9 + Value_10 + Value_11 + Value_12){
			Sin_Trafo(Random , theta1 , theta_Begin  , theta_End );
			Sin_Trafo(Random , theta2 , theta_Begin  , theta_End );
			Cos_Trafo(Random , Phi ,  Phi_Begin , Phi_End );
			TimeCut = Time_Trafo(Random , Time ,
					   1 , CPV_D , CPV_C , -CPV_S ,
					   Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 16);
			cut = 12;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6 + Value_7 + Value_8 + Value_9 + Value_10 + Value_11 + Value_12 + Value_13){

		}
		Cos_Trafo(Random , Test_Cos , 0, 2*PI , 1 , 1);
		Sin_Trafo(Random , Test_Sin , 0  , 2*PI , 1 , 1);
		SinSquare(Random , Test_SinSquare , -PI , PI , 1 , 1);
		CosSquare(Random , Test_CosSquare , -PI , PI , 1 , 1);
		ExpSin(Random , Test_ExpSin , 1 , 1 , 1 , 0 , 5 );
		ExpCos(Random , Test_ExpCos , 1 , 1 , 1 , 0 , 5 );

		cos_theta1 = TMath::Cos(theta1);
		cos_theta2 = TMath::Cos(theta2);
		cos_Phi = TMath::Cos(Phi);
		sin_Phi = TMath::Sin(Phi);

	   	B_s0_f1 = 4*TMath::Cos(theta1)*TMath::Cos(theta1)*TMath::Cos(theta2)*TMath::Cos(theta2);
	   	B_s0_f2 = 2*TMath::Sin(theta1)*TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Sin(theta2)*TMath::Cos(Phi)*TMath::Cos(Phi);
	   	B_s0_f3 = 2*TMath::Sin(theta1)*TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Sin(theta2)*TMath::Sin(Phi)*TMath::Sin(Phi);
	   	B_s0_f4 = -2*TMath::Sin(theta1)*TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Sin(theta2)*TMath::Sin(2*Phi);	 
	   	B_s0_f5 = TMath::Sqrt(2)*TMath::Sin(2*theta1)*TMath::Sin(2*theta2)*TMath::Cos(Phi);
	   	B_s0_f6 = -TMath::Sqrt(2)*TMath::Sin(2*theta1)*TMath::Sin(2*theta2)*TMath::Sin(Phi);

	   	/*B_s0_f7 = 4/9;
	   	B_s0_f8 = 4/3 * (TMath::Cos(theta1) + TMath::Cos(theta2))*(TMath::Cos(theta1) + TMath::Cos(theta2));
	   	B_s0_f9 = 8/(3*TMath::Sqrt(3))*(TMath::Cos(theta1) + TMath::Cos(theta2));
	   	B_s0_f10= 8/3*TMath::Cos(theta1)*TMath::Cos(theta2);
	   	B_s0_f11= 4*TMath::Sqrt(2)/3*TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Cos(Phi);
	   	B_s0_f12=-4*TMath::Sqrt(2)/3*TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Sin(Phi);
	   	B_s0_f13=8/TMath::Sqrt(3)*TMath::Cos(theta1)*TMath::Cos(theta2)*(TMath::Cos(theta1)+TMath::Cos(theta2));
	   	B_s0_f14=4*TMath::Sqrt(2/3)*TMath::Sin(theta1)*TMath::Sin(theta2)*(TMath::Cos(theta1)+TMath::Cos(theta2))*TMath::Cos(Phi);
	   	B_s0_f15=-4*TMath::Sqrt(2/3)*TMath::Sin(theta1)*TMath::Sin(theta2)*(TMath::Cos(theta1)+TMath::Cos(theta2))*TMath::Sin(Phi);*/

		tree->Fill();
	}
   	root_file->Write("P2VV");

   	auto end_time  = std::chrono::high_resolution_clock::now();
 	
   	double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();

   	std::cout << "Time: " << duration << " ms"<< std::endl;

   	root_file->Close();
	return 0;
}