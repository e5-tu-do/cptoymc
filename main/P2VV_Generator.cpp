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
	if (TMath::Abs(Omega *( min - max ) ) < 2*TMath::Pi() && sign != 1 ){
		double Alpha = TMath::Cos(Omega * min);
		double Beta = TMath::Cos(Omega * max);
		obs = 1/Omega * TMath::ACos((Beta - Alpha )*Random.Rndm()-Alpha);
		return true;
	}else {
		while( true ){
			double RandomY = sign * (1 - 2*Random.Rndm());
			//double RandomY = sign * Random.Rndm();
			double RandomX = max - (max - min) * Random.Rndm();
			if( RandomY <= sign * TMath::Sin(  RandomX * Omega)){
				obs = RandomX;
				return true;
			}
		}
	}
	return false;
}
bool Cos_Trafo(TRandom3& Random , double& obs, double min , double max , double sign = 1 , double Omega = 1){
	if (TMath::Abs(Omega *( min - max ) ) < 2*TMath::Pi() && sign != 1 ){
		double Alpha = TMath::Sin(Omega * min);
		double Beta = TMath::Sin(Omega * max);
		obs = 1/Omega * TMath::ASin((Beta - Alpha )*Random.Rndm()-Alpha);
		return true;
	}else {
		while( true ){
			double RandomY = sign * (1 - 2 * Random.Rndm());
			//double RandomY = sign *Random.Rndm();
			double RandomX = max - (max - min) * Random.Rndm();
			if(  RandomY <=   sign *TMath::Cos(RandomX * Omega)){
				obs = RandomX;
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
		if(RandomY <= sign * TMath::Cos(Omega * RandomX) * TMath::Exp(Lambda * RandomX) ){
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
		RandomX = -1/(Lambda + Omega)*TMath::Log( Alpha - (Alpha - Beta) * Random.Rndm());
		RandomY = sign *Random.Rndm() * TMath::Exp( ( Lambda + Omega ) * RandomX);
		if(RandomY <=  sign*TMath::SinH( Omega * RandomX) * TMath::Exp(-Lambda * RandomX) ){
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
		//double Amplitude_S = 0;
		//double Amplitude_SS = 0;
		/*double Amplitude_0 = 0.701;
		double Amplitude_Parallel = 0.506;
		double Amplitude_Vertical = 0.502;*/

		double Phase_0 = 0 ;
		double Phase_Parallel = 2.71;
		double Phase_Vertical = 2.39;
		//double Phase_S = 0;
		//double Phase_SS = 0;
		/*double Phase_Parallel = 2.40;
		double Phase_Vertical = 2.39;*/

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

		/*double Value_1 = 4  			* Amplitude_0 		 * Amplitude_0 			*PI*PI*PI/2;
		double Value_2 = 2   			* Amplitude_Parallel * Amplitude_Parallel 	*PI*PI*PI/4;
		double Value_3 = 2   			* Amplitude_Vertical * Amplitude_Vertical 	*PI*PI*PI/4;
		double Value_4 = 2  			* Amplitude_Parallel * Amplitude_Vertical 	*PI*PI;
		double Value_5 = TMath::Sqrt(2) * Amplitude_Parallel * Amplitude_0			*16;
		double Value_6 = TMath::Sqrt(2) * Amplitude_0 		 * Amplitude_Vertical 	*16;*/

		/*double Value_1 = 4 * 			  Amplitude_0 		 * Amplitude_0 			;
		double Value_2 = 2 * 			  Amplitude_Parallel * Amplitude_Parallel 	;
		double Value_3 = 2 * 			  Amplitude_Vertical * Amplitude_Vertical 	;
		double Value_4 = 2 *			  Amplitude_Parallel * Amplitude_Vertical 	;
		double Value_5 = TMath::Sqrt(2) * Amplitude_Parallel * Amplitude_0			;
		double Value_6 = TMath::Sqrt(2) * Amplitude_0 		 * Amplitude_Vertical 	;*/

		double Value_1 = 	 			  Amplitude_0 		 * Amplitude_0 			;
		double Value_2 = 	 			  Amplitude_Parallel * Amplitude_Parallel 	;
		double Value_3 = 	 			  Amplitude_Vertical * Amplitude_Vertical 	;
		double Value_4 = 				  Amplitude_Parallel * Amplitude_Vertical 	;
		double Value_5 = 				  Amplitude_Parallel * Amplitude_0			;
		double Value_6 = 				  Amplitude_0 		 * Amplitude_Vertical 	;

		/*double Value_7 = 4/9 * Amplitude_SS*Amplitude_SS;
		double Value_8 = 4/3 * Amplitude_S *Amplitude_S;
		double Value_9 = 8/(3*TMath::Sqrt(3))* Amplitude_S*Amplitude_SS;
		double Value_10= 8/3 * Amplitude_0 * Amplitude_SS;
		double Value_11= 4* TMath::Sqrt(2)/3 * Amplitude_Parallel * Amplitude_SS;
		double Value_12= 4* TMath::Sqrt(2)/3 * Amplitude_Vertical * Amplitude_SS;*/


		double Sum = Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6;

	//Creation of the Tree

   		TFile *root_file= new TFile("P2VV.root", "RECREATE");
   		TTree *tree = new TTree("P2VV" , "P2VV");
   		double Test1;
   		tree->Branch("Test1" , &Test1 , "Test1/D");
   		double Test2;
   		tree->Branch("Test2" , &Test2 , "Test2/D");
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

	for (int i = 0 ; i < Daten ; ++i)
	{
		double RndmValue = Random.Rndm();

		if(( Sum * RndmValue < Value_1)){
			CosSquare(Random , theta1 ,  0  , PI );
			CosSquare(Random , theta2 ,  0  , PI );
			ConstTrafo(Random , Phi   , -PI , PI );
			TimeCut = Time_Trafo(Random , Time ,
								 1 , CPV_D , CPV_C , -CPV_S ,
								 Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);
			cut = 1;
		}else if( Sum * RndmValue < Value_1 + Value_2){
			SinSquare(Random , theta1 ,  0  , PI );
			SinSquare(Random , theta2 ,  0  , PI );
			CosSquare(Random , Phi    , -PI , PI );
			TimeCut = Time_Trafo(Random , Time ,
								 1 , CPV_D , CPV_C , -CPV_S ,
								 Gamma_s , Delta_Gamma_s , Delta_Mass_s,  0 , 18);
			cut = 2;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3){
			SinSquare(Random , theta1 ,  0  , PI );
			SinSquare(Random , theta2 ,  0  , PI );
			SinSquare(Random , Phi    , -PI , PI );
			TimeCut = Time_Trafo(Random , Time ,
								 1 , -CPV_D , CPV_C , CPV_S ,
								 Gamma_s    , Delta_Gamma_s , Delta_Mass_s,  0 , 18);
			cut = 3;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4){
			SinSquare(Random , theta1 ,  0  , PI );
			SinSquare(Random , theta2 ,  0  , PI );
			Sin_Trafo(Random , Phi    , -PI , PI , -1 , 2);
			//Sin_Trafo(Random , Phi    , -PI , PI , 1 , 2);
			TimeCut = Time_Trafo(Random , Time ,
								 CPV_C * TMath::Sin(Phase_1) , CPV_S * TMath::Sin(Phase_1), TMath::Sin(Phase_1) , CPV_D*TMath::Sin(Phase_1) ,
								 Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);
			cut = 4;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5){
			/*Sin_Trafo(Random , theta1 ,  0  , PI , 1 , 2);
			Sin_Trafo(Random , theta2 ,  0  , PI , 1 , 2);*/
			Sin_Trafo(Random , theta1 ,  0  , PI , -1 , 2);
			Sin_Trafo(Random , theta2 ,  0  , PI , -1 , 2);
			Cos_Trafo(Random , Phi    , -PI , PI , 1 , 1);
			TimeCut = Time_Trafo(Random , Time ,
							 	 1,CPV_D , 1 , -CPV_S  ,
						   		 Gamma_s, Delta_Gamma_s , Delta_Mass_s ,  0 , 18);
			cut = 5;
		}else if(Sum * RndmValue < Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6){
			/*Sin_Trafo(Random , theta1 , 0   , PI ,  -1 , 2);
			Sin_Trafo(Random , theta2 , 0   , PI ,  -1 , 2);
			Sin_Trafo(Random , Phi    , -PI , PI ,  -1 , 1);*/
			Sin_Trafo(Random , theta1 , 0   , PI ,  1 , 2);
			Sin_Trafo(Random , theta2 , 0   , PI ,  1 , 2);
			Sin_Trafo(Random , Phi    , -PI , PI ,  1 , 1);
			TimeCut = Time_Trafo(Random , Time ,
					   CPV_C * TMath::Sin(Phase_2) , CPV_S * TMath::Cos(Phase_2) , TMath::Sin(Phase_2) , CPV_D * TMath::Cos(Phase_2) ,
					   Gamma_s, Delta_Gamma_s , Delta_Mass_s,  0 , 18);
			cut = 6;
		}
		SinSquare(Random , Test1 , 0 , PI , 2 , 2);
		SinSquare(Random , Test2 , 0 , PI , 1 , 2);

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
   	std::cout << "Value_1/sum:" << Value_1 / Sum << std::endl;
   	std::cout << "Value_2/sum:" << Value_2 / Sum << std::endl;
   	std::cout << "Value_3/sum:" << Value_3 / Sum << std::endl;
   	std::cout << "Value_4/sum:" << Value_4 / Sum << std::endl;
   	std::cout << "Value_5/sum:" << Value_5 / Sum << std::endl;
   	std::cout << "Value_6/sum:" << Value_6 / Sum << std::endl;
   	std::cout << (Value_1 + Value_2 + Value_3 + Value_4 + Value_5 + Value_6)/Sum << std::endl;

   	auto end_time  = std::chrono::high_resolution_clock::now();
 	
   	double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();

   	std::cout << "Time: " << duration << " ms"<< std::endl;

   	root_file->Close();
	return 0;
}