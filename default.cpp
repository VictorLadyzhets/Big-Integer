#include"BigInt.h"


using namespace std;


int main()
{
	BigInt a("13");
	BigInt b("7");
	BigInt c("48");
	BigInt d("8");
	
	/*for (int i = 0; i < 27; i++)
	{
		a = a.KaratsubaMultiply(a);
	}*/
	a = a.PowModN(a, BigInt("28"), BigInt("47"));
	
	cout << a.GetNumber();


	system("pause");
	return 0;

}



	   