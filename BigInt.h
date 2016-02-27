#pragma once
#include"DefinesAndConstatns.h"
using namespace std;
class BigInt
{
public:
	BigInt(string Number0); 
	BigInt();
	~BigInt();
	
	//BigInt(vll v);
	
	BigInt(vs v);
	
	BigInt operator +(const BigInt&);  
	BigInt operator -(const BigInt&);	
	BigInt operator *(const BigInt&);	
	BigInt operator /(const BigInt&);	
	BigInt operator %(const BigInt&);	
	BigInt operator -();

	BigInt operator+ (ll);
	BigInt operator- (ll);
	BigInt operator* (ll);
	BigInt operator/ (ll);
	ll operator% (ll);
	void operator+= (ll);
	void operator-= (ll);
	void operator*= (ll);
	void operator/= (ll);
	void operator%= (ll);


	void operator ++();
	void operator --();
	void operator =(const BigInt&);
	void operator =(ll value);
	void operator +=(const BigInt&);	// return number in basis of first argument
	void operator -=(const BigInt&);	// return number in basis of first argument
	void operator *=(const BigInt&);	// return number in basis of first argument
	void operator /=(const BigInt&);	// return number in basis of first argument
	void operator %=(const BigInt&);
	
	bool operator >(const BigInt&) const ;
	bool operator >=(const BigInt&) const;
	bool operator <(const BigInt&) const;
	bool operator <=(const BigInt&) const;
	bool operator ==(const BigInt&)const ;
	bool operator !=(const BigInt&)const ;


	bool operator >( ll) const;
	bool operator >=( ll) const;
	bool operator <( ll) const;
	bool operator <=( ll) const;
	bool operator ==( ll)const;
	bool operator !=( ll)const;

	BigInt operator << (int v);
	BigInt operator >>(int v);
	
	void input(string s);
	string	GetNumber();
	void  PrintBin();
	
	BigInt abs() const;
	void trim();
	int FirstOne() const;
	bool isZero() const;
	
	BigInt Cut(ll k);

	static vs KaratsubaMiltiplication(const vs &a, const vs &b);
	

	void fft(vc & a, bool inverse);
	void fft_multiply(const vector<ll> & a, const vector<ll> & b, vector<ll> & res);
	BigInt DivKuka(BigInt a, BigInt b);
	BigInt Shtrassen(BigInt a, BigInt b);
	BigInt KaratsubaMultiply(const BigInt&);
	BigInt Dixon(BigInt N);
	BigInt ToomaKuka(BigInt a, BigInt b);
	BigInt Shouhage(const BigInt& a, const BigInt& b);
	BigInt Toom(BigInt a, BigInt b);
	vector<BigInt> PolardP_1(BigInt N);
	vector<BigInt> PollardP(BigInt N);
	BigInt CookInverse(BigInt B);

	BigInt PowModN(BigInt a, BigInt m, BigInt N);
	BigInt PowN(BigInt a, BigInt N);
	BigInt gcd(BigInt a, BigInt b);
	BigInt sqrt(BigInt);
	bool issquare(BigInt);
	double GetTime();
//private:
	
	int CookInverse(BigInt v, BigInt& result, int c);
	static const long long Standart_basis;
	static const int MaxLengthOfNumber=1000000;
	static const int DecBasis = 10;
	float AfterPoint;
	int sign;
	int DischargeLenght;
	ll basis;
	static const int LenghOfFB = 303;
	static const int FactorBase_2000[LenghOfFB];
	int PointPos;
	vs Number;
	vll DecNumber;
	double TimeOfLast;

	BigInt PollardP_get_factor(BigInt N);
	BigInt PolardP_1_get_factor(BigInt K,BigInt a);
	ll gcdex(ll, ll, ll&, ll&);
	vll DecMult(vll a, vll b);

	vs SetInBinaryBasis();

	friend pair<BigInt, BigInt> DivMode(const BigInt &, const BigInt &);
	ll ConvertToLL(vs);
	vs ConvertToBin(ll);
	void MakeDecNumber();
	vll DecPlus(vll a, vll b);
	vll DecMinus(vll a, vll b);
	vll DivToLL(vll,ll);
	int ModLL(vll, ll);
	BigInt FastModule(ll v) const ;
	ll gcd(ll a, ll b, ll & x, ll & y);
	ll GetRes(ll a, ll mod);

	BigInt HelpMult(ll , ll , ll , BigInt );
	
	

	bool is_Smooth(BigInt A, vector<int> FBase,vector<vector<int> > &a, vector<vector<int> > &eps);
	


	bool isPrime(BigInt N);
	
	BigInt Dixon_get_factor(BigInt N);
	void MakeFactorBase(vector<int> &B, int Base);
	vector<int> Eratosfen(int B);
	BigInt LCM(int K);
	vector<int> Polard_Help_Eratosfen(int B1, int B2);
	int gauss(vector < vector<int> > a, vector<int> & ans);
	BigInt GerRandomBettwen(BigInt L, BigInt R);
	
};






