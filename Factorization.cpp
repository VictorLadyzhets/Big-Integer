#include "BigInt.h"
/*
BigInt BigInt::Dixon_get_factor(BigInt N)
{
	
	int B = exp(0.5*std::sqrt(N.Number.size()*log(N.Number.size())));
	vector<int> fBase;
	MakeFactorBase(fBase, B);
	vector<vector<int > >a, eps;
	vector<BigInt> BI;
	int h = fBase.size();
	vector< vector<int> > a, eps;
	for (int i = 0; i < h + 1;)
	{
		BigInt b = GerRandomBettwen(sqrt(N), N);
		BigInt A = b.KaratsubaMultiply(b) % N;
		if (is_Smooth(A, fBase, a, eps))
		{
			i++;
			BI.push_back(A);
		}
	}
	vector<int> x;
	int NumRes = gauss(eps, x);
	bool isOne = false;
	BigInt res("1");
	for (int i = 0; i < x.size(); i++)
	{
		if (x[i])
		{
			isOne = true;
			break;
		}
	}
	if (isOne)
	{
		BigInt y("1");
		for (int i = 0; i < BI.size(); i++)
		{
			if (x[i])
			{
				res = res*BI[i];
				
			}

		}
		
	}
	else
	{
	

	}
}


BigInt BigInt::Dixon(BigInt N)
{
	
}

*/
BigInt BigInt::GerRandomBettwen(BigInt L, BigInt R)
{
	srand(time(NULL));

	
	while (true)
	{
		int index = rand() % L.Number.size();
		if (L.Number[index] == 0)
		{
			L.Number[index] = 1;
			if (L < R)
				return L;
			else
			{
				L.Number[index] = 0;
			}
		}
	}
}

vector<BigInt> BigInt::PolardP_1(BigInt N)
{
	ll start = clock();
	vector<BigInt> res;
	BigInt a("2");
	BigInt One("1");
	while (!isPrime(N) && N!=BigInt("1"))
	{
		res.push_back(PolardP_1_get_factor(N, a));
		N = N / res.back();
		a = a + One;
	}
	if (isPrime(N))
		res.push_back(N);
	ll finish = clock();
	TimeOfLast = static_cast<double>(finish - start) / 1000;
	return res;
}


BigInt BigInt::PowModN(BigInt a, BigInt k, BigInt n)
{
	BigInt A = a, b = BigInt("1");
	
	if (k == 0) return b;
	if (k.Number[0] == 1) b = a;


	for (ll i = 1; i < k.Number.size(); i++) {
		A = A.KaratsubaMultiply(A) % n;
		if (k.Number[i]) b = A.KaratsubaMultiply(b) % n;
	}
	
	return b;
}


BigInt BigInt::PowN(BigInt a, BigInt N)
{
	if (N == BigInt("0"))
		return BigInt("1");
	if (N % 2 == 1)
		return PowN(a, N - 1)* a;
	else {
		BigInt b = PowN(a, N / 2);
		return b*b;
	}
}

bool BigInt::is_Smooth(BigInt A,vector<int> FBase,vector<vector<int> > &a,vector<vector<int> > &eps)
{
	vector<int> a_tmp, eps_tmp;
	a_tmp.resize(FBase.size());
	eps_tmp.resize(FBase.size());
	for (int i = 0; i<FBase.size() && A>1; i++)
	{
		int counter = 0;
		while (A%FBase[i] == 0)
		{
			A = A / FBase[i];
			counter++;
		}
		a_tmp[i] = counter;
		eps_tmp[i] = counter % 2;
	}
	if (A == 1)
	{
		a.push_back(a_tmp);
		eps.push_back(eps_tmp);
		return true;
	}
	return false;
}


bool BigInt::isPrime(BigInt N)
{
	ll start = clock();
	BigInt One("1");
	if (N == One)
		return false;
	for (BigInt i("2"); i < sqrt(N); i = i + One)
	{
		if (N%i == BigInt("0"))
			return false;
	}
	ll finish = clock();
	TimeOfLast = static_cast<double>(finish - start) / 1000;
	return true;
}


BigInt BigInt::PolardP_1_get_factor(BigInt N, BigInt a)
{
	BigInt d, k, t;
	ll K=10 ;
	k = LCM(K);
	int max_iter = 234;
	int i = 1;
	for (;; i++)
	{
		if (gcd(a, N) > 1) return gcd(a, N);
		t = PowModN(a, k, N);
		d = gcd(t - 1, N);
		if (d > 1 && d< N)
			return d;
		else
		{
			vector<int> deltas = Polard_Help_Eratosfen(K, K * K-K+2);
			BigInt q0;
			q0 = (long long)deltas.back();
			deltas.pop_back();
			BigInt c = PowModN(t, q0, N);
			for (int i = 1; i < deltas.size(); i++)
			{
				d = gcd(c-1, N);
				if (d == 1)
				{
					BigInt q;
					q = (long long)deltas[i - 1];
					c = c*PowModN(t, q, N);
				}
				else
					return d;
			}
		}
		
	}
}

vector<int> BigInt::Polard_Help_Eratosfen(int B1, int B2)
{
	vector<int> b2 = Eratosfen(B2);
	vector<int> primes;
	for (int i = 0; i < b2.size(); i++)
	{
		if (b2[i]>=B1)
			primes.push_back(b2[i]);
	}

	vector<int> res;

	int size = primes.size();
	for (int i = 0; i < size; i++)
	{
		if (i < size - 2)
		res.push_back(primes[i + 1] - primes[i]);
	}
	res.push_back(primes[0]);
	return res;
}

void BigInt::MakeFactorBase(vector<int> &B, int Base)
{
	B = Eratosfen(Base);
	int size = B.size();
	for (int i = 0; i < size; i++)
	{
		int k = B[i]*B[i];
		while (k < Base)
		{
			B.push_back(k);
			k = k*B[i];
		}
	}
}

vector<int> BigInt::Eratosfen(int B)
{
	
	vector<int> lp,res;
	lp.resize(B + 1);
	lp.assign(B+1,0);
	for (int i = 2; i <= B; ++i) 
	{
		if (lp[i] == 0) 
		{
			lp[i] = i;
			res.push_back(i);
		}
		for (int j = 0; j < (int)res.size() && res[j] <= lp[i] && i*res[j] <= B; ++j)
			lp[i * res[j]] = res[j];
	}
	return res;
}

BigInt BigInt::LCM(int K)
{
	vector<int> B;
	MakeFactorBase(B, K);
	BigInt res("1");
	for (int i = 0; i < B.size(); i++)
	{
		res = res*(long long)B[i];
	}

	return res;
}

vector<BigInt> BigInt::PollardP(BigInt N)
{
	ll start = clock();
	vector<BigInt> res;
	BigInt One("1");
	while (!isPrime(N) && N != One)
	{
		res.push_back(PollardP_get_factor(N));
		N = N / res.back();
	}
	if (isPrime(N))
		res.push_back(N);
	ll finish = clock();
	TimeOfLast = static_cast<double>(finish - start) / 1000;
	return res;
}



BigInt	BigInt::PollardP_get_factor(BigInt N)
{
	srand(time(NULL));

	ll x0_tmp = rand();
	BigInt x;
	x = x0_tmp;
	x = x%N;
	BigInt y("1");
	ll stage = 2;
	ll i = 0;
	while (gcd(N, x - y) == BigInt("1"))
	{
		if (i == stage)
		{
			y = x;
			stage *= 2;
		}
		x = (x.KaratsubaMultiply(x) + BigInt("1")) % N;
		i++;
	}
	return gcd(N, x - y);
}