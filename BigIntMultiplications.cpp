#include"DefinesAndConstatns.h"
#include"BigInt.h"

using namespace std;


BigInt BigInt::KaratsubaMultiply(const BigInt& arg) 
{
	ll start = clock();
	BigInt res;
	vs a(this->Number);
	vs b(arg.Number);
	while (a.size() > b.size())
	{
		b.push_back(0);
	}
	while (a.size() < b.size())
	{
		a.push_back(0);
	}
	while (a.size() & (a.size() - 1))
	{
		a.push_back(0);
		b.push_back(0);
	}
	vs resNum = KaratsubaMiltiplication(a, b);
	int size = resNum.size();

	for (int i = 0, carry = 0; i < (int)resNum.size(); i++) {
		long long cur = resNum[i] + carry;
		res.Number.push_back((int)(cur % basis));
		carry = (int)(cur / basis);
	}
	res.trim();
	res.sign = sign * arg.sign;
	ll finish = clock();
	TimeOfLast = static_cast<double>(finish - start) / 1000;
	return res;
}



vs BigInt::KaratsubaMiltiplication(const vs &a, const vs &b)
{

	int n = a.size();
	vs res(n + n);
	if (n <= 32)
	{
		for (int i = 0; i < a.size();i++)
		for (int j = 0; j < b.size(); j++)
			res[i + j] += a[i] * b[j];
		return res;
	}

	int k = n >> 1;
	vs a1(a.begin(), a.begin() + k);
	vs a2(a.begin() + k, a.end());
	vs b1(b.begin(), b.begin() + k);
	vs b2(b.begin() + k, b.end());

	vs a1b1 = KaratsubaMiltiplication(a1, b1);
	vs a2b2 = KaratsubaMiltiplication(a2, b2);

	for (int i = 0; i < k; i++)
		a2[i] += a1[i];
	for (int i = 0; i < k; i++)
		b2[i] += b1[i];

	vs r = KaratsubaMiltiplication(a2, b2);
	for (int i = 0; i < (int)a1b1.size(); i++)
		r[i] -= a1b1[i];
	for (int i = 0; i < (int)a2b2.size(); i++)
		r[i] -= a2b2[i];

	for (int i = 0; i < (int)r.size(); i++)
		res[i + k] += r[i];
	for (int i = 0; i < (int)a1b1.size(); i++)
		res[i] += a1b1[i];
	for (int i = 0; i < (int)a2b2.size(); i++)
		res[i + n] += a2b2[i];
	return res;
}


BigInt BigInt::ToomaKuka(BigInt a, BigInt b)
{
	
	return Toom(a, b);
}



BigInt BigInt::Toom(BigInt a, BigInt b)
{
	ll start = clock();
	ll size = __max(a.Number.size(), b.Number.size());
	const int N = 5;
	BigInt ForNewton[N];
	for (int i = 0; i < N; i++)
		ForNewton[i] = BigInt("0");

	const int parts = 3;
	ll k = 1;
	while (k < size)
	{
		k *= parts;
	}
	size = k;
	BigInt U[parts], V[parts];
	BigInt tmpA = a, tmpB = b;
	
	tmpA.Number.resize(size);
	tmpB.Number.resize(size);
	
	U[0]=tmpA >> (size / 3 * 2);
	U[0]=U[0].Cut(size / 3);
	
	U[1]=tmpA >> (size / 3);
	U[1]=U[1].Cut(size / 3);
	
	U[2]=tmpA.Cut(size / 3);

	V[0] = tmpB >> (size / 3 * 2);
	V[0] = V[0].Cut(size / 3);

	V[1] = tmpB >> (size / 3);
	V[1] = V[1].Cut(size / 3);

	V[2] = tmpB.Cut(size / 3);


	if (size / parts > 27)
	{
		for (int i = 0; i < N; i++)
		{
			BigInt a1("0"), b1("0");
			/*int multcoef = 1;
			for (int j = 0; j < parts; j++)
			{
				a1 = a1 + U[parts-j-1] * multcoef;
				b1 = b1 + V[parts-j-1] * multcoef;
				multcoef *= i;
			}
			*/
			a1 = U[0] * i*i + U[1] * i + U[2];
			b1 = V[0] * i*i + V[1] * i + V[2];
		a1.trim();
		b1.trim();
			ForNewton[i] = Toom(a1, b1);
		}
	}
	else
	{
		for (int i = 0; i < N; i++)
		{
			BigInt a1("0"), b1("0");
			/*int multcoef = 1;
			for (int j = 0; j < parts; j++)
			{
				a1 = a1 + U[parts - j - 1] * multcoef;
				b1 = b1 + V[parts - j - 1] * multcoef;
				multcoef *= i;
			}
			*/
			a1 = U[0] * i * i + U[1] * i + U[2];
			b1 = V[0] * (i * i) + V[1] * i + V[2];
			ForNewton[i] = a1.KaratsubaMultiply(b1);
			bool NeedCut = true;
			for (int j = 27; NeedCut && j < ForNewton[i].Number.size(); ++j)
			if (ForNewton[i].Number[j])
				NeedCut = false;
			if (NeedCut)
				ForNewton[i]=ForNewton[i].Cut(27);
		}
	}


	BigInt coef[N];
	for (int i = 0; i < N; i++)
		coef[i] = BigInt("0");
	for (int i = 1; i <= N; ++i)
	{
		coef[i - 1] = ForNewton[0];

		for (int j = 0; j < N - i; ++j)
		{
			ForNewton[j] = (ForNewton[j + 1] - ForNewton[j]);
			ForNewton[j] = ForNewton[j] / i;
		}
	}


	BigInt res = coef[N - 1];
	BigInt tmp = res;

	for (int i = 0; i < N; i++)
	{
		tmp = res;
		res = res << size / parts;
		for (int j = N - i - 1; j > 0; --j)
			res -= tmp;
		res += coef[N - i - 1];
	}


	res.trim();
	ll finish = clock();
	TimeOfLast = static_cast<double>(finish - start) / 1000;
	return res;
}


BigInt BigInt::Shouhage(const BigInt& a, const BigInt& b)
{		
	ll start = clock();
	
	ll NeededSize = __max(a.FirstOne(), b.FirstOne());
	ll q0 = 1;
	ll p0 = 18 * q0 + 8;
	ll sizeForRecurs = p0;
	while (p0 < NeededSize)
	{
		sizeForRecurs = p0;
		q0 = 3 * q0 - 1;
		p0 = 18 * q0 + 8;
	}

	if (p0 == 26)
		return const_cast<BigInt&>(a)*b;
	const int N = 6;

	ll m[N] = { 6 * q0 - 1, 6 * q0 + 1, 6 * q0 + 2, 6 * q0 + 3, 6 * q0 + 5, 6 * q0 + 7 };
	BigInt TwoInN[N];
	BigInt U[N];
	BigInt V[N];
	BigInt W[N];

	/*for (int i = 0; i < N; i++)
		TwoInN[i] = TwoInPower(m[i]);
		*/

	for (int i = 0; i < N; i++)
	{
		U[i] = a.FastModule(m[i]);
		while (U[i].Number.size()>sizeForRecurs)
			U[i].Number.pop_back();
		V[i] = b.FastModule(m[i]);
		while (V[i].Number.size()>sizeForRecurs)
			V[i].Number.pop_back();
	}

	for (int i = 0; i < N; i++)
	{
		W[i] = Shouhage(U[i], V[i]).FastModule(m[i]);
		while (W[i].Number.size()>sizeForRecurs)
			W[i].Number.pop_back();
	}

	int B[N - 1][N - 1]={0};
	/*
	for (int i = 0; i < N - 1;i++)
	for (int j = 0; j < N - 1; j++)
		B[i][j] = 0;*/

	for (int i = 0; i < N-1; i++)
	{
		for (int j = 0; j < i+1; j++)
		{
			B[j][i] = GetRes(m[j], m[i + 1]);
		}
	} 


	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			W[i] = HelpMult(B[j][i - 1],m[j], m[i], (W[i] < W[j] ? (BigInt("1") << m[i]) + W[i] - W[j] - BigInt("1") : W[i] - W[j]));
		}
	}

	BigInt res("0");

	for (int i = N-1; i>=0; i--)
	{
		res = (res << m[i]) - res;
		res += W[i];
	}


	res.trim();
	ll finish = clock();
	TimeOfLast = static_cast<double>(finish - start) / 1000;
	return res;
}

BigInt BigInt::HelpMult(ll b, ll e, ll f, BigInt u)
{
	ll a0 = e;
	ll d0 = (b % 2)*e;
	BigInt u0 = u;
	BigInt v0 = u*(b%2);
	b /= 2;
	while (b)
	{
	
		u0 = (u0 + (u0 << a0)).FastModule(f);
		v0 = (b % 2 ? (v0 + (u0 << d0)) : v0).FastModule(f);
		a0 = (2 * a0) % f;
	
		d0 = (d0 + (b % 2)*a0) % f;
	
		b = b/2;
	}
	
	return v0;
}


BigInt BigInt::Shtrassen(BigInt a, BigInt b)
{
	ll start = clock();
	BigInt res;
	ll two = 2;
	while (two < a.Number.size())
		two = two * 2;
	while (a.Number.size() < two)
		a.Number.push_back(0);
	
	two = 2;
	while (two < b.Number.size())
		two = two * 2;
	while (b.Number.size() < two)
		b.Number.push_back(0);
	
	
	
	
	fft_multiply(a.Number, b.Number, res.Number);
	ll carry = 0;
	for (int i = 0; i < res.Number.size(); i++)
	{
		res.Number[i] += carry;
		carry = res.Number[i] / basis;
		res.Number[i] -=carry*basis;
	}
	res.trim();
	ll finish = clock();
	TimeOfLast = static_cast<double>(finish - start) / 1000;

	return res;
}





