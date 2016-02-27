#include"BigInt.h"




BigInt BigInt::DivKuka(BigInt a, BigInt b)
{
	ll start = clock();
	if (b.isZero())
		return BigInt();
	if (b.Number.size() == 1 && b.Number.back() == 1)
		return a;
	if (b == a)
		return BigInt("1");
	if (a.isZero())
		return a;
	
	BigInt res;
	ll c = a.FirstOne() >= 1000 ? a.FirstOne() : 1000;
	int pointPos = CookInverse(b, res, c);
	
	res = res.KaratsubaMultiply(a);
	float two = 0.5;
	float afterPoint = 0;
	
	//ll k = (pointPos >= (res.Number.size() - 1)) ?  (res.Number.size() - 1) : pointPos;
	for (int i = pointPos-1; i >=0; i--)
	{
		afterPoint += two*(i<res.Number.size()?res.Number[i]:0);
		two = two / 2;
	}
	res.PointPos = pointPos;
	res = res >> pointPos;
	ll finish = clock();
	TimeOfLast = static_cast<double>(finish - start) / 1000;
	res.AfterPoint = afterPoint;
	return res;
}


int BigInt::CookInverse(BigInt v,BigInt &res,int c)
{
	ll start = clock();
	
	if (c <= 0)
	{
		res = 0;
		return 0;
	}

	int argPoint = v.FirstOne();
	res = BigInt("0");

	int resPoint = 2;
	res.Number.resize(v.Number.size()+resPoint);
	ll k = 2 * ((argPoint - 2 >= 0) ? v.Number[argPoint - 2] : 0) + ((argPoint - 3 >= 0) ? v.Number[argPoint - 3] : 0);
	if (k == 0)
	{
		res.Number[2 + resPoint - 1] = 1;
	}
	else if (k == 1)
	{
		res.Number[resPoint] = res.Number[resPoint - 1] = 1;
	}
	else if (k == 2)
	{
		res.Number[resPoint] = res.Number[resPoint - 2] = 1;
	}
	else
	{
		res.Number[resPoint] = 1;
	}
	
	ll n = 1;
	BigInt V = v;
	V.PointPos = argPoint;

	while (n < c)
	{
		V = (v >> (argPoint - 2 * n - 3)) << (argPoint - 2 * n - 3 < 0 ? 2 * n + 3 - argPoint : 0);
		V.PointPos = 2 * n + 3;
		res = (res << (resPoint + V.PointPos + 1)) - V.KaratsubaMultiply(res.KaratsubaMultiply(res));
		bool flag = false;
		for (int i = 0; !flag && i < 2 * resPoint; ++i)
		if (res.Number[i])
			flag = true;

		res = res >> (2 * resPoint);
		if (flag)
			res = res + BigInt("1");
		resPoint = V.PointPos;
		n = 2 * n;

	}

	res.PointPos = resPoint;
	res.trim();
	float two = 0.5;
	res.AfterPoint = 0;
	//ll w = (resPoint >= (res.Number.size() - 1)) ? (res.Number.size() - 1) : resPoint;
	for (int i = resPoint-1; i >= 0; i--)
	{
		res.AfterPoint += two*(i<res.Number.size() ? res.Number[i] : 0);
		two = two / 2;
	}
	ll finish = clock();
	TimeOfLast = static_cast<double>(finish - start) / 1000;
	return res.PointPos + argPoint;
}



BigInt BigInt::CookInverse(BigInt B)
{
	ll start = clock();
	BigInt res;
	res = DivKuka(BigInt("1"), B);
	ll finish = clock();
	TimeOfLast = static_cast<double>(finish - start) / 1000;
	return res;
}