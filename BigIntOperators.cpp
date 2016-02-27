#include"BigInt.h"
#include"DefinesAndConstatns.h"
#include<math.h>
#include<cmath>
using namespace std;


// Comp operators 


bool BigInt::operator==(const BigInt& arg0) const 
{
	return !(*this < arg0) && (!(arg0 < *this));
}

bool BigInt::operator <(const BigInt& arg) const
{
	BigInt thisArg = *this;
	BigInt Arg2 = arg;
	thisArg.trim();
	Arg2.trim();
	if (thisArg.sign != Arg2.sign)
		return thisArg.sign < Arg2.sign;
	/*if (this->isZero() && arg.isZero())
	{
		return false;
	}
	if (arg.isZero())
		return (sign>-1);
	if (this->isZero())
		return arg.sign == -1;*/
	ll thisSize = thisArg.Number.size();
	ll ArgSize = Arg2.Number.size();
	if (thisSize != ArgSize)
		return thisSize*sign < ArgSize*sign;
	else
	{
		for (ll i = thisSize - 1; i >= 0; i--)
		{
			if (thisArg.Number[i] != Arg2.Number[i])
				return thisArg.Number[i]*sign < Arg2.Number[i]*sign;
		}
	}
	return false;
}

bool BigInt::operator > (const BigInt& arg0) const
{
	return arg0 < *this;
}


bool BigInt:: operator!=(const BigInt& arg0) const
{
	return !(*this == arg0);
}

bool BigInt::operator<= (const BigInt& arg0) const 
{
	return !(arg0 < *this);
}

bool BigInt::operator>= (const BigInt& arg0) const 
{
	return !(arg0 > *this);
}

bool BigInt::operator >( ll a) const
{
	BigInt a1;
	a1 = a;
	return *this > a1;
}
bool BigInt::operator >=( ll a) const
{
	BigInt a1;
	a1 = a;
	return *this >= a1;
}
bool BigInt::operator <( ll a) const
{
	BigInt a1;
	a1 = a;
	return *this < a1;
}
bool BigInt::operator <=( ll a) const
{
	BigInt a1;
	a1 = a;
	return *this <= a1;
}
bool BigInt::operator ==( ll a)const
{
	BigInt a1;
	a1 = a;
	return *this == a1;
}
bool BigInt::operator !=( ll a)const
{
	BigInt a1;
	a1 = a;
	return *this != a1;
}


//


 /// Operators BigInt OP BigInt 
BigInt BigInt::operator-()
{
	BigInt res = *this;
	res.sign = -sign;
	return res;
}
BigInt BigInt::operator+(const BigInt& arg)
{

	BigInt thisArg = *this;
	BigInt arg2 = arg;
	thisArg.trim();
	arg2.trim();
	BigInt res;
	if (thisArg.sign == arg2.sign)
	{

		res.Number.resize(__max(thisArg.Number.size(), arg2.Number.size())+1);
		int r = 0;
		for (int i = 0; i < (int)res.Number.size() || r; i++)
		{
			res.Number[i] = (i<thisArg.Number.size() ? thisArg.Number[i] : 0) +(i< arg2.Number.size() ? arg2.Number[i] :0)+r;
			r = res.Number[i] / basis;
			res.Number[i] -= r*basis;
		}
		//res.Normalize();
		res.PointPos = __max(thisArg.PointPos, arg2.PointPos);
		return res;
	}
	return *this - (-const_cast<BigInt&>(arg));
}



BigInt BigInt::operator-(const BigInt& arg)
{
	BigInt thisArg = *this;
	BigInt arg2 = arg;
	thisArg.trim();
	arg2.trim();
	BigInt res;
	if (thisArg.sign == arg2.sign)
	{
		if (thisArg.abs() >= arg2.abs())
		{
			res = thisArg;
			for (int i = 0, carry = 0; i < (int)arg2.Number.size() || carry; i++)
			{
				res.Number[i] -= carry + (i < arg2.Number.size() ? arg2.Number[i]:0);
				carry = res.Number[i] < 0;
				if (carry)
				{
					res.Number[i] += basis;
				}
			}
			res.trim();
			return res;
		}
		return -(const_cast<BigInt&>(arg) -const_cast<BigInt&>( *this));
	}
	return *this + (-const_cast<BigInt&>(arg));
}

BigInt BigInt::operator* (const BigInt& arg)
{
	
	BigInt res;

	ll size = this->Number.size() + arg.Number.size() + 10;
	res.Number.clear();
	res.Number.resize(size);
	for (int i = 0; i < (int)this->Number.size(); i++)
	{
		ll r = 0;
		for (int j = 0; j < (int)arg.Number.size(); j++)
		{
			res.Number[i + j] += this->Number[i] * arg.Number[j];
			
		}
	}
	ll r = 0;
	for (int i = 0; i<res.Number.size(); i++)
	{
		if (res.Number[i] >= basis)
		{
			r = res.Number[i] / basis;
			res.Number[i] -= basis * r;
			res.Number[i + 1] += r;
		}
	}



	//res.Normalize();
	res.trim();
	res.sign = sign*arg.sign;
	res.PointPos = this->PointPos + arg.PointPos;
	return res;
}


BigInt BigInt::operator <<(int v)
{
	if (v <= 0)
		return *this;
	BigInt res;
	trim();
	res.Number.resize(v+this->Number.size());
	int i = this->Number.size() + v;
	while (i > v)
	{
		res.Number[i - 1] = this->Number[i - v - 1];
		--i;
	}

	while (i > 0)
	{
		res.Number[i - 1] = 0;
		--i;
	}
	
	return res;
}

BigInt BigInt::operator >>(int v)
{
	if (v <= 0)
		return *this;
	if (v>Number.size())
		return BigInt("0");
	BigInt res;
	res.Number.resize(Number.size() - v);
	for (int i = 0; i < Number.size() - v; i++)
		res.Number[i] = Number[i + v];
	//res.trim();
	return res;
}

void BigInt::operator +=(const BigInt& arg)
{
	*this = *this + arg;
}

void BigInt::operator -=(const BigInt& arg)
{
	*this = *this - arg;
}
void BigInt::operator *=(const BigInt& arg)
{
	*this = *this * arg;
}

void BigInt::operator /=(const BigInt& arg)
{
	*this = *this / arg;
}
void BigInt::operator %=(const BigInt& arg)
{
	*this = *this % arg;
}


 pair<BigInt, BigInt> DivMode(const BigInt &arg1_0, const BigInt &arg2_0)
{
	int norm = arg2_0.basis / (arg2_0.Number.back() + 1);
	BigInt a = arg1_0.abs()*norm;
	BigInt b = arg2_0.abs()*norm;

	BigInt q, r;
	q.Number.resize(a.Number.size());
	r.Number.resize(a.Number.size());
	
	for (int i = a.Number.size() - 1; i >= 0; i--)
	{
		r *= a.basis;
		r += a.Number[i];

		ll s1 = r.Number.size() <= b.Number.size() ? 0 : r.Number[b.Number.size()];
		ll s2 = r.Number.size() <= b.Number.size() - 1 ? 0 : r.Number[b.Number.size() - 1];

		ll d = (a.basis * s1 + s2) / b.Number.back();
		r -= b*d;
		while (r < 0)
		{
			r += b;
			d--;
		}
		q.Number[i] = d;
	}
	q.sign = arg1_0.sign * arg2_0.sign;
	r.sign = arg1_0.sign;
	q.trim();
	r.trim();
	r = const_cast<BigInt&>(arg1_0)-(q*arg2_0);
	return make_pair(q, r/norm);
}

 BigInt BigInt::operator /(const BigInt & arg)
 {
	 return DivMode(*this, arg).first;
 }

 BigInt BigInt::operator%(const BigInt& arg)
 {
	 return DivMode(*this, arg).second;
 }

void BigInt::operator=(ll value)
{
	sign = 1;
	if (value < 0)
	{
		sign = -1;
		value = -value;
	}

	for (; value>0; value = value / basis)
	{
		Number.push_back(value%basis);
	}
}


void BigInt::operator=(const BigInt& arg)
{
	sign = arg.sign;
	Number.resize(arg.Number.size());
	//TimeOfLast = arg.TimeOfLast;
	AfterPoint = arg.AfterPoint;
	PointPos = arg.PointPos;
	for (int i = 0; i <(int) arg.Number.size(); i++)
		Number[i] = arg.Number[i];
	basis = Standart_basis;
}
BigInt BigInt::operator+(ll arg)
{
	BigInt argB;
	argB = arg;
	BigInt res;	
	res = *this + argB;
	return res;
}

BigInt BigInt::operator-(ll arg)
{
	BigInt argB;
	argB = arg;
	BigInt res;

	res = *this - argB;
	return res;

}
BigInt BigInt::operator*(ll arg)
{
	BigInt argB;
	argB = arg;
	BigInt res;

	res = *this * argB;
	return res;
}
BigInt BigInt::operator/(ll arg)
{
	BigInt res;
	res = *this;
	if (arg < 0)
		res.sign = -res.sign, arg = -arg;
	for (int i = (int)res.Number.size() - 1, rem = 0; i >= 0; --i) {
		long long cur = res.Number[i] + rem * basis;
		res.Number[i] = (int)(cur / arg);
		rem = (int)(cur % arg);
	}
	res.trim();

	return res;
}

ll BigInt::operator%(ll arg)
{
	if (arg < 0)
		arg = -arg;
	int m = 0;
	for (int i = Number.size() - 1; i >= 0; --i)
		m = (Number[i] + m * (long long)basis) % arg;
	return m ;
}

void BigInt::operator%=(ll arg)
{
	BigInt argB;
	argB = arg;	
	*this = *this % argB;
}

void BigInt::operator+=(ll arg)
{
	BigInt argB;
	argB = arg;
	*this = *this + argB;
	
}
void BigInt::operator-=(ll arg)
{
	BigInt argB;
	argB = arg;
	*this = *this - argB;

}

void BigInt::operator*=(ll arg)
{
	BigInt argB;
	argB = arg;
	*this = *this * argB;
}

void BigInt::operator/=(ll arg)
{
	BigInt argB;
	argB = arg;
	*this = *this / argB;
}




