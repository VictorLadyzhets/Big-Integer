#include"BigInt.h"
#include"DefinesAndConstatns.h"
#define _USE_MATH_DEFINES
#include<cmath>
#include<math.h>
using namespace std;

const ll BigInt::Standart_basis = 2;



// Prime numbers below 2000 using for Dixons algorithm

const int BigInt::FactorBase_2000[BigInt::LenghOfFB] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599,
601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691,
701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887,
907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997,
1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097,
1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193,
1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297,
1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399,
1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499,
1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597,
1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699,
1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789,
1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999 };

BigInt::BigInt(string Number0) : basis(Standart_basis)
{
	sign = 1;
	DischargeLenght = 0;
	ll basis0 = DecBasis;
	while (basis0 != 0)
	{

		basis0 = basis0 / 10;
		DischargeLenght++;
	}
	DischargeLenght--;
	input(Number0);
	Number = SetInBinaryBasis();
}

/*BigInt::BigInt(vll v) : basis(Standart_basis)
{
	sign = 1;
	if (v.back() < 0)
	{
		sign = -1;
		ll val = v.back();
		v.pop_back();
		v.push_back(-val);
	}

	ll r = 0;
	for (int i = 0; i < v.size(); i++)
	{
		ll a = v[i]+r;
		r = a / basis;
		a -= r*basis;
		DecNumber.push_back(a);
	}
	
	if (r)
		DecNumber.push_back(r);

	DischargeLenght = 0;
	ll basis0 = DecBasis;
	while (basis0 != 0)
	{

		basis0 = basis0 / 10;
		DischargeLenght++;
	}
	DischargeLenght--;
	

	Number = SetInBinaryBasis();
}
*/
BigInt::BigInt(vs v) : basis(Standart_basis)
{
	sign = 1;
	ll r = 0;
	for (int i = 0; i < v.size(); i++)
	{
		Number.push_back(v[i]);
	}



	DischargeLenght = 0;
	ll basis0 = DecBasis;
	while (basis0 != 0)
	{

		basis0 = basis0 / 10;
		DischargeLenght++;
	}
	DischargeLenght--;

	//Normalize();
}

BigInt BigInt::abs() const
{
	BigInt res = *this;
	res.sign *= res.sign;
	return res;
}


bool BigInt::isZero()const
{
	return Number.empty() || (Number.size() == 1 && Number[0] == 0);
}

void BigInt::PrintBin()
{
	for (int i = 0; i < Number.size(); i++)
	{
		cout << Number[i];
	}
}
 
vs BigInt::SetInBinaryBasis()
{

	vll tmp = this->DecNumber;
	vs v;
	if (tmp.size() == 1 && tmp.back() == 0)
	{
		v.push_back(0);
		return v;
	}
	/*while (tmp.size() == 1 && tmp.back() >= 2)
	{
		short r = ModLL(tmp, 2);
		tmp = DivToLL(tmp, 2);
		v.push_back(r);
	}*/
	while (tmp.size()>1 || tmp.back()>=2)
	{
		short r = ModLL(tmp, 2);
		tmp = DivToLL(tmp, 2);
		v.push_back(r);
	
	}

	
	if(!tmp.empty())
		v.push_back(1);
	return v;
}





vll BigInt::DecPlus(vll a, vll b)
{
	vll res;
	res.resize(__max(a.size(), b.size()) + 1);
	int r = 0;
	int sa = a.size(), sb = b.size();
	for (int i = 0; i < res.size() | r; i++)
	{
		res[i] = (i < sa ? a[i] : 0) + (i < sb ? b[i] : 0) + r;
		if (res[i] >= DecBasis)
		{
			res[i] -= DecBasis;
			r = 1;
		}
		else
		{
			r = 0;
		}
	}

	return res;
}
vll BigInt::DecMult(vll a, vll b)
{
	vll res;
	ll size = a.size() + b.size() + 1;
	res.clear();
	res.resize(size);
	
	for (int i = 0; i < a.size(); i++)
	{
		
		for (int j = 0; j < (int)b.size() ; j++)
		{
			res[i + j] += a[i]*b[j];
		}
	}

	ll r = 0;
	for (int i = 0; i<res.size(); i++)
	{
		if (res[i] >= DecBasis)
		{
			r = res[i] /DecBasis;
			res[i] -= DecBasis * r;
			res[i + 1] += r;
		}
	}

	while (res.size()>1 && res.back() == 0)
		res.pop_back();
	/*while (r)
	{
		res.push_back(r % 2);
		r = r / 2;
	}*/
	//res.push_back(r % 2);
	return res;
}


void BigInt::MakeDecNumber()
{
	DecNumber.clear();
	vll Two;
	Two.push_back(1);
	for (int i = 0; i < Number.size(); i++)
	{
		vll bit;
		bit.push_back(Number[i]);
		DecNumber = DecPlus(DecNumber, DecMult(bit, Two));
		bit.clear();
		bit.push_back(2);
		Two = DecMult(Two, bit);
	}
	/*if (DecNumber.size() != 1)
	{
		double two = 0.5;
		AfterPoint = 0;
		for (int i = 0; i <= PointPos; i++)
		{
			AfterPoint += Number[i] * two;
			two = two / 2;
		}
	}
	*/
}

ll BigInt::ConvertToLL(vs a)
{
	ll res = 0;
	int two = 1;
	ll size = a.size();
	for (int i = 0; i < size; i++)
	{
		res += a.back() * two;
		a.pop_back();
		two *= 2;
	}
	return res;
}

vs BigInt::ConvertToBin(ll a)
{
	vs res;
	while (a >= 2)
	{
		res.push_back(a % 2);
		a /= 2;
	}

	if (a)
	res.push_back(a);

	return res;
}
BigInt::BigInt()
{
	sign = 1;
	basis = Standart_basis;

}


BigInt::~BigInt()
{
	Number.clear();
	DecNumber.clear();
}

// Input and output
//
void BigInt::input(string Number0)
{

	sign = 1;
	DecNumber.clear();
	int pos = 0;
	while (pos < (int)Number0.size() && (Number0[pos] == '-' || Number0[pos] == '+'))
	{
		if (Number0[pos] == '-')
			sign = -sign;
		++pos;
	}
	//ll r = 0;
	for (int i = Number0.size() - 1; i >= pos; i -= DischargeLenght)
	{
		int start = i - DischargeLenght + 1;
		if (start < 0) start = 0;
		string str = Number0.substr(start, i - start + 1);
		ll num = atoi(str.c_str());

		/*r = num / basis;
		num -= r*basis;*/
		DecNumber.push_back(num);
	}
	trim();
}

void BigInt::trim() 
{
	while (Number.size() >1  && !Number.back())
		Number.pop_back();
	if (Number.empty())
		sign = 1;
}



string BigInt::GetNumber()
{
	string res;
	if (!sign)
		res += '-';
	MakeDecNumber();
	while (DecNumber.size() > 1 && DecNumber.back() == 0)
		DecNumber.pop_back();
	res += to_string(DecNumber.back());
	for (int i = DecNumber.size() - 2; i >= 0; i--)
	{
		string tmp = to_string(DecNumber[i]);
		int k = DischargeLenght - tmp.size();
		while (k--> 0)
		{
			res += "0";

		}
		res += tmp;
	}
	
	if (AfterPoint != 0)
	{
		res += ".";
		string afterPoint = to_string(AfterPoint);
			afterPoint=afterPoint.substr(2,afterPoint.size());
		res += afterPoint;
	}
	return res;
}


vll BigInt::DivToLL(vll v, ll l)
{
	vll res;
	ll  ost = 0;
	res.resize(v.size());
	for (int i = res.size() - 1; i >= 0; i--)
	{
		ll cur = ost * DecBasis + v[i];
		res[i] = cur / l;
		ost = cur % l;
	}

	while (!res.empty() && !res.back())
		res.pop_back();
	return res;
}


int BigInt::ModLL(vll v, ll l)
{
	if (l < 0)
		l = -l;
	int m = 0;
	for (int i = v.size() - 1; i >= 0; --i)
		m = (v[i] + m * (long long)basis) % l;
	return m;
}




BigInt BigInt::FastModule(ll v) const
{
	if (v <= 0 || this->isZero())
		return *this;
	ll n = this->Number.size();
	BigInt res;
	BigInt help;
	res = *this;
	BigInt number;
	number = *this;
	BigInt mod;
	mod.Number.resize(v);
	for (int i = 0; i<mod.Number.size(); i++)
	{
		mod.Number[i]=1;
	}

	while (number > mod)
	{
		res = BigInt("0");
		int size = number.Number.size();
		for (int i = 0; i < size; i += v)
		{
			help = (number >> i);
			while (help.Number.size()>v)
				help.Number.pop_back();

			res = res + help;
		}
		number = res;
	}


	if (res == mod)
		return BigInt("0");

	ll size = 1;
	while (size < v)
		size *= 2;
	while (res.Number.size()>size)
			res.Number.pop_back();
	return res;


}



ll BigInt::gcd(ll a, ll b, ll & x, ll & y) {
	if (a == 0) {
		x = 0; y = 1;
		return b;
	}
	ll x1, y1;
	ll d = gcd(b%a, a, x1, y1);
	x = y1 - (b / a) * x1;
	y = x1;
	return d;
}

ll BigInt::gcdex(ll a, ll b, ll &x, ll&y)
{
	
		if (a == 0) {
			x = 0; y = 1;
			return b;
		}
		ll x1, y1;
		ll d = gcdex(b%a, a, x1, y1);
		x = y1 - (b / a) * x1;
		y = x1;
		return d;
	
}

ll BigInt::GetRes(ll x0, ll MAX)
{

	ll x, y;
	ll g = gcdex(x0, MAX, x, y);
	
		x = (x % MAX + MAX) % MAX;
		return x;
	
}



int BigInt::FirstOne()const
{
	if (Number.size()==1 && Number[0]==0)
		return 0;
	for (int i = Number.size() - 1; i >= 0;i--)
	if (Number[i]) return i + 1;

	return 0;
}



BigInt BigInt::Cut(ll k)
{
	BigInt res=*this;
	while (res.Number.size() > k)
		res.Number.pop_back();
	return res;
}



BigInt BigInt::sqrt(BigInt N)
{
	BigInt res;
	ll pos = (N.Number.size() + 1) / 2;
	res.Number.resize(pos+1);
	BigInt One("1");
	pos--;
		while (pos >= 0)
	{
			int l = 0, r = basis;
			int curDigit = 0;
			while (l <= r)
			{
				int m = (l + r) >> 1;
				res.Number[pos] = m;
				if (res*res <= N)
				{
					curDigit = m;
					l = l + 1;
				}
				else
					r = m - 1;
			}
			res.Number[pos] = curDigit;
			pos--;
	}

		res.trim();
		return res;
}


bool BigInt::issquare(BigInt N)
{
	BigInt k;
	k = sqrt(N);
	return (k*k) == N;
}

BigInt BigInt::gcd(BigInt a, BigInt b)
{
	return b.isZero() ? a : gcd(b, a % b);
}

double BigInt::GetTime()
{
	return TimeOfLast;
}





void BigInt::fft(vc &a, bool invert)
{
	int n = (int)a.size();
	if (n == 1)  return;

	vc a0(n / 2), a1(n / 2);
	for (int i = 0, j = 0; i<n; i += 2, ++j)
	{
		a0[j] = a[i];
		a1[j] = a[i + 1];
	}
	fft(a0, invert);
	fft(a1, invert);

	double ang = 2 * M_PI / n * (invert ? -1 : 1);
	complex<double> w(1), wn(cos(ang), sin(ang));
	for (int i = 0; i<n / 2; ++i) 
	{
		a[i] = a0[i] + w * a1[i];
		a[i + n / 2] = a0[i] - w * a1[i];
		if (invert)
		{
			a[i] /= 2;
			a[i + n / 2] /= 2;
		}
		w *= wn;
	}
}


void BigInt::fft_multiply(const vector<ll> & a, const vector<ll> & b, vector<ll> & res) {
	vc fa(a.begin(), a.end()), fb(b.begin(), b.end());
	ll n = 1;
	while (n < __max(a.size(), b.size()))  n <<= 1;
	n <<= 1;
	fa.resize(n); fb.resize(n);

	fft(fa, false);
	fft(fb, false);
	for (ll i = 0; i<n; ++i)
		fa[i] *= fb[i];
	fft(fa, true);

	res.resize(n);
	for (ll i = 0; i<n; ++i)
		res[i] = (long long)(fa[i].real() + 0.5);
}


int BigInt::gauss(vector < vector<int> > a, vector<int> & ans)
{
	int n = (int)a.size();
	int m = (int)a[0].size() - 1;

	vector<int> where(m, -1);
	for (int col = 0, row = 0; col<m && row<n; ++col) {
		int sel = row;
		for (int i = row; i<n; ++i)
		if (std::abs(a[i][col]) > std::abs(a[sel][col]))
			sel = i;
		if (std::abs(a[sel][col]) < EPS)
			continue;
		for (int i = col; i <= m; ++i)
			swap(a[sel][i], a[row][i]);
		where[col] = row;

		for (int i = 0; i<n; ++i)
		if (i != row) 
		{
			double c = a[i][col] / a[row][col];
			for (int j = col; j <= m; ++j)
				a[i][j] -= a[row][j] * c;
		}
		++row;
	}

	ans.assign(m, 0);
	for (int i = 0; i<m; ++i)
	if (where[i] != -1)
		ans[i] = a[where[i]][m] / a[where[i]][i];
	for (int i = 0; i<n; ++i) {
		double sum = 0;
		for (int j = 0; j<m; ++j)
			sum += ans[j] * a[i][j];
		if (std::abs(sum - a[i][m]) > EPS)
			return 0;
	}

	for (int i = 0; i<m; ++i)
	if (where[i] == -1)
		return -1;
	return 1;
}