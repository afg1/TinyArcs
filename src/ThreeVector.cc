
#include <algorithm>
#include <iostream>
#include <cmath>

#include "ThreeVector.hh"

ThreeVector::ThreeVector()
{
    for(int i=0; i<3; i++)
    {
        elements.push_back(static_cast<long double>(0.0));
    }
}

ThreeVector::ThreeVector(long double x0, long double x1, long double x2)
{
    elements.push_back(x0);
    elements.push_back(x1);
    elements.push_back(x2);
}

ThreeVector::ThreeVector(const ThreeVector& vec)
{
    if(vec == (*this))
    {
        return;
    }
    else if (elements.size() == 3)
    {
        for(int i=0; i<3; i++)
        {
            elements[i] = vec.GetElem(i);
        }
    }
    else
    {
        for(int i=0; i<3; i++)
        {
            elements.push_back(vec.GetElem(i));
        }
    }
}

long double ThreeVector::GetElem(int i) const
{
    if(i < 3 && i >= 0)
    {
        return elements[i];
    }
    else
    {
        return 0;
    }
}

void ThreeVector::SetElem(int i, long double val)
{
    elements[i] = val;
}

ThreeVector ThreeVector::Cross(ThreeVector& vec)
{
    long double temp0(0), temp1(0), temp2(0);
    temp0 = elements[1]*vec.GetElem(2) - elements[2]*vec.GetElem(1);
    temp1 = -elements[0]*vec.GetElem(2) - elements[2]*vec.GetElem(0);
    temp2 = elements[0]*vec.GetElem(1) - elements[1]*vec.GetElem(0);
    ThreeVector retval(temp0, temp1, temp2);
    return retval;
}

long double ThreeVector::Dot(ThreeVector& vec)
{
    return (elements[0]*vec.GetElem(0) + elements[1]*vec.GetElem(1) + elements[2]*vec.GetElem(2));
}


long double ThreeVector::Magnitude()
{
    long double sum(0);
    for(int i=0; i<3; ++i)
    {
        sum += (elements[i]*elements[i]);
    }
    return sqrtl(sum);
}

ThreeVector ThreeVector::Normalise()
{
    ThreeVector rval;
    rval = *this;
    if(rval.Magnitude() != 0)
    {
        return rval /= Magnitude();
    }
    else
    {
        std::cerr << "WARNING: tried to normalise a vector with zero magnitude, returned the unaltered vector!" << std::endl;
        return rval;
    }
}

// Here be many overloads...

std::ostream& operator<<(std::ostream& os, const ThreeVector& v)
{
  os << v.elements[0] << "\t" << v.elements[1] << "\t" << v.elements[2];
  return os;
}

void ThreeVector::swap(ThreeVector& first, ThreeVector& second)
{
    using std::swap;
    swap(first.elements, second.elements);
}

ThreeVector& ThreeVector::operator=(const ThreeVector& rhs)
{
    if(rhs.elements.size() == 3 && this->elements.size() == 3)
    {
        for(int i=0; i<3; ++i)
        {
            this->elements[i] = rhs.elements[i];
        }
    }
    return *this;
}

bool ThreeVector::operator==(const ThreeVector& rhs) const
{
    if(this->elements.size() == 0 || rhs.elements.size() == 0)
    {
        return false;
    }
    for(int i=0; i<3; ++i)
    {
        if((this->elements[i] != rhs.elements[i]))
        {
            return false;
        }
    }
    return true;
}
    
ThreeVector& ThreeVector::operator+=(const ThreeVector& rhs)
{
    for(int i=0; i<3; ++i)
    {
        this->elements[i] += rhs.elements[i];
    }
    return *this;
}

ThreeVector& ThreeVector::operator-=(const ThreeVector& rhs)
{
    for(int i=0; i<3; ++i)
    {
        this->elements[i] -= rhs.elements[i];
    }
    return *this;
}

ThreeVector& ThreeVector::operator*=(const long double rhs)
{
    for(int i=0; i<3; i++)
    {
        this->elements[i] *= rhs;
    }
    return *this;
}

ThreeVector& ThreeVector::operator/=(const long double rhs)
{
    for(int i=0; i<3; ++i)
    {
        this->elements[i] /= rhs;
    }
    return (*this);
}

ThreeVector& operator+(ThreeVector lhs, const ThreeVector& rhs)
{
    lhs += rhs;
    return lhs;
}

ThreeVector& operator-(ThreeVector lhs, const ThreeVector& rhs)
{
    lhs -= rhs;
    return lhs;
}

ThreeVector& operator*(ThreeVector lhs, long double rhs)
{
    lhs *= rhs;
    return lhs;
}

ThreeVector& operator/(ThreeVector lhs, long double rhs)
{
    lhs /= rhs;
    return lhs;
}

bool ThreeVector::operator!=(const ThreeVector& rhs) const
{
    return !operator==(rhs);
}

bool ThreeVector::operator>(const ThreeVector& rhs) const
{
    bool rval(true);
    for(int i=0; i < elements.size(); ++i)
    {
        rval *= (elements[i] > rhs.GetElem(i));
    }
    return rval;
}

bool ThreeVector::operator<(const ThreeVector& rhs) const
{
    bool rval(true);
    for(int i=0; i < elements.size(); ++i)
    {
        rval *= (elements[i] < rhs.GetElem(i));
    }
    return rval;
}

bool ThreeVector::operator<=(const ThreeVector& rhs) const
{
    bool rval(true);
    for(int i=0; i < elements.size(); ++i)
    {
        rval *= (elements[i] <= rhs.GetElem(i));
    }
    return rval;
}

bool ThreeVector::operator>=(const ThreeVector& rhs) const
{
    bool rval(true);
    for(int i=0; i < elements.size(); ++i)
    {
        rval *= (elements[i] >= rhs.GetElem(i));
    }
    return rval;
}