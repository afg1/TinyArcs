#ifndef THREEVECTOR_H
#define THREEVECTOR_H 1

#include <ostream>
#include <vector>

class ThreeVector
{
    public:
        ThreeVector();
        ThreeVector(long double x0, long double x1, long double x2);
        ThreeVector(const ThreeVector&);

        ThreeVector Cross(ThreeVector&);
        long double Dot(ThreeVector&);
        
        long double GetElem(int) const;
        void SetElem(int, long double);
        void SetAll(long double val)
        {
            for(int i=0;i<3;i++)
            {
                elements[i]=val;
            }
        }
        long double Magnitude();
        ThreeVector Normalise();
    
        void BinaryWriteToStream(std::ostream &os);
    
        // Overload a whole bunch of operators!
        friend std::ostream& operator<<(std::ostream &os, const ThreeVector& v);//
    
        ThreeVector& operator=(const ThreeVector& rhs);// done?
        bool operator==(const ThreeVector& rhs) const;//done?
        bool operator!=(const ThreeVector& rhs) const;//
        bool operator>(const ThreeVector& rhs) const;//
        bool operator<(const ThreeVector& rhs) const;//
        bool operator<=(const ThreeVector& rhs) const;//
        bool operator>=(const ThreeVector& rhs) const;//All done! Note, these comparisons work on an element-by-element basis, not magnitude
    
        ThreeVector& operator+=(const ThreeVector& rhs);//
        ThreeVector& operator-=(const ThreeVector& rhs);//
        ThreeVector& operator*=(const long double rhs);//
        ThreeVector& operator/=(const long double rhs);// All done?
        friend ThreeVector& operator+(ThreeVector lhs, const ThreeVector& rhs);//
        friend ThreeVector& operator-(ThreeVector lhs, const ThreeVector& rhs);//done?
        friend ThreeVector& operator*(ThreeVector lhs, long double rhs);//
        friend ThreeVector& operator/(ThreeVector lhs,long double rhs);//
    
        void swap(ThreeVector& first, ThreeVector& second);// done?
    
    

    private:
        std::vector<long double> elements;
};
#endif

