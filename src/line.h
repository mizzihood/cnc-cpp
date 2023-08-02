#include <iostream>
#include <map>
#include <array>

namespace Cnc
{

typedef std::array<double, 3> TVector;
typedef std::array<TVector, 3> TMatrix;
    
class CncCorrection 
{
public:
    CncCorrection();
    void loadSettings(std::istream& is);
    void parseLine(std::string& line, std::ostream& os);
    void printSettings(std::ostream& os);
    void calculateMatrix(double angleAB, double angleAC, double angleBC);
    TVector transform(const TVector& currentPosition);
    void parseFile(std::istream& is, std::ostream& os);

private:
    std::map<std::string, double> m_settings;
    TMatrix m_matrix;
    TVector m_currentPosition;
    TVector m_currentPositionOut;

};

}